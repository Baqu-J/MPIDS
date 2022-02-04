/***************************************************************************
    metaheuristic.cpp
    (C) 2021 by C. Blum & M. Blesa

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <numeric>
#include <random>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "Timer.h"
#include "Random.h"
#include "basics.h"
#include "greedy.hpp"

#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <vector>
#include <set>
#include <limits>
#include <iomanip>
#include <deque>
#include <map>
#include <thread>

using namespace DetailImpl;
using namespace GreedyImpl;
using namespace std;
// global variables concerning the random number generator (in case needed)
time_t t;
class Xoroshiro {
    uint64_t shuffle_table[4];
public:
    Xoroshiro() {
        std::random_device rd;
        for (auto &val : shuffle_table) {
            val = rd();
        }
    }
    uint64_t next()
    {
        uint64_t s1 = shuffle_table[0];
        uint64_t s0 = shuffle_table[1];
        uint64_t result = s0 + s1;
        shuffle_table[0] = s0;
        s1 ^= s1 << 23;
        shuffle_table[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
        return result;
    }
};

float randFloat() {
    static thread_local Xoroshiro mt{};
    return float(mt.next())/std::numeric_limits<uint64_t>::max();
}

auto timeSeconds() {
    thread_local static auto start = std::chrono::system_clock::now();
    return std::chrono::duration<double, std::ratio<1,1>>(std::chrono::system_clock::now()-start).count();
}

// Data structures for the problem data
int n_of_nodes;
int n_of_arcs;
Graph neighbors;

// string for keeping the name of the input file
string inputFile;

// computing time limit for each application of the metaheuristic
double time_limit = 10.0;

// number of applications of the metaheuristic
int n_apps = 1;

// dummy parameters as examples for creating command line parameters
// (see function read_parameters(...))
int dummy_integer_parameter = 0;
int dummy_double_parameter = 0.0;


inline int stoi(string &s) {

  return atoi(s.c_str());
}

inline double stof(string &s) {

  return atof(s.c_str());
}

void read_parameters(int argc, char **argv) {

    int iarg = 1;
    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i") == 0) inputFile = argv[++iarg];
        // reading the computation time limit
        // from the command line (if provided)
        else if (strcmp(argv[iarg],"-t") == 0) time_limit = atoi(argv[++iarg]);
        // reading the number of applications of the metaheuristic
        // from the command line (if provided)
        else if (strcmp(argv[iarg],"-n_apps") == 0) n_apps = atoi(argv[++iarg]);
        // example for creating a command line parameter
        // param1 -> integer value is stored in dummy_integer_parameter
        else if (strcmp(argv[iarg],"-param1") == 0) {
            dummy_integer_parameter = atoi(argv[++iarg]);
        }
        // example for creating a command line parameter
        // param2 -> double value is stored in dummy_double_parameter
        else if (strcmp(argv[iarg],"-param2") == 0) {
            dummy_double_parameter = atof(argv[++iarg]);
        }
        iarg++;
    }
}

struct State {
    int total;
    Solution sol;
    Hs hs;
    std::unordered_map<uint32_t, uint32_t> bannedStates;
};

int evalState(Graph const&g, int64_t diff, State const& state) {
    int res = state.total;

    if (diff == -1) {
    } else if (state.sol[diff]) {
        --res;
    } else {
        ++res;
    }
    return res;
}

inline int64_t computeMaxItWithoutImprovement(int64_t i) {
    return 500 * (i/(1024*32)+1);
}

void makeNextState(Graph const&g, int &node, State const& current, int64_t i, int32_t iterationsWithoutImprovement, Solution const&mustBeTrue) {
    int score = std::numeric_limits<int>::max();
    auto randMin = 0.5;
    for (uint32_t i = 0; i < g.size(); ++i) {
        if(mustBeTrue[i] || current.bannedStates.find(i) != current.bannedStates.end())
            continue;
        if (!current.sol[i] || std::none_of(g[i].begin(), g[i].end(), [&](auto u){return current.hs[u] == 0;})) {
            if (int candEval = !current.sol[i]; score > candEval) {
                score = candEval;
                node = i;
            } else if (score == candEval) {
                if (!candEval) {
                    if (g[i].size() < g[node].size() ||
                        randFloat() >= randMin)
                        node = i;
                } else {
                    if (g[i].size() > g[node].size() ||
                        randFloat() >= randMin)
                        node = i;
                }
            }
        }
    }
    if (score == std::numeric_limits<int>::max()) {
        for (auto [i, _] : current.bannedStates) {
            if(mustBeTrue[i])
                continue;
            if (!current.sol[i] || std::none_of(g[i].begin(), g[i].end(), [&](auto u){return current.hs[u] == 0;})) {
                if (int candEval = !current.sol[i]; score > candEval) {
                    score = candEval;
                    node = i;
                } else if (score == candEval) {
                    if (!candEval) {
                        if (g[i].size() < g[node].size())
                            node = i;
                    } else {
                        if (g[i].size() > g[node].size())
                            node = i;
                    }
                }
            }
        }
    }
}

inline int is_mpids_assisted(Graph const& neighbors, State& sub) {
    //ara falta comprovar si és minimal!
    for(int i = 0; i < sub.hs.size(); ++i) { // |V|
        if (sub.sol[i]) {
            if (std::all_of(neighbors[i].begin(), neighbors[i].end(), [&](int val){return sub.hs[val] < 0;})) {
                return i;
            }
        }
    }
    return -1;
}

State getInitialState(Graph const& g) {
    State stateCurrent;
    auto res = greedy_genmpid_impl(g);
    stateCurrent.sol = std::move(res.first);
    stateCurrent.hs = std::move(res.second);

    stateCurrent.total = std::accumulate(stateCurrent.sol.begin(), stateCurrent.sol.end(), 0);
    return stateCurrent;
}

void applyToState(Graph const&g, int &diffNode, State& stateCurrent) {
    if (stateCurrent.sol[diffNode]) {
        --stateCurrent.total;
        stateCurrent.sol[diffNode] = false;
        for (auto u : g[diffNode]) {
            ++stateCurrent.hs[u];
        }
    } else {
        ++stateCurrent.total;
        stateCurrent.sol[diffNode] = true;
        for (auto u : g[diffNode]) {
            --stateCurrent.hs[u];
        }
    }
}

inline void reduce(State& s, Graph const& g) {
    for (int v = 0; v < s.sol.size(); ++v) {
        if (!s.sol[v]) continue;

        if (std::all_of(g[v].begin(), g[v].end(), [&](auto u){
            return s.hs[u] < 0;
        })) {
            s.sol[v] = false;
            --s.total;
            for (auto i : g[v]) {
                ++s.hs[i];
            }
        }
    }
}

__attribute__((noinline)) State modifiedTabu(Graph const& g) {
    State stateCurrent = getInitialState(g);
    auto mustBeTrue = graphPrunning(g);
    int diffNode = -1;
    int valueCurrent = evalState(g, diffNode, stateCurrent);

    State bestCandidate = stateCurrent;
    int valueBest = valueCurrent;

    int valueNext;
    int iterationsWithoutImprovement = 0;
    for (int64_t i = 0; timeSeconds() < time_limit; ++i) {
        makeNextState(g, diffNode, stateCurrent, i, iterationsWithoutImprovement, mustBeTrue);
        valueNext = evalState(g, diffNode, stateCurrent);

        if (valueNext < valueCurrent) {
            applyToState(g, diffNode, stateCurrent);
            for (auto it = stateCurrent.bannedStates.begin(); it != stateCurrent.bannedStates.end();)
                (--(it->second) == 0) ? it = stateCurrent.bannedStates.erase(it) : ++it;
        } else {
            if (valueCurrent <= valueBest) {
                bestCandidate.sol = stateCurrent.sol;
                bestCandidate.hs = stateCurrent.hs;
                bestCandidate.total = stateCurrent.total;
                valueBest = valueCurrent;
                iterationsWithoutImprovement = -1;
            }
            if (iterationsWithoutImprovement >= computeMaxItWithoutImprovement(i)) {
                stateCurrent = bestCandidate;
                iterationsWithoutImprovement = -1;
            } else {
                applyToState(g, diffNode, stateCurrent);
                stateCurrent.bannedStates.emplace(diffNode, 1 * (i/(1024*32)+20));
            }
        }
        iterationsWithoutImprovement++;
        valueCurrent = valueNext;
    }
    if (valueCurrent < valueBest) {
        bestCandidate.sol = stateCurrent.sol;
        bestCandidate.hs = stateCurrent.hs;
        bestCandidate.total = stateCurrent.total;
        valueBest = valueCurrent;
    }

    reduce(bestCandidate, g);

    return bestCandidate;
}

/**********
Main function
**********/

int main( int argc, char **argv ) {

    read_parameters(argc,argv);

    // setting the output format for doubles to 2 decimals after the comma
    std::cout << std::setprecision(2) << std::fixed;

    // vectors for storing the result and the computation time
    // obtained by the <n_apps> applications of the metaheuristic
    vector<double> results(n_apps, std::numeric_limits<int>::max());
    vector<double> times(n_apps, 0.0);

    // opening the corresponding input file and reading the problem data
    ifstream indata;
    indata.open(inputFile.c_str());
    if (not indata) { // file couldn't be opened
        cout << "Error: file could not be opened" << endl;
    }

    indata >> n_of_nodes;
    indata >> n_of_arcs;
    neighbors = Graph(n_of_nodes);
    int u, v;
    while(indata >> u >> v) {
        neighbors[u-1].push_back(v-1);
        neighbors[v-1].push_back(u-1);
    }
    indata.close();

    // We don't want to overstress our CPU, so we'll at most run half the available threads at the same time.
    for (int rep = 0; rep < n_apps; rep+=std::thread::hardware_concurrency()/2) {
        auto runTotal = std::min(n_apps - rep,int(std::thread::hardware_concurrency()/2));
        std::vector<std::thread> parallelApps;
        // main loop over all applications of the metaheuristic
        for (int na = 0; na < runTotal; ++na) {
            parallelApps.push_back(std::thread([na, rep, &results, &times]{
                //cout << "start application " << na + 1 << endl;
                timeSeconds();
                State res = modifiedTabu(neighbors);

                int total = std::accumulate(res.sol.begin(), res.sol.end(), 0);
                std::stringstream output;
                output << "data: " << is_mpids(neighbors, res.sol) << " " << total << endl;
                // HERE GOES YOUR METAHEURISTIC

                // For implementing the metaheuristic you probably want to take profit
                // from the greedy heuristic and/or the local search method that you
                // already developed.
                //
                // Whenever the best found solution is improved, first take the
                // computation time as explained above. Say you store it in variable ct.
                // Then, write the following to the screen:
                output << "value " << total;
                output << "\ttime " << timeSeconds() << endl;
                //
                // Store the value of the new best found solution in vector results:
                results[na] = total;
                times[na] = timeSeconds();
                //
                // Stop the execution of the metaheuristic
                // once the time limit "time_limit" is reached.

                output << "end application " << na + 1 << endl;
                cout << output.str();
            }));
        }

        for (auto &t : parallelApps)
            t.join();
    }
    // calculating the average of the results and computation times,
    // and their standard deviations, and write them to the screen
    double r_mean = 0.0;
    int r_best = std::numeric_limits<int>::max();
    double t_mean = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        r_mean = r_mean + results[i];
        if (int(results[i]) < r_best) r_best = int(results[i]);
        t_mean = t_mean + times[i];
    }
    r_mean = r_mean/double(results.size());
    t_mean = t_mean/double(times.size());
    double rsd = 0.0;
    double tsd = 0.0;
    for (int i = 0; i < int(results.size()); i++) {
        rsd = rsd + pow(results[i]-r_mean,2.0);
        tsd = tsd + pow(times[i]-t_mean,2.0);
    }
    rsd = rsd/double(results.size());
    if (rsd > 0.0) {
        rsd = sqrt(rsd);
    }
    tsd = tsd/double(results.size());
    if (tsd > 0.0) {
        tsd = sqrt(tsd);
    }
    // printing statistical information
    cout << r_best << "\t" << r_mean << "\t" << rsd << "\t";
    cout << t_mean << "\t" << tsd << endl;
}

