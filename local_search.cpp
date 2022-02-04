/***************************************************************************
    local_search.cpp
    (C) 2021 by C.Blum & M.Blesa

 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

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
#include <limits>
#include <iomanip>
#include <type_traits>
#include <algorithm>
#include <chrono>
#include <random>

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

Xoroshiro mt{};

float randFloat() {
    return float(mt.next())/std::numeric_limits<uint64_t>::max();
}

int randInt(int max) {
    return mt.next() % max;
}

// Data structures for the problem data
int n_of_nodes;
int n_of_arcs;
vector< vector<int> > neighbors;

// string for keeping the name of the input file
string inputFile;

// number of applications of local search
int n_apps = 1;

// dummy parameters as examples for creating command line parameters ->
// see function read_parameters(...)
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
        if (strcmp(argv[iarg],"-i")==0) inputFile = argv[++iarg];

        // reading the number of applications of local search
        // from the command line (if provided)
        else if (strcmp(argv[iarg],"-n_apps")==0) n_apps = atoi(argv[++iarg]);

        // example for creating a command line parameter param1 ->
        // integer value is stored in dummy_integer_parameter
        else if (strcmp(argv[iarg],"-param1")==0) {
            dummy_integer_parameter = atoi(argv[++iarg]);
        }
        // example for creating a command line parameter param2 ->
        // double value is stored in dummy_double_parameter
        else if (strcmp(argv[iarg],"-param2")==0) {
            dummy_double_parameter = atof(argv[++iarg]);
        }
        iarg++;
    }
}
constexpr auto maxIterations = 100000;
constexpr auto k = 6;
constexpr auto lambda = 0.2;

struct State {
    int total;
    Solution sol;
    Hs hs;
};

int evalState(Graph const&g, int diff, State const& state) {
    int res = state.total;

    if (diff == -1) {
    } else if (state.sol[diff]) {
        --res;
    } else {
        ++res;
    }

    return res;
}

void makeNextState(Graph const&g, int &node, State const& current) {
    for (bool success = false; !success;) {
        success = true;
        node = randInt(current.sol.size()-1);
        //std::cout << int(current.sol[node]) << " " << g[node].size() << std::endl;
        if (current.sol[node]) {
            for (auto u : g[node]) {
                if (current.hs[u] == 0) {
                    success = false;
                    break;
                }
            }
        }
    }
}

float pValue(int current, int next, float temperature) {
    if (current > next)
        return 1;
    return (current/float(next))*temperature;
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

/*
inline int is_mpids_assisted(Graph const& neighbors, State& sub) {
    //ara falta comprovar si és minimal!
    for(int i = 0; i < sub.hs.size(); ++i) { // |V|
        if (sub.sol[i]) {
            sub.sol[i] = false;
            for (auto && u : neighbors[i])
                ++sub.hs[u];
            bool is_sub_pids = is_pids_assisted(neighbors, sub.hs);
            sub.sol[i] = true;
            for (auto && u : neighbors[i])
                --sub.hs[u];
            if(is_sub_pids)
                return i; //O(|V| + |E|)
        }
    }
    return -1;
}
*/
State getInitialState(Graph const& g) {
    State stateCurrent;
    auto res = greedy_genmpid_impl(g);
    stateCurrent.sol = std::move(res.first);
    stateCurrent.hs = std::move(res.second);
    //stateCurrent.sol = Solution(g.size(), 1);
    //stateCurrent.hs = compute_hs(g, stateCurrent.sol);
    stateCurrent.total = std::accumulate(stateCurrent.sol.begin(), stateCurrent.sol.end(), 0);
    return stateCurrent;
}

__attribute__((noinline)) State simulAneal(Graph const& g) {
    auto tempLoss = 1/float(maxIterations);
    State stateCurrent = getInitialState(g);
    int diffNode = -1;
    int valueCurrent = evalState(g, diffNode, stateCurrent);
    //std::cout << printSol(stateCurrent.sol) << " " << valueCurrent << endl;

    int valueNext;
    for (int i = 0; i < maxIterations; ++i) {
        float tempPercent = (exp(-k*i*tempLoss)*lambda);
        makeNextState(g, diffNode, stateCurrent);
        valueNext = evalState(g, diffNode, stateCurrent);

        if (pValue(valueCurrent, valueNext, tempPercent) >= randFloat()) {
            //std::cout << "vc: " << valueCurrent << " vn: " << valueNext << " p: " << pValue(valueCurrent, valueNext, tempPercent) << " sol: " << diffNode << endl;
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
            valueCurrent = valueNext;
        }
    }

    int idx;
    while((idx = is_mpids_assisted(g, stateCurrent)) >= 0) {
        stateCurrent.sol[idx] = false;
        for (auto u : g[idx]) {
            ++stateCurrent.hs[u];
        }
        --stateCurrent.total;
    }

    return stateCurrent;
}

/**********
Main function
**********/

int main( int argc, char **argv ) {

    read_parameters(argc,argv);

    // setting the output format for doubles to 2 decimals after the comma
    std::cout << std::setprecision(2) << std::fixed;

    // vectors for storing the result and the computation time
    // obtained by the <n_apps> applications of local search
    vector<double> results(n_apps, std::numeric_limits<int>::max());
    vector<double> times(n_apps, 0.0);

    // opening the corresponding input file and reading the problem data
    ifstream indata;
    indata.open(inputFile.c_str());
    if(not indata) { // file couldn't be opened
        cout << "Error: file could not be opened" << endl;
        return 0;
    }

    indata >> n_of_nodes;
    indata >> n_of_arcs;
    neighbors = vector< vector<int> >(n_of_nodes);
    int u, v;
    while(indata >> u >> v) {
        neighbors[u-1].push_back(v-1);
        neighbors[v-1].push_back(u-1);
    }
    indata.close();

    // main loop over all applications of local search
    for (int na = 0; na < n_apps; ++na) {

        cout << "start application " << na + 1 << endl;
        Timer timer;
        State res = simulAneal(neighbors);
        double st = timer.elapsed_time(Timer::VIRTUAL);

        int total = std::accumulate(res.sol.begin(), res.sol.end(), 0);
        cout << "data: " << is_mpids(neighbors, res.sol) << " " << total << endl;

        // The starting solution for local search may be randomly generated,
        // or you may incorporate your greedy heuristic in order to produce
        // the starting solution.

        // Whenever you move to a new solution, first take the computation
        // time as explained above. Say you store it in variable ct.
        // Then, write the following to the screen:
        cout << "value " << total;
        cout << "\ttime " << st << endl;

        // When a local minimum is reached, store the value of the
        // corresponding solution in vector results:
        results[na] = evalState(neighbors, -1, res);

        // Finally store the needed computation time (that is, the time
        // measured once the local minimum is reached) in vector times:
        times[na] = st;

        cout << "end application " << na + 1 << endl;
    }

    // calculating the average of the results and computation times, and
    // their standard deviations, and write them to the screen
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
    cout << r_best << "\t" << r_mean << "\t" << rsd << "\t";
    cout << t_mean << "\t" << tsd << endl;
}

