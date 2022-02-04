/***************************************************************************
    greedy.cpp
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

using namespace DetailImpl;
using namespace GreedyImpl;

// global variables concerning the random number generator (in case needed)
time_t t;
Random* rnd;

// Data structures for the problem data
int n_of_nodes;
int n_of_arcs;
vector< vector<int> > neighbors;

// string for keeping the name of the input file
string inputFile;

// dummy parameters as examples for creating command line parameters
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

        // example for creating a command line parameter param1
        //-> integer value is stored in dummy_integer_parameter
        else if (strcmp(argv[iarg],"-param1")==0)
            dummy_integer_parameter = atoi(argv[++iarg]);

        // example for creating a command line parameter param2
        //-> double value is stored in dummy_double_parameter
        else if (strcmp(argv[iarg],"-param2")==0)
            dummy_double_parameter = atof(argv[++iarg]);

        iarg++;
    }
}

/************
Main function
*************/

int main( int argc, char **argv ) {

    read_parameters(argc,argv);

    // setting the output format for doubles to 2 decimals after the comma
    std::cout << std::setprecision(2) << std::fixed;

    // initializing the random number generator.
    // A random number in (0,1) is obtained with: double rnum = rnd->next();
    rnd = new Random((unsigned) time(&t));
    rnd->next();

    // variables for storing the result and the computation time
    // obtained by the greedy heuristic
    double results = std::numeric_limits<int>::max();
    double time = 0.0;

    // opening the corresponding input file and reading the problem data
    ifstream indata;
    indata.open(inputFile.c_str());
    if(not indata) { // file couldn't be opened
        cout << "Error: file could not be opened" << endl;
    }

    indata >> n_of_nodes;
    indata >> n_of_arcs;
    cout << n_of_nodes << " nodes" << endl;
    cout << n_of_arcs << " arcs" << endl;
    neighbors = Graph(n_of_nodes);
    int u, v;
    while(indata >> u >> v) {
        neighbors[u-1].push_back(v-1);
        neighbors[v-1].push_back(u-1);
    }
    indata.close();
    /*for(auto& neighbour : neighbors) {
        cout << "NEIGHBOR " << endl;
        for(int x : neighbour) {
            cout <<  x << " ";
        }
        cout << endl;
    }*/
    // the computation time starts now
    {

    // Example for requesting the elapsed computation time at any moment:
    // double ct = timer.elapsed_time(Timer::VIRTUAL);
    auto t_start = std::chrono::high_resolution_clock::now();
    // at greedy.hpp
    auto res = greedy_genmpid(neighbors);
    auto t_end1 = std::chrono::high_resolution_clock::now();
    // at base.hpp
    auto valid = is_mpids(neighbors, res);
    auto t_end2 = std::chrono::high_resolution_clock::now();
    int total = std::accumulate(res.begin(), res.end(), 0);
    cout << "RES: " << total << endl;
    double elapsed_time_milli = std::chrono::duration<double, std::milli>(t_end1-t_start).count();
    double elapsed_time_milli2 = std::chrono::duration<double, std::milli>(t_end2-t_end1).count();

    cout << elapsed_time_milli << "ms is " << (valid ? "true" : "false") << " validation: " << elapsed_time_milli2 << "ms" << endl;
    }

}

