#pragma once

#include <algorithm>
#include <string>
#include <vector>
#include <sstream>

namespace DetailImpl {

using Graph = std::vector<std::vector<int>>;
using Solution = std::vector<uint8_t>;
using Hs = std::vector<int32_t>;

template<class T, class U>
inline auto ceildiv(T i, U j) {
    return i/j + (i % j != 0);
}

//Comprovem si sub es PIDS (NO MINIMAL) de neighbors (Funció auxiliar)
//COST: O(|V| + |E|) (recorrem vertex i arestes 1 vegada)
inline bool is_pids(Graph const& neighbors, Solution& sub) {
//Recorrem tots els vèrtex, i per cada adjacència comprovem si pertany a sub.
    //Si menys de la meitat de les adj pertanyen a sub per algun node, no es PIDS
    for(int i = 0; i < neighbors.size(); ++i) {
        auto & node = neighbors[i];
        int adj_in_set = 0;
        for(auto adj : node) {
            if(sub[adj]) ++adj_in_set;
            if(adj_in_set >=  ceildiv(node.size(),2)) break; //Per estalviar
                                                    // iteracions innecessàries
        }
        if(adj_in_set < ceildiv(node.size(),2)) {
            //cout << i << endl;
            return false;
        }
    }
    return true;
}

inline bool is_pids_assisted(Graph const& neighbors, Hs const& hs) {
    for(int i = 0; i < neighbors.size(); ++i) {
        if (hs[i] > 0) return false;
    }
    return true;
}

inline Hs compute_hs(Graph const& neighbors, Solution& sub) {
    Hs hs(neighbors.size());
    for (int i = 0; i < neighbors.size(); ++i) {
        hs[i] = ceildiv(neighbors[i].size(),2);
        for (auto && u : neighbors[i])
            hs[i] -= sub[u];
    }
    return hs;
}

//Comprovem si sub es MPIDS de neighbors, retorna ()
// COST: O(|V|^2 + |V|*|E|)
inline bool is_mpids(Graph const& neighbors, Solution& sub) {
    Hs hs = compute_hs(neighbors, sub);
    bool pids = is_pids_assisted(neighbors, hs); //O(|V| + |E|)
    //ara falta comprovar si és minimal!
    if(pids) {
        for(int i = 0; i < sub.size(); ++i) { // |V|
            if (sub[i]) {
                if (std::all_of(neighbors[i].begin(), neighbors[i].end(), [&](int val){return hs[val] < 0;})) {
                    return i;
                }
            }
        }
    } else {
        return false;
    }
    return true;
}

template<class T>
inline std::string printSol(T const& sol) {
    std::stringstream s;
    for (auto i : sol) {
        s << int(i);
    }
    return s.str();
}

}