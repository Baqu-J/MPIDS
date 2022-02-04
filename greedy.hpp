#pragma once

#include <algorithm>
#include <cstring>
#include <numeric>
#include <string>
#include <vector>
#include <sstream>
#include "basics.h"

namespace GreedyImpl {
using namespace DetailImpl;

inline auto nsSize(int v, Graph const& g, Solution const& s) {
    uint32_t res = 0;
    for (auto && u : g[v])
        res += s[u];
    return res;
}

inline auto hsCompute(int v, Graph const& g, Solution const& s) {
    return int64_t(ceildiv(g[v].size(), 2)) - int64_t(nsSize(v, g, s));
}

// O(|V|)
inline Solution graphPrunning(Graph const& graph) {
    //O(|V|)
    Solution s(graph.size(), false);

    // O(|V|)
    for (int v = 0; v < graph.size(); ++v) {
        if (graph[v].size() != 1)
            continue;
        int u = *graph[v].begin();

        if (!s[u]) {
            s[u] = true;
        }

        if (graph[u].size() == 2 && !s[v]) {
            int w = [&]{
                auto itu = graph[u].begin();
                if (*itu == v) ++itu;
                return *itu;
            }(); // w != v, segon vei de u

            if (!s[w]) s[w] = true;
        }
    }

    return s;
}

inline uint64_t cover_degree(int v, Graph const& g, Solution const& s, Hs const& hs) {
    uint64_t res = 0;
    for (int u : g[v]) {
        if (hs[u] > 0)
            ++res;
    }
    return res;
}

inline uint64_t need_degree(int v, Graph const& g, Solution const& s, Hs const& hs) {
    uint64_t res = 0;
    for (int u : g[v]) {
        res += std::max(hs[u], int32_t(0));
    }
    return res;
}

inline void reduce(Solution& s, Graph const& g, Hs &hs) {
    for (int v = 0; v < s.size(); ++v) {
        if (!s[v]) continue;

        if (std::all_of(g[v].begin(), g[v].end(), [&](auto u){
            return hs[u] < 0;
        })) {
            s[v] = false;
            for (auto i : g[v]) {
                ++hs[i];
            }
        }
    }
}

inline void implOfficial(Solution &s, Hs &hs, std::vector<uint32_t> &c, std::vector<uint32_t> &cover_degrees, Graph const& graph) {
    std::sort(c.begin(), c.end(), [&graph](uint32_t a,uint32_t b){
        return graph[a].size() < graph[b].size();
    });

    auto viIt = c.begin();
    while (!is_pids_assisted(graph, hs)) {
        //cout << "IT: " << printSol(s) << endl;
        int vi, p;
        while (viIt != c.end()) {
            if ((p = hs[*viIt]) > 0) {
                vi = *viIt;
                break;
            }
            ++viIt;
        }

        //cout << "vi: " << vi << endl;
        //cout << "p: " << p << endl;

        for (int j = 0; j < p; ++j) {
            std::memset(&cover_degrees[0], 0xFF, graph[vi].size()*sizeof(uint32_t));
            {
                int i = 0;
                for (int u : graph[vi]) {
                    if (!s[u]) {
                        //cout << "u cover_degrees: " << u << endl;
                        cover_degrees[i] = cover_degree(u, graph, s, hs);
                    }
                    ++i;
                }
            };
            /*cout << "cover_degrees: ";
            for (auto cd : cover_degrees) {
                cout << cd << ",";
            }
            //cout << endl;*/

            uint32_t cdmax = [&]{
                uint32_t cdmax = 0;
                for (int i = 0; i < graph[vi].size(); ++i) {
                    if (cover_degrees[i] != -1)
                        cdmax = std::max(cdmax, cover_degrees[i]);
                }
                return cdmax;
            }();

            //cout << "cdmax: " << cdmax << endl;
            uint32_t u = [&] {
                uint32_t umax = 0;
                uint32_t needumax = 0;

                auto uit = graph[vi].begin();
                int i = 0;
                while (uit != graph[vi].end()) {
                    if (cover_degrees[i] != -1 && cover_degrees[i] == cdmax) {
                        auto needu = need_degree(*uit, graph, s, hs);
                        //cout << "uit: " << *uit << " needu " << needu << endl;
                        if (needumax < needu){
                            needumax = needu;
                            umax = *uit;
                        }
                    }

                    ++i;
                    ++uit;
                }
                return umax;
            }();
            //cout << "u: " << u << endl;

            if (!s[u]) {
                s[u] = true;
                for (auto i : graph[u]) {
                    --hs[i];
                }
            }

        }
        ++viIt;
    }
    reduce(s, graph, hs);
}

inline void implVariation1(Solution &s, Hs &hs, std::vector<uint32_t> &c, std::vector<uint32_t> &cover_degrees, Graph const& graph) {
    std::sort(c.begin(), c.end(), [&graph](uint32_t a,uint32_t b){
        return graph[a].size() < graph[b].size();
    });

    auto viIt = c.begin();
    while (!is_pids_assisted(graph, hs)) {
        //cout << "IT: " << printSol(s) << endl;
        int vi, p;
        while (viIt != c.end()) {
            if ((p = hs[*viIt]) > 0) {
                vi = *viIt;
                break;
            }
            ++viIt;
        }

        //cout << "vi: " << vi << endl;
        //cout << "p: " << p << endl;

        for (int j = 0; j < p; ++j) {
            std::memset(&cover_degrees[0], 0xFF, graph[vi].size()*sizeof(uint32_t));
            {
                int i = 0;
                for (int u : graph[vi]) {
                    if (!s[u]) {
                        //cout << "u cover_degrees: " << u << endl;
                        cover_degrees[i] = cover_degree(u, graph, s, hs);
                    }
                    ++i;
                }
            };
            /*cout << "cover_degrees: ";
            for (auto cd : cover_degrees) {
                cout << cd << ",";
            }
            //cout << endl;*/

            uint32_t cdmax = [&]{
                uint32_t cdmax = 0;
                for (int i = 0; i < graph[vi].size(); ++i) {
                    if (cover_degrees[i] != -1)
                        cdmax = std::max(cdmax, cover_degrees[i]);
                }
                return cdmax;
            }();

            //cout << "cdmax: " << cdmax << endl;
            uint32_t u = [&] {
                uint32_t umax = 0;
                uint32_t needumax = 0;

                auto uit = graph[vi].begin();
                int i = 0;
                while (uit != graph[vi].end()) {
                    if (cover_degrees[i] != -1 && cover_degrees[i] == cdmax) {
                        auto needu = need_degree(*uit, graph, s, hs);
                        //cout << "uit: " << *uit << " needu " << needu << endl;
                        if (needumax < needu || (needumax == needu && graph[*uit].size() > graph[umax].size())){
                            needumax = needu;
                            umax = *uit;
                        }
                    }

                    ++i;
                    ++uit;
                }
                return umax;
            }();
            //cout << "u: " << u << endl;

            if (!s[u]) {
                s[u] = true;
                for (auto i : graph[u]) {
                    --hs[i];
                }
            }

        }
        ++viIt;
    }
    reduce(s, graph, hs);
}

inline void implVariation2(Solution &s, Hs &hs, std::vector<uint32_t> &c, std::vector<uint32_t> &cover_degrees, Graph const& graph) {
    std::sort(c.begin(), c.end(), [&hs](uint32_t a,uint32_t b){
        return hs[a] > hs[b];
    });

    auto viIt = c.begin();
    while (!is_pids_assisted(graph, hs)) {
        //cout << "IT: " << printSol(s) << endl;
        int vi, p;
        while (viIt != c.end()) {
            if ((p = hs[*viIt]) > 0) {
                vi = *viIt;
                break;
            }
            ++viIt;
        }

        //cout << "vi: " << vi << endl;
        //cout << "p: " << p << endl;

        for (int j = 0; j < p; ++j) {
            std::memset(&cover_degrees[0], 0xFF, graph[vi].size()*sizeof(uint32_t));
            {
                int i = 0;
                for (int u : graph[vi]) {
                    if (!s[u]) {
                        //cout << "u cover_degrees: " << u << endl;
                        cover_degrees[i] = cover_degree(u, graph, s, hs);
                    }
                    ++i;
                }
            };
            /*cout << "cover_degrees: ";
            for (auto cd : cover_degrees) {
                cout << cd << ",";
            }
            //cout << endl;*/

            uint32_t cdmax = [&]{
                uint32_t cdmax = 0;
                for (int i = 0; i < graph[vi].size(); ++i) {
                    if (cover_degrees[i] != -1)
                        cdmax = std::max(cdmax, cover_degrees[i]);
                }
                return cdmax;
            }();

            //cout << "cdmax: " << cdmax << endl;
            uint32_t u = [&] {
                uint32_t umax = 0;
                uint32_t needumax = 0;

                auto uit = graph[vi].begin();
                int i = 0;
                while (uit != graph[vi].end()) {
                    if (cover_degrees[i] != -1 && cover_degrees[i] == cdmax) {
                        auto needu = need_degree(*uit, graph, s, hs);
                        //cout << "uit: " << *uit << " needu " << needu << endl;
                        if (needumax < needu || (needumax == needu && graph[*uit].size() > graph[umax].size())){
                            needumax = needu;
                            umax = *uit;
                        }
                    }

                    ++i;
                    ++uit;
                }
                return umax;
            }();
            //cout << "u: " << u << endl;

            if (!s[u]) {
                s[u] = true;
                for (auto i : graph[u]) {
                    --hs[i];
                }
            }

        }
        ++viIt;
    }
    reduce(s, graph, hs);
}

inline void implVariation3(Solution &s, Hs &hs, std::vector<uint32_t> &c, std::vector<uint32_t> &cover_degrees, Graph const& graph) {
    std::sort(c.begin(), c.end(), [&hs, &graph](uint32_t a,uint32_t b){
        return hs[a] > hs[b] || (hs[a] == hs[b] && graph[a].size() < graph[b].size());
    });

    auto viIt = c.begin();
    while (!is_pids_assisted(graph, hs)) {
        //cout << "IT: " << printSol(s) << endl;
        int vi, p;
        while (viIt != c.end()) {
            if ((p = hs[*viIt]) > 0) {
                vi = *viIt;
                break;
            }
            ++viIt;
        }

        //cout << "vi: " << vi << endl;
        //cout << "p: " << p << endl;

        for (int j = 0; j < p; ++j) {
            std::memset(&cover_degrees[0], 0xFF, graph[vi].size()*sizeof(uint32_t));
            {
                int i = 0;
                for (int u : graph[vi]) {
                    if (!s[u]) {
                        //cout << "u cover_degrees: " << u << endl;
                        cover_degrees[i] = cover_degree(u, graph, s, hs);
                    }
                    ++i;
                }
            };
            /*cout << "cover_degrees: ";
            for (auto cd : cover_degrees) {
                cout << cd << ",";
            }
            //cout << endl;*/

            uint32_t cdmax = [&]{
                uint32_t cdmax = 0;
                for (int i = 0; i < graph[vi].size(); ++i) {
                    if (cover_degrees[i] != -1)
                        cdmax = std::max(cdmax, cover_degrees[i]);
                }
                return cdmax;
            }();

            //cout << "cdmax: " << cdmax << endl;
            uint32_t u = [&] {
                uint32_t umax = 0;
                uint32_t needumax = 0;

                auto uit = graph[vi].begin();
                int i = 0;
                while (uit != graph[vi].end()) {
                    if (cover_degrees[i] != -1 && cover_degrees[i] == cdmax) {
                        auto needu = need_degree(*uit, graph, s, hs);
                        //cout << "uit: " << *uit << " needu " << needu << endl;
                        if (needumax < needu || (needumax == needu && graph[*uit].size() > graph[umax].size())){
                            needumax = needu;
                            umax = *uit;
                        }
                    }

                    ++i;
                    ++uit;
                }
                return umax;
            }();
            //cout << "u: " << u << endl;

            if (!s[u]) {
                s[u] = true;
                for (auto i : graph[u]) {
                    --hs[i];
                }
            }

        }
        ++viIt;
    }
    reduce(s, graph, hs);
}

// https://mdpi-res.com/d_attachment/algorithms/algorithms-14-00079/article_deploy/algorithms-14-00079-v2.pdf
inline std::pair<Solution, Hs> greedy_genmpid_impl(Graph const& graph) {
    Solution s = graphPrunning(graph); // O(|V|)
    //cout << "BEGIN: " << printSol(s) << endl;

    Hs hs(s.size());
    for (int i = 0; i < hs.size(); ++i) {
        hs[i] = hsCompute(i, graph, s);
    }

    std::vector<uint32_t> cover_degrees(graph.size()+1);
    std::vector<uint32_t> c;
    c.reserve(graph.size());  // O(|V|)
    for (int i = 0; i < s.size(); ++i) {  // O(|V|)
        if (hs[i] > 0) {
            c.push_back(i);
        }
    }

    Solution currBestS(s.size());
    Hs currentBestHs(hs.size());
    int count = 0;
    auto sc = s;
    auto hsc = hs;
    #define CPY(v1, v2) std::memcpy(v1.data(), v2.data(), v2.size()*sizeof(decltype(v2)::value_type))
    implOfficial(sc, hsc, c, cover_degrees, graph);
    CPY(currBestS, sc);
    CPY(currentBestHs, hsc);
    count = std::accumulate(currBestS.begin(), currBestS.end(), 0);
    {
    CPY(sc, s);
    CPY(hsc, hs);
    implVariation1(sc, hsc, c, cover_degrees, graph);
    if (int countCand = std::accumulate(sc.begin(), sc.end(), 0);
        countCand < count) {
            count = countCand;
            CPY(currBestS, sc);
            CPY(currentBestHs, hsc);
        }
    }
    {
    CPY(sc, s);
    CPY(hsc, hs);
    implVariation2(sc, hsc, c, cover_degrees, graph);
    if (int countCand = std::accumulate(sc.begin(), sc.end(), 0);
        countCand < count) {
            count  = countCand;
            CPY(currBestS, sc);
            CPY(currentBestHs, hsc);
        }
    }
    {
    CPY(sc, s);
    CPY(hsc, hs);
    implVariation3(sc, hsc, c, cover_degrees, graph);
    if (int countCand = std::accumulate(sc.begin(), sc.end(), 0);
        countCand < count) {
            count  = countCand;
            CPY(currBestS, sc);
            CPY(currentBestHs, hsc);
        }
    }
    #undef CPY
    return {std::move(currBestS), std::move(currentBestHs)};
}

inline Solution greedy_genmpid(Graph const& g) {return greedy_genmpid_impl(g).first;}
}