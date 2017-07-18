//
//  Created by Jithin K Sreedharan on 2/20/17.
//  Copyright © 2017 Jithin K Sreedharan. All rights reserved.
//

#ifndef cover_problem_hpp
#define cover_problem_hpp

#include <stdio.h>
#include "Snap.h"
#include <iostream>
#include <vector>
#include <random>
#include <iterator>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <cmath>
#include <queue>
#include <set>


#undef max
#undef min

typedef THash<TUInt64, TUInt64>  TUInt64UInt64H;
typedef THash<TUInt64, TUInt64V> TUInt64UInt64VH;

void print_TIntV(TIntV v, char* S);
void print_TUInt64V(TUInt64V v, char* S);
template<typename Iter, typename RandomGenerator> Iter select_randomly(Iter start, Iter end, RandomGenerator& g);
template<typename Iter> Iter select_randomly(Iter start, Iter end);
int rangeRandom(int min, int max);
PNEANet GenPrefAttachGeneral(const int& time_n, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma);
PNEANet GenPrefAttachGeneral_undirected(const int& time_n, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma, const bool& self_loops_allowed);
void write_dot(PUNGraph G,const char *FName_t,const char *Desc);
PUNGraph GenErdosRenyi(const int& n, const double& p);
double GetWeightedShortestPathCover(const PUNGraph G, const int& no_nodes, const int& no_edges, const int& Attr);

// Template definitions
template <class PGraph>
inline void PrintGStats(const char s[], PGraph Graph) {
        printf("graph %s, nodes %d, edges %d, empty %s\n",
               s, Graph->GetNodes(), Graph->GetEdges(),
               Graph->Empty() ? "yes" : "no");
}

template<typename Iter, typename RandomGenerator>
inline Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
        std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
        std::advance(start, dis(g));
        return start;
}

template<typename Iter>
inline Iter select_randomly(Iter start, Iter end) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        return select_randomly(start, end, gen);
}

#endif /* cover_problem_hpp */
