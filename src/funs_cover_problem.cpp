#include "cover_problem.hpp"
#include<unordered_map>
#include<hash>

struct SimpleHash
{
    size_t operator()(const pair<int,int>&x)const
    {
        return hash<long long>()(((long long)x.first)^(((long long)x.second)<<32));
    }
};

void write_dot(PNEANet G,const char *FName_t,const char *Desc)
{
        // Output node IDs as numbers
        TIntStrH NIdLabelH;
        // Generate labels for random graph
        for (TNEANet::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
                NIdLabelH.AddDat(NI.GetId(), TStr::Fmt("Node%d", NI.GetId()));
        }
        TSnap::SaveGViz(G, FName_t, Desc,NIdLabelH);
}

// template <class PGraph>
// void write_dot(PGraph G,const char *FName_t,const char *Desc)
// {
//     // Output node IDs as numbers
//     TIntStrH NIdLabelH;
//     // Generate labels for random graph
//     for (typename PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
//       NIdLabelH.AddDat(NI.GetId(), TStr::Fmt("Node%d", NI.GetId()));
//     }
//     TSnap::SaveGViz(G, FName_t, Desc,NIdLabelH);
// }

void write_dot(PUNGraph G,const char *FName_t,const char *Desc)
{
        // Output node IDs as numbers
        TIntStrH NIdLabelH;
        // Generate labels for random graph
        for (TUNGraph::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
                NIdLabelH.AddDat(NI.GetId(), TStr::Fmt("Node%d", NI.GetId()));
        }
        TSnap::SaveGViz(G, FName_t, Desc,NIdLabelH);
}

uint64_t longest_distance(const PNEANet & G_DAG, const uint64_t &G_no_nodes, uint64_t node_u)
{
        uint64_t node_v = 0;
        uint64_t len_long_path = 0, len_long_path_t = 0;
        uint64_t max;
        TNEANet::TNodeI NI_temp;
        if (node_u == 0)
        {
                len_long_path = 1;
                return len_long_path;
        }
        else
        {
                NI_temp = G_DAG->GetNI(node_u);
                max = 0;
                for(int i = 0; i < NI_temp.GetOutDeg(); i++)
                {
                        node_v = (uint64_t)NI_temp.GetOutNId(i);
                        len_long_path_t = longest_distance(G_DAG, G_no_nodes, node_v) + 1;
                        if (max < len_long_path_t)
                                max = len_long_path_t;
                }
                len_long_path += max;
                return len_long_path;
        }
}

void print_TIntV(TIntV v, char* S)
{
        printf("%s: ", S);
        for (int i =0; i < v.Len(); i++)
        {
                printf("%d ", v[i].Val);
        }
        printf("\n");
}

void print_TUInt64V(TUInt64V v, char* S)
{
        printf("%s: ", S);
        for (int i =0; i < v.Len(); i++)
        {
                printf("%d ", v[i].Val);
        }
        printf("\n");
}


int rangeRandom(int min, int max)
{
        int n = max - min + 1;
        int remainder = RAND_MAX % n;
        int x;
        do {
                x = rand();
        } while (x >= RAND_MAX - remainder);
        return min + x % n;
}


PNEANet GenPrefAttachGeneral(const int& time_n, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2,
                             const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma)
{
        PNEANet G = PNEANet::New();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> rv_m_new(vec_p_1,vec_p_2); //unif distbn is a class
        std::uniform_int_distribution<> rv_m_old(vec_q_1,vec_q_2);
        std::bernoulli_distribution rv_alpha(pr_alpha);
        std::bernoulli_distribution rv_beta(pr_beta);
        std::bernoulli_distribution rv_delta(pr_delta);
        std::bernoulli_distribution rv_gamma(pr_gamma);


        std::vector<int> target_list;
        std::vector<int> target_list_temp;
        int target;
        uint16_t m;

        G->AddNode(0);
        m = rv_m_new(gen);
        for(int i = 0; i < m; i++)
        {
                G->AddEdge(0,0);
                target_list.push_back(0);
                target_list.push_back(0);
        }
        G->AddNode(1);
        m = rv_m_new(gen);
        for(int i = 0; i < m; i++)
        {
                G->AddEdge(1,0);
                target_list.push_back(0);
                target_list.push_back(1);
        }
        int source = 2;
        int source_exist;

        TNEANet::TNodeI NI_temp;


        for (int time_i = 2; time_i < time_n; time_i++)
        {
                if(rv_alpha(gen))
                {
                        target_list_temp.clear();
                        G->AddNode(source);
                        m = rv_m_new(gen);
                        for(int i =0; i<m; i++)
                        {
                                if(rv_beta(gen))
                                {
                                        target = *select_randomly(target_list.begin(), target_list.end());
                                        G->AddEdge(source,target);
                                        target_list_temp.push_back(source);
                                        target_list_temp.push_back(target);
                                }
                                else{
                                        //                    target = rangeRandom(0,(source-1));
                                        std::uniform_int_distribution<int> random_node(0,(source-1)); // guaranteed unbiased
                                        target = random_node(gen);
                                        G->AddEdge(source,target);
                                        target_list_temp.push_back(source);
                                        target_list_temp.push_back(target);
                                }
                        }
                        source++;
                        target_list.insert(target_list.end(), target_list_temp.begin(), target_list_temp.end());
                }
                else{
                        target_list_temp.clear();
                        m = rv_m_old(gen);
                        for(int i =0; i<m; i++)
                        {
                                if (rv_delta(gen))
                                {
                                        source_exist = *select_randomly(target_list.begin(), target_list.end());
                                        //* is added because it returns a pointer Iter
                                }
                                else{
                                        std::uniform_int_distribution<int> random_node(0,(source-1)); // guaranteed unbiased
                                        source_exist = random_node(gen);
                                }
                                if (rv_gamma(gen))
                                {
                                        target = *select_randomly(target_list.begin(), target_list.end());
                                }
                                else{
                                        std::uniform_int_distribution<int> random_node(0,(source-1)); // guaranteed unbiased
                                        target = random_node(gen);
                                }
                                G->AddEdge(source_exist,target);
                                target_list_temp.push_back(source_exist);
                                target_list_temp.push_back(target);
                        }
                        target_list.insert(target_list.end(), target_list_temp.begin(), target_list_temp.end());
                }
        }
        return G;
}

PNEANet GenPrefAttachGeneral_undirected(const int& time_n, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma, const bool& self_loops_allowed)
{
        PNEANet G = PNEANet::New();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> rv_m_new(vec_p_1,vec_p_2); //unif distbn is a class
        std::uniform_int_distribution<> rv_m_old(vec_q_1,vec_q_2);
        std::bernoulli_distribution rv_alpha(pr_alpha);
        std::bernoulli_distribution rv_beta(pr_beta);
        std::bernoulli_distribution rv_delta(pr_delta);
        std::bernoulli_distribution rv_gamma(pr_gamma);


        std::vector<int> target_list;
        std::vector<int> target_list_temp;
        int target;
        uint16_t m;

        G->AddNode(0);
        m = rv_m_new(gen);
        if (self_loops_allowed)
        {
                for(int i = 0; i < m; i++)
                {
                        G->AddEdge(0,0);
                        target_list.push_back(0);
                        target_list.push_back(0);
                }
        }
        G->AddNode(1);
        m = rv_m_new(gen);
        for(int i = 0; i < m; i++)
        {
                G->AddEdge(1,0);
                G->AddEdge(0,1);
                target_list.push_back(0);
                target_list.push_back(1);
        }
        int source = 2;
        int source_exist;

        TNEANet::TNodeI NI_temp;


        for (int time_i = 2; time_i < time_n; time_i++)
        {
                if(rv_alpha(gen))
                {
                        target_list_temp.clear();
                        G->AddNode(source);
                        m = rv_m_new(gen);
                        for(int i =0; i<m; i++)
                        {
                                if(rv_beta(gen))
                                {
                                        target = *select_randomly(target_list.begin(), target_list.end());
                                        G->AddEdge(source,target);
                                        G->AddEdge(target,source);
                                        target_list_temp.push_back(source);
                                        target_list_temp.push_back(target);
                                }
                                else
                                {
                                        //                    target = rangeRandom(0,(source-1));
                                        std::uniform_int_distribution<int> random_node(0,(source-1)); // guaranteed unbiased
                                        target = random_node(gen);
                                        G->AddEdge(source,target);
                                        G->AddEdge(target,source);
                                        target_list_temp.push_back(source);
                                        target_list_temp.push_back(target);
                                }
                        }
                        source++;
                        target_list.insert(target_list.end(), target_list_temp.begin(), target_list_temp.end());
                }
                else
                {
                        target_list_temp.clear();
                        m = rv_m_old(gen);
                        for(int i =0; i<m; i++)
                        {
                                if (rv_delta(gen))
                                {
                                        source_exist = *select_randomly(target_list.begin(), target_list.end());
                                        //* is added because it returns a pointer Iter
                                }
                                else
                                {
                                        std::uniform_int_distribution<int> random_node(0,(source-1)); // guaranteed unbiased
                                        source_exist = random_node(gen);
                                }
                                if (rv_gamma(gen))
                                {
                                        target = *select_randomly(target_list.begin(), target_list.end());
                                }
                                else
                                {
                                        std::uniform_int_distribution<int> random_node(0,(source-1)); // guaranteed unbiased
                                        target = random_node(gen);
                                }
                                G->AddEdge(source_exist,target);
                                G->AddEdge(target,source_exist);
                                target_list_temp.push_back(source_exist);
                                target_list_temp.push_back(target);
                        }
                        target_list.insert(target_list.end(), target_list_temp.begin(), target_list_temp.end());
                }
        }
        return G;
}

std::pair<int,int> my_make_pair(int a, int b)
{
        if ( a < b ) return std::pair<int,int>(a,b);
        else return std::pair<int,int>(b,a);
}

PUNGraph GenErdosRenyiAttr(const int& n, const double& p,
                           std::map<std::pair<int,int>,float>& Attr)
{
        std::random_device rd;
        std::mt19937 gen(rd());
        std::bernoulli_distribution rv_b(p);
        std::uniform_real_distribution<> rv_u(0,1);
        PUNGraph G = PUNGraph::New();
        for(int i = 0; i < n; i++)
                G->AddNode(i);

        for(int i = 0; i < n; ++i)
        {
                for(int j = i+1; j < n; ++j)
                {
                        if(rv_b(gen)) {
                                G->AddEdge(i, j);
                                Attr[my_make_pair(i, j)] = rv_u(gen);
                        }
                }
        }
        return G;
}

int findMinimum_jks(TIntV& Frontier, TIntFltH& NIdDistH) {
        TFlt minimum = TInt::Mx;
        int min_index = 0;
        for (int i = 0; i < Frontier.Len(); i++) {
                int NId = Frontier.GetVal(i);
                if (NIdDistH[NId] < minimum) {
                        minimum = NIdDistH[NId];
                        min_index = i;
                }
        }
        const int NId = Frontier.GetVal(min_index);
        Frontier.Del(min_index);
        return NId;
}

//My implementation of SNAP function
int GetWeightedShortestPath_jks(
        const PUNGraph Graph, const int& no_nodes, const int& SrcNId,
        std::set<std::pair<int,int> >& covered_edges,
        std::map<std::pair<int,int>,float>& Attr) {
        TIntV frontier;
        TIntH prev;
        TIntFltH NIdDistH;
        NIdDistH.Clr(false); NIdDistH.AddDat(SrcNId, 0);
        frontier.Add(SrcNId); //NOTE: Is it needed? Commented it
        // int Attr = 1; //NOTE: Changed Attr vector to a scalar in jks implementation
        TFlt attr_temp;
        while (!frontier.Empty()) {
                const int NId = findMinimum_jks(frontier, NIdDistH);
                const PUNGraph::TObj::TNodeI NodeI = Graph->GetNI(NId);
                for (int v = 0; v < NodeI.GetOutDeg(); v++) {
                        int DstNId = NodeI.GetOutNId(v);
                        // int EId = NodeI.GetOutEId(v);
                        if (DstNId > SrcNId) {
                                // if (1) {
                                attr_temp = Attr[my_make_pair(NId,DstNId)];
                                if (!NIdDistH.IsKey(DstNId)) {
                                        NIdDistH.AddDat(DstNId, NIdDistH.GetDat(NId)+ attr_temp);
                                        frontier.Add(DstNId);
                                        prev.AddDat(DstNId,NId);
                                }
                                else {
                                        if (NIdDistH.GetDat(DstNId) > NIdDistH.GetDat(NId) + attr_temp) {
                                                NIdDistH.AddDat(DstNId,NIdDistH.GetDat(NId) + attr_temp);
                                                prev.AddDat(DstNId,NId);
                                        }
                                }
                        }
                }
        }
        //**TEST**
        // int Key,Value;
        // for (TIntIntH::TIter It = prev.BegI(); It < prev.EndI(); It++) {
        //         // get the key
        //         Key = It.GetKey();
        //         // get the value
        //         Value = It.GetDat();
        //         printf("(%d,%d) ", Key,Value);
        // }
        // printf("\n" );
        //********
        int u,v;
        for (int uu = SrcNId+1; uu < no_nodes; ++uu) {
                u = uu;
                while(prev.IsKey(u)) {
                        v = prev.GetDat(u);
                        covered_edges.insert(my_make_pair(u,v));
                        u = v;
                }
        }
        return 0;
}

/*int getWeightedShortestPath_jks(const PUNGraph Graph, const int& no_nodes, const int& SrcNId,
                                std::unordered_set<std::pair<int,int>, SimpleHash >& covered_edges,
                                std::unordered_map<std::pair<int,int>,float, SimpleHash>& Attr)
{
    auto set<std::pair<float, int>> Q; //dist, nodeID
    auto unordered_map<std::pair<int,int>,SimpleHash> prev;
    auto unordered_map<std::pair<int,int>,SimpleHash> dist;
    dist.insert(SrcNId, 0);
    Q.insert(SrcNId, 0);
    double attr_temp;

    while(!Q.Empty())
    {
        std::set<std::pair<pair, int>::iterator it = Q.begin();
        const int NId = *it->second();
        const PUNGraph::TObj::TNodeI NodeI = Graph->GetNI(NId);
        for (int v = 0; v < NodeI.GetOutDeg(); v++)
        {
            int DstNId = NodeI.GetOutNId(v);
            if(DstNId > SrcNId)
            {
                attr_temp = Attr[(NId ^ DstNId)];
                std::unordered_map<int,int>::const_iterator got = dist.find(DstNId);
                if (!(got == dist.end()))
                {
                    dist.insert(my_make_pair(DstNId, NId + attr_temp));
                    Q.insert(DstNId);
                    prev.insert(DstNId,NId);
                }
                else if(got->first  > got->second + attr_temp)
                {
                    dist.insert(my_make_pair())
                    prev.insert(my_make_pair(DstNId, NId));
                }

            }
        }
    }

    int u ,v;
    for (int uu = SrcNId+1; uu < no_nodes; ++uu)
    {
        u = uu;
        std::unordered_map<int,int>::const_iterator got = dist.find(u)
        std::unordered_map<int,int>::const_iterator previ = prev.find(u)
        while(!(got == dist.end()))
        {
            v = previ->second;
            covered_edges.insert(my_make_pair(u,v));
            u = v;
        }
    }

} */

std::pair<double,double> GetWeightedShortestPathCover(const PUNGraph G, const int& no_nodes, const int& no_edges,
                                                      std::map<std::pair<int,int>,float>& Attr ) {
        std::set<std::pair<int,int> > covered_edges;
        for(int SrcNId = 0; SrcNId < no_nodes; ++SrcNId) {
                GetWeightedShortestPath_jks(G, no_nodes, SrcNId, covered_edges, Attr);
        }
        //**TEST**
        // printf("covered_edges.size():%d, no_edges:%d\n",covered_edges.size(),no_edges);
        double cover_size = covered_edges.size();
        double cover_size_prop = (double)covered_edges.size()/no_edges;
        return std::make_pair(cover_size_prop,cover_size);
}
