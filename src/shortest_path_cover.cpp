#include "cover_problem.hpp"
#include<unordered_map>

int main(int argc, char* argv[])
{

        Env = TEnv(argc, argv, TNotify::StdNotify);
        Env.PrepArgs(TStr::Fmt("cpm. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
        TExeTm ExeTm;
        const int n = Env.GetIfArgPrefixInt("-n:", 100, "No of nodes");
        const double p =  Env.GetIfArgPrefixFlt("-p:", 0.5,"probability");
        const int no_runs = Env.GetIfArgPrefixInt("-noruns:", 100, "No. of runs");

        PUNGraph G;
        int no_nodes = n, no_edges;
        double avg_cover_prop = 0, avg_cover_prop_temp = 0, avg_cover = 0, avg_cover_temp = 0;
        std::unordered_map <std::pair<int,int>,float> Attr;
        for(int i = 0; i < no_runs; ++i)
        {
                printf("run: %d\n",i);
                G = GenErdosRenyiAttr(n,p, Attr);
                printf("IsConnected(G): %d\n",TSnap::IsConnected(G));
                // write_dot(G,"G.dot","G");
                // for (TUNGraph::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++) {
                //         printf("edge (%d, %d)\n", EI.GetSrcNId(), EI.GetDstNId());
                // }
                // G = TSnap::GenRndPowerLaw(n, 3, true);
                no_edges = G->GetEdges();
                // std::cout<<no_edges<<std::endl;
                std::tie(avg_cover_prop_temp, avg_cover_temp)
                        = GetWeightedShortestPathCover(G, no_nodes, no_edges, Attr);
                avg_cover_prop += avg_cover_prop_temp;
                avg_cover += avg_cover_temp;
        }
        avg_cover /= no_runs;
        avg_cover_prop /= no_runs;
        printf("===============================================================\n");
        printf("n: %d, p: %f, no_runs: %d\n",n,p,no_runs);
        printf("Cover_prop: %f, Cover: %f, Expected_no_edges: %f\n",avg_cover_prop,
               avg_cover,((double)(n*(n-1))*p/2));
        printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TExeTm::GetCurTm());
        return 0;
}
