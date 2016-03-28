#include <atomic>
#include <bits/stdc++.h>
#include <unordered_map>

#include "labeling.h"
#include "graph_parser.h"
#include "labeling/parallel-pruned-labeling.h"
#include "labeling/parallel-bit-parallel-labeling_memory-light.h"
#include "labeling/parallel-bit-parallel-labeling.h"

#ifdef BENCH_STANDARD_PLL
	#include "pruned-landmark-labeling/src/pruned_landmark_labeling.h"
#endif

using namespace std;

typedef long long ll;

namespace {
const bool kBenchMemoryLight = 1;

int number_of_nodes;
vector<pair<int,int> > edges;
vector<int> original_order_index;

// Used for making the compiler not optimize out benchmarking.
// We use this variable to sum all the distances. This way,
// the compiler won't remove the computation, which is not
// creating any side-effects.
int sum_negate_c_opt;

#ifndef BENCH_STANDARD_PLL
	#define kNumBitParallelRoots 100
#endif 

#ifdef BENCH_STANDARD_PLL
	PrunedLandmarkLabeling<kNumBitParallelRoots> pll;

	// Compares the output of the PLL and your parallel PLL.
	// The test is done over randomly generated pairs of nodes.
	bool TestCorrectness(Labeling *labeling) {
		bool test_pass = 1;
		for (int t = 0; t < 100000; ++t) {
			int a = rand()%number_of_nodes;
			int b = rand()%number_of_nodes;
			int d1 = pll.QueryDistance(original_order_index[a],
							           original_order_index[b]);
			int d2 = labeling->QueryDistanceCacheEfficient(a,b);
			if (d1 != d2) {test_pass = false; break;}
		}
		return test_pass;
	}

	// Returns the average time a query is taking within PLL.
	double GetAvgQueryTimePLL() {
		int cnt = 100000;
		double ts = -utils::GetCurrentTimeSec();
		for (int t = 0; t < cnt; ++t) {
			int a = rand()%number_of_nodes;
			int b = rand()%number_of_nodes;
			int d = pll.QueryDistance(a,b);
			sum_negate_c_opt += d;
		}
		ts += utils::GetCurrentTimeSec();
		return ts / cnt;
	}
#endif /** BENCH_STANDARD_PLL **/

// Returns the average time a query is taking within our parallel PLL.
double GetAvgQueryTime(Labeling *labeling) {
	int cnt = 100000;
	double ts = -utils::GetCurrentTimeSec();
	for (int t = 0; t < cnt; ++t) {
		int a = rand()%number_of_nodes;
		int b = rand()%number_of_nodes;
		int d = labeling->QueryDistanceCacheEfficient(a,b);
		sum_negate_c_opt += d;
	}
	ts += utils::GetCurrentTimeSec();
	return ts / cnt;
}
}  // namespace


int parallel_main(int argc, char *argv[]) {
	double exec_start_time = utils::GetCurrentTimeSec();
	int inputBPRootCount;
	ReadGraph(number_of_nodes, edges, inputBPRootCount);

	// No specific number of bit-parallel roots provided? Take the default one.
	//if (kNumBitParallelRoots == -1) kNumBitParallelRoots = kNumBitParallelRoots_ARRSZ;

	utils::RemoveDoubleEdges(edges);
	vector<pair<int,int> > bi_edges = edges;
	utils::AddBackEdges(bi_edges);
	cout << "GRAPH: nodes(" << number_of_nodes << "), edges("<<bi_edges.size()<<")"<<endl;
	cout << "CHECKPOINT(READING DONE): "
		 << (utils::GetCurrentTimeSec() - exec_start_time) << " seconds" << endl;
	cout.flush();

	
	// CREATING ADJACENCY GRAPH //
	vector<int> offset;
	vector<int> adjacency;
	utils::EdgeArrayToAdjacencyGraph(bi_edges, offset, adjacency);


	original_order_index = vector<int>(number_of_nodes+5,0);
	utils::RearrangeGraph(offset, adjacency, original_order_index);
	
#ifdef BENCH_STANDARD_PLL
		cout << "----------------------- STANDARD PLL SCORE -------------------- " << endl;
		pll.ConstructIndex(edges);
		pll.PrintStatistics();
		cout << "AVG QUERY: "<< GetAvgQueryTimePLL()*1000000 << " micro sec"<< endl;
		cout << "----------------------------------------------------------- " << endl;
#endif

	// timeer reset for exec analysis
	exec_start_time = utils::GetCurrentTimeSec();

	cout << "Number of bit-parallel roots: "<< kNumBitParallelRoots << endl << endl;
	double startTime = utils::GetCurrentTimeSec();
	Labeling *first_stage;

	if (kBenchMemoryLight)
	  first_stage = new ParallelBitParallelLabeling_light(offset, adjacency);
	else
		first_stage = new ParallelBitParallelLabeling(offset, adjacency);


	cout << "First stage initialization: " << (utils::GetCurrentTimeSec()-exec_start_time) << endl;
	first_stage->BuildLabels();
	first_stage->clearTempMem();
	cout << "First stage done: "<< (utils::GetCurrentTimeSec()-exec_start_time) << endl;


	Labeling *second_stage = new ParallelPrunedLabeling(offset, adjacency, false);
	cout << endl<< "Second stage initialization: "
		 << (utils::GetCurrentTimeSec() - exec_start_time) << endl;
	second_stage->takeLabelData(first_stage);
	second_stage->BuildLabels();
	cout << "Second stage done: "<< (utils::GetCurrentTimeSec()-startTime)
		 << " (avg. label size: " << second_stage->GetAvgLabelSize() << ")"<<endl;
	double build_time = utils::GetCurrentTimeSec() - exec_start_time;
	double avg_label_size = second_stage->GetAvgLabelSize();
	second_stage->fixVertexOrder();


	cout << endl << "Final-time: "<< build_time << endl;
	cout << "Avg-label-size: " << avg_label_size <<endl;
	cout << "Query-speed: " << GetAvgQueryTime(second_stage)*1000000 <<endl; // micro-sec
	
#ifndef BENCH_STANDARD_PLL
		cout << "Testing skipped... (make kBenchStandardPLL=true for testing)!" << endl;
#else
		cout << "Correctness: "
		     << TestCorrectness(second_stage) << endl;	
#endif
	cout << "(Ignore this value: " << sum_negate_c_opt << ")" << endl;
	return 0;
}