#ifndef _H_PPLL_LABELING_PARALLEL_BIT_PARALLEL_LABELING_LIGHT_H_
#define _H_PPLL_LABELING_PARALLEL_BIT_PARALLEL_LABELING_LIGHT_H_

#include <mutex>
#include <cilk/cilk.h>

class ParallelBitParallelLabeling_light: public Labeling{
public:	

	explicit ParallelBitParallelLabeling_light(vector<int> &node_offset,
										 vector<int> adjacencyList,
										 bool ownLabel = 1)
		: Labeling(node_offset, adjacencyList, 1) {}


	void BuildLabels() {
		AllocateMemory();
		
		double start_time = -utils::GetCurrentTimeSec();
		index_ = (index_t*)memalign(64, nodes * sizeof(index_t));
		cout << "Time taken to allocate the bp-index: "
			 << start_time + utils::GetCurrentTimeSec() << endl;
	   
	    cilk_for (int j = 0; j < nodes; ++j)
		    for (int i = 0; i < kNumBitParallelRoots; ++i)
	    		index_[j].bpspt_d[i] = MAX_DIST;
	    cout << "Time taken to set-up the bp-index: "
	    	 << start_time + utils::GetCurrentTimeSec() << endl;

	    int root_giver = 0;
	    int rootsLeft = kNumBitParallelRoots;
	    std::mutex lck,lck2;

	    bool minit[kNumBitParallelRoots];
	    vector<uint32_t> mtmp_d[kNumBitParallelRoots];
	    vector<std::pair<uint64_t, uint64_t> > mtmp_s[kNumBitParallelRoots];
	    vector<int> mque[kNumBitParallelRoots];

	    stack<int> mem_id;
	    for (int i = kNumBitParallelRoots-1; i >= 0; --i) {
	    	mem_id.push(i);
	    	minit[i] = false;
	    }


	    int g = min(nodes, kNumBitParallelRoots * 64 + 5);
	    cilk_for(int run = 0; run < g; ++run) {
		    lck.lock();
		    int i_bpspt = kNumBitParallelRoots - rootsLeft;
		    rootsLeft--;
		    while (root_giver < nodes && used[root_giver]) root_giver++;
		    int r = root_giver;
		    root_giver++;
		    if (rootsLeft<0)rootsLeft=0;
		    printf("\r[%2.1lf%%]",(kNumBitParallelRoots-rootsLeft)*100.0/kNumBitParallelRoots);
		    
		    lck.unlock();

		    if (r == nodes) continue;
		    if (rootsLeft <= 0) continue;

		    lck2.lock();
			int thread_mem_id = mem_id.top(); mem_id.pop();
		    lck2.unlock();

		    if (minit[thread_mem_id] == false) {
		    	minit[thread_mem_id] = 1;
		    	mtmp_d[thread_mem_id] = vector<uint32_t>(nodes);
		    	mtmp_s[thread_mem_id] = vector<std::pair<uint64_t, uint64_t> >(nodes);
		    	mque[thread_mem_id] = vector<int>(nodes);
		    }
		   
	      	used[r] = true;

	      	std::vector<uint32_t> &tmp_d = mtmp_d[thread_mem_id];
	      	std::vector<std::pair<uint64_t, uint64_t> > &tmp_s = mtmp_s[thread_mem_id];
	      	std::vector<int> &que = mque[thread_mem_id];
	      
	      	cilk_for(int i = 0; i < nodes; ++i) tmp_d[i] = MAX_DIST, tmp_s[i] = std::make_pair(0,0);

	      	int que_t0 = 0, que_t1 = 0, que_h = 0;
	      	que[que_h++] = r;
	     	tmp_d[r] = 0;
	      	que_t1 = que_h;

	      	int ns = 0;
	     	int cnte = offset[r+1] - offset[r];
	     	sort(adjacency.begin() + offset[r], adjacency.begin() + offset[r] + cnte);

	     	lck.lock();
	      	for (size_t i = 0; i < cnte; ++i) {
	        	int v = adjacency[offset[r] + i];
	        	if (!used[v]) {
	          		used[v] = true;
	          		//que[que_h++] = v;
	          		//tmp_d[v] = 1;
	          		tmp_s[v].first = 1ULL << ns;
	          		if (++ns == 64) break;
	        	}
	      	}
	      	lck.unlock();

	      	for (int d = 0; que_t0 < que_h; ++d) {
	      		assert(d<MAX_DIST);
	        	int numsibling_es = 0, numchild_es = 0;

	        	for (int que_i = que_t0; que_i < que_t1; ++que_i) {
	          		int v = que[que_i];

	          		int cnte = offset[v+1] - offset[v];
	          		for (size_t i = 0; i < cnte; ++i) {
	            		int tv = adjacency[offset[v] + i];//
	            		int td = d + 1;

	            		if (d > tmp_d[tv]);
	            		else if (d == tmp_d[tv]) {
	              			if (v < tv) {
	              				int w = tv;
	              				tmp_s[v].second |= tmp_s[w].first;
	          					tmp_s[w].second |= tmp_s[v].first;
	              			}
	            		} else {
	              			if (tmp_d[tv] == MAX_DIST) {
	                			que[que_h++] = tv;
	                			tmp_d[tv] = td;
	              			}
	              		}
	            	}
	        	}


	        	for (int que_i = que_t0; que_i < que_t1; ++que_i) {
	          		int v = que[que_i];

	          		int cnte = offset[v+1] - offset[v];
	          		for (size_t i = 0; i < cnte; ++i) {
	            		int tv = adjacency[offset[v] + i];//
	            		int td = d + 1;

	            		if (tmp_d[tv] == td) {
	            			int c = tv;
	            			tmp_s[c].first  |= tmp_s[v].first;
	          				tmp_s[c].second |= tmp_s[v].second;
	            		}
	          		}
	        	}

	        	que_t0 = que_t1;
	        	que_t1 = que_h;
	      	}
	      

	      	cilk_for (int v = 0; v < nodes; ++v) {
	        	index_[v].bpspt_d[i_bpspt] = tmp_d[v];
	        	index_[v].bpspt_s[i_bpspt][0] = tmp_s[v].first;
	        	index_[v].bpspt_s[i_bpspt][1] = tmp_s[v].second & ~tmp_s[v].first;
	      	}

	      	lck2.lock();
	      	mem_id.push(thread_mem_id);
	      	lck2.unlock();
	    }printf("\n");

	}
};

#endif  // _H_PPLL_LABELING_PARALLEL_BIT_PARALLEL_LABELING_LIGHT_H_
