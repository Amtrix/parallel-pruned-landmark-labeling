#ifndef _H_PPLL_LABELING_PARALLEL_BIT_PARALLEL_LABELING_H_
#define _H_PPLL_LABELING_PARALLEL_BIT_PARALLEL_LABELING_H_

#include <mutex>
#include <cilk/cilk.h>

class ParallelBitParallelLabeling: public Labeling{
public:

	explicit ParallelBitParallelLabeling(vector<int> &node_offset,
										 vector<int> adjacencyList,
										 bool ownLabel = 1)
		: Labeling(node_offset, adjacencyList, 1) {
		}


	void BuildLabels() {
		AllocateMemory();
		double start_time = -utils::GetCurrentTimeSec();
		index_ = (index_t*)memalign(64, nodes * sizeof(index_t));
		cout << "Time taken to allocate the bp-index: "
			 << start_time + utils::GetCurrentTimeSec() << endl;
		
	    cilk_for (int j = 0; j < nodes; ++j) {
			index_[j] = index_t();
		    for (int i = 0; i < kNumBitParallelRoots; ++i)
	    		index_[j].bpspt_d[i] = MAX_DIST;
	    }
	    cout << "Time taken to set-up the bp-index: "
	    	 << start_time + utils::GetCurrentTimeSec() << endl;

	    int root_giver = 0;
	    int rootsLeft = kNumBitParallelRoots + 1;
	    std::mutex lck,lck2;

	    bool m_init[kNumBitParallelRoots];
	    vector<uint32_t> m_tmp_d[kNumBitParallelRoots];
	    vector<std::pair<uint64_t, uint64_t> > m_tmp_s[kNumBitParallelRoots];
	    vector<int> m_que[kNumBitParallelRoots];
	    vector<std::pair<int,int> > m_sibling_es[kNumBitParallelRoots];
	    vector<std::pair<int,int> > m_child_es[kNumBitParallelRoots];

	    stack<int> mem_ids;
	    for (int i = kNumBitParallelRoots-1; i >= 0; --i) {
	    	mem_ids.push(i);
	    	m_init[i] = false;
	    }


	    int g = min(nodes, kNumBitParallelRoots * 64 + 5);
	    cilk_for(int run = 0; run < g; ++run) {
		    lck.lock();
		    int i_bpspt = kNumBitParallelRoots - rootsLeft + 1;
		    rootsLeft--;
		    if (rootsLeft<0)rootsLeft=0;
		    while (root_giver < nodes && used[root_giver])
		    	root_giver++;
		    int r = root_giver;
		    root_giver++;
		    printf("\r[%2.1lf%%]",(
		    	kNumBitParallelRoots-rootsLeft)*100.0/kNumBitParallelRoots);
		    lck.unlock();

		    if (r >= nodes) continue;
		    if (i_bpspt >= kNumBitParallelRoots) continue;

		    lck2.lock();
			int thread_mem_id = mem_ids.top(); mem_ids.pop();
		    lck2.unlock();

		    if (m_init[thread_mem_id] == false) {
		    	m_init[thread_mem_id] = 1;
		    	m_tmp_d[thread_mem_id] = vector<uint32_t>(nodes);
		    	m_tmp_s[thread_mem_id] = vector<std::pair<uint64_t, uint64_t> >(nodes);
		    	m_que[thread_mem_id] = vector<int>(nodes);
		    	m_sibling_es[thread_mem_id] = vector< pair<int,int> >(edges);
		    	m_child_es[thread_mem_id] = vector< pair<int,int> >(edges);
		    }
		   
	      	used[r] = true;

	      	std::vector<uint32_t> &tmp_d = m_tmp_d[thread_mem_id];
	      	std::vector<std::pair<uint64_t, uint64_t> > &tmp_s = m_tmp_s[thread_mem_id];
	      	std::vector<int> &que = m_que[thread_mem_id];
	      	std::vector<std::pair<int, int> > &sibling_es = m_sibling_es[thread_mem_id];
	      	std::vector<std::pair<int, int> > &child_es = m_child_es[thread_mem_id];

	      	cilk_for(int i = 0; i < nodes; ++i)
	      		tmp_d[i] = MAX_DIST, tmp_s[i] = std::make_pair(0,0);

	      	int que_t0 = 0, que_t1 = 0, que_h = 0;
	      	que[que_h++] = r;
	     	tmp_d[r] = 0;
	      	que_t1 = que_h;

	      	int ns = 0;
	     	int cnte = offset[r+1] - offset[r];
	     	sort(adjacency.begin() + offset[r],
	     		 adjacency.begin() + offset[r] + cnte);

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
	        	int num_sibling_es = 0, num_child_es = 0;

	        	for (int que_i = que_t0; que_i < que_t1; ++que_i) {
	          		int v = que[que_i];

	          		int cnte = offset[v+1] - offset[v];
	          		for (size_t i = 0; i < cnte; ++i) {
	            		int tv = adjacency[offset[v] + i];//
	            		int td = d + 1;

	            		if (d > tmp_d[tv]);
	            		else if (d == tmp_d[tv]) {
	              			if (v < tv) {
	                			sibling_es[num_sibling_es].first  = v;
	                			sibling_es[num_sibling_es].second = tv;
	                			++num_sibling_es;
	              			}
	            		} else {
	              			if (tmp_d[tv] == MAX_DIST) {
	                			que[que_h++] = tv;
	                			tmp_d[tv] = td;
	              			}
	              			child_es[num_child_es].first  = v;
	              			child_es[num_child_es].second = tv;
	              			++num_child_es;
	            		}
	          		}
	        	}
	       

	       	 	for (int i = 0; i < num_sibling_es; ++i) {
	          		int v = sibling_es[i].first, w = sibling_es[i].second;
	          		tmp_s[v].second |= tmp_s[w].first;
	          		tmp_s[w].second |= tmp_s[v].first;
	       	 	}
	       		for (int i = 0; i < num_child_es; ++i) {
	          		int v = child_es[i].first, c = child_es[i].second;
	          		tmp_s[c].first  |= tmp_s[v].first;
	          		tmp_s[c].second |= tmp_s[v].second;
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
	      	mem_ids.push(thread_mem_id);
	      	lck2.unlock();
	    }
	    printf("\n");
	}
};

#endif  // _H_PPLL_LABELING_PARALLEL_BIT_PARALLEL_LABELING_H_
