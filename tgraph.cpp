#include "tgraphreader.h"
#include "tgraph.h"


#include "cppUtils.h"
#include "debug.h"

#include "arraysort.h"
#include <cstdlib>
#include <cstring>



using namespace std;

void tovector(btree_set<uint> &s, vector<uint> &v) {
	for(btree_set<uint>::iterator it=s.begin(); it != s.end(); ++it) {
		v.push_back(*it);
	}
}

void encodediff(vector<uint> &t) {
        uint val, old;
        old = t[0];
        for(uint j=1; j < t.size(); j++) {                        
                val = t[j];
                t[j] -= old;
                old = val;
        }
}

void decodediff(uint *k, uint size) {
        if (size == 0) return;
        
        for(uint i=1; i < size; i++) {
                k[i] += k[i-1];
        }
}


void TGraph::save(ofstream &out) {
        cds_utils::saveValue<TGraph>(out, this, 1);
        
        cds_utils::saveValue<TGraphEventList>(out, tgraph, nodes);
        
        for(uint i=0; i < nodes; i++) {
		if (tgraph[i].neighbors > 0) {
			//printf("node: %u\n", i);
			//printf("neighbors: %u\n", tgraph[i].neighbors);
			//printf("csize_edges: %u\n", tgraph[i].csize_cedges);
			//printf("csize_cchanges: %u\n", tgraph[i].csize_cchanges);
			//printf("csize_ctime: %u\n", tgraph[i].csize_ctime);
			//printf("csize_cedgetimesize: %u\n", tgraph[i].csize_cedgetimesize);
                cds_utils::saveValue<uint>(out, tgraph[i].cedges, tgraph[i].csize_cedges);
                cds_utils::saveValue<uint>(out, tgraph[i].cchanges, tgraph[i].csize_cchanges);
                cds_utils::saveValue<uint>(out, tgraph[i].ctime, tgraph[i].csize_ctime);
                cds_utils::saveValue<uint>(out, tgraph[i].cedgetimesize, tgraph[i].csize_cedgetimesize);
		}
        }

        //save reverse graph
        cds_utils::saveValue<TGraphReverse>(out, reverse, nodes);
        for(uint i=0; i < nodes; i++) {
        	if (reverse[i].size > 0) {
        		cds_utils::saveValue<uint>(out, reverse[i].clist, reverse[i].csize);
        	}
        }
}

TGraph* TGraph::load(ifstream &in) {
        TGraph *tg;
        
        tg = cds_utils::loadValue<TGraph>(in, 1);
        
        LOG("nodes: %u", tg->nodes);
        LOG("edges: %u", tg->edges);
        LOG("changes: %u", tg->changes);
        LOG("maxtime: %u", tg->maxtime);
        
        tg->tgraph = cds_utils::loadValue<TGraphEventList>(in, tg->nodes);
        
        for(uint i=0; i < tg->nodes; i++) {
		if (tg->tgraph[i].neighbors > 0) {
			//printf("node: %u\n", i);
			//printf("csize_edges: %u\n", tg->tgraph[i].csize_cedges);
			
                tg->tgraph[i].cedges = cds_utils::loadValue<uint>(in, tg->tgraph[i].csize_cedges);
                tg->tgraph[i].cchanges = cds_utils::loadValue<uint>(in, tg->tgraph[i].csize_cchanges);
                tg->tgraph[i].ctime = cds_utils::loadValue<uint>(in, tg->tgraph[i].csize_ctime);
                tg->tgraph[i].cedgetimesize = cds_utils::loadValue<uint>(in, tg->tgraph[i].csize_cedgetimesize);
		}
        }
	
        //read reverse graph
        tg->reverse = cds_utils::loadValue<TGraphReverse>(in, tg->nodes);
        for(uint i=0; i < tg->nodes; i++) {
                if (tg->reverse[i].size > 0) {
                	tg->reverse[i].clist =  cds_utils::loadValue<uint>(in, tg->reverse[i].csize);
                }
        }


        tg->loadpolicy();
        
        return tg;
}






void TGraph::create(TGraphReader &tgr) {
        //uint *nodesbuffer = new uint[(BUFFER/BLOCKSIZE+1)*BLOCKSIZE];
        //uint *timebuffer = new uint[(BUFFER/BLOCKSIZE+1)*BLOCKSIZE];
        uint *ccedgebuffer = new uint[(BUFFER/BLOCKSIZE+1)*BLOCKSIZE];
        uint *cctimebuffer = new uint[(BUFFER/BLOCKSIZE+1)*BLOCKSIZE];
        //uint *timebuffer = new uint[(BUFFER/BLOCKSIZE+1)*BLOCKSIZE];
        
        nodes = tgr.nodes;
        edges = tgr.edges;
        changes = tgr.changes;
        maxtime = tgr.maxtime;
        
        tgraph = new TGraphEventList[nodes];
        
        LOG("nodes: %u", nodes);
        LOG("edges: %u", edges);
        LOG("changes: %u", changes);
        LOG("maxtime: %u", maxtime);

        this->loadpolicy(); //load compression policy for time
        
		vector<uint> readset;
        
        uint neighbors; //list of edges
	
	vector<uint> curr_edges;
	vector<uint> curr_changes;
	vector<uint> curr_edgetimesize;
	
	uint csize_time;
        
	uint csize=0;
        for(uint i=0; i < nodes; i++) {
                if (i%1000==0) fprintf(stderr, "Compressing: %0.2f%%\r", (float)i*100/nodes);

		tgraph[i].csize_ctime = 0;
		tgraph[i].csize_cedges = 0;
                tgraph[i].csize_cchanges = 0;
		tgraph[i].csize_cedgetimesize = 0;

		tgraph[i].neighbors = 0;

		tgraph[i].ctime = NULL;
		tgraph[i].cedges = NULL;
                tgraph[i].cchanges = NULL;
		tgraph[i].cedgetimesize = NULL;
		
        neighbors = tgr.tgraph[i].size();
		
		btree_map< uint, btree_set<uint> >::iterator it;
		if (neighbors == 0) continue;
		
                curr_edges.clear();
                curr_changes.clear();
                curr_edgetimesize.clear();
		
                uint j=0;
                csize_time=0;
               
                for (it=tgr.tgraph[i].begin(); it!=tgr.tgraph[i].end(); ++it) {
					readset.clear();
					tovector(it->second, readset);
                        encodediff(readset);
                        
                        csize = cc->Compress(readset.data(), ccedgebuffer, readset.size());
                        
                        memcpy(&cctimebuffer[csize_time], ccedgebuffer, csize * sizeof(uint));
                        
                        curr_edges.push_back(it->first);
						curr_changes.push_back(readset.size());
                        curr_edgetimesize.push_back(csize_time);
                        
                        j++;
                        
                        csize_time += csize;
                        
                        (it->second).clear();
                }
		
		tgraph[i].neighbors = neighbors;

		tgraph[i].csize_ctime = csize_time;
        tgraph[i].ctime = new uint [tgraph[i].csize_ctime];
        memcpy(tgraph[i].ctime, cctimebuffer, csize_time * sizeof(uint));
		

		encodediff(curr_edges);
		tgraph[i].csize_cedges = cc->Compress(curr_edges.data(), ccedgebuffer, neighbors);
		tgraph[i].cedges = new uint[tgraph[i].csize_cedges];
		memcpy(tgraph[i].cedges, ccedgebuffer, tgraph[i].csize_cedges * sizeof(uint));
		
		//printf("node: %u\n", i);
		//printf("csize_edges: %u\n", tgraph[i].csize_cedges);

		tgraph[i].csize_cchanges = cc->Compress(curr_changes.data(), ccedgebuffer, neighbors);
		tgraph[i].cchanges = new uint[tgraph[i].csize_cchanges];
		memcpy(tgraph[i].cchanges, ccedgebuffer, tgraph[i].csize_cchanges * sizeof(uint));
		
		tgraph[i].csize_cedgetimesize = cc->Compress(curr_edgetimesize.data(), ccedgebuffer, neighbors);
		tgraph[i].cedgetimesize = new uint[tgraph[i].csize_cedgetimesize];
		memcpy(tgraph[i].cedgetimesize, ccedgebuffer, tgraph[i].csize_cedgetimesize * sizeof(uint));
		
		
		//clean used variable
		tgr.tgraph[i].clear();
        }
        fprintf(stderr, "\n");
        
		tgr.tgraph.clear();

        // Creating reverse structure
        reverse = new TGraphReverse[nodes];
        uint size;
        for(uint i=0; i < nodes; i++) {
                if (i%1000==0) fprintf(stderr, "Compressing reverse graph: %0.2f%%\r", (float)i*100/nodes);

                reverse[i].size = 0;
                reverse[i].csize = 0;
                reverse[i].clist = NULL;

                size = tgr.revgraph[i].size();

                if ( size == 0 ) {
                	continue;
                }

				readset.clear();
				tovector(tgr.revgraph[i], readset);
                encodediff(readset);
                csize = cc->Compress(readset.data(), ccedgebuffer, size);

                reverse[i].size = size;
                reverse[i].csize = csize;
                reverse[i].clist = new uint[csize];
                memcpy(reverse[i].clist, ccedgebuffer, csize * sizeof(uint));

                //vector<uint>().swap(tgr.revgraph[i].neighbors);
                tgr.revgraph[i].clear();

        }

		tgr.revgraph.clear();


        delete [] ccedgebuffer;
	delete [] cctimebuffer;
}




void TGraph::decodetime(uint u, uint v, uint *edgetimesize, uint *changes, uint *res) {
        //if (tgraph[v].changes == 0) return;
         
        cc->Decompress(&tgraph[u].ctime[ edgetimesize[v] ], res, changes[v]);
        decodediff(res, changes[v]);
}

void TGraph::decodereverse(uint v, uint *res) {
        //if (tgraph[v].changes == 0) return;

        cc->Decompress(reverse[v].clist, res, reverse[v].size);
        decodediff(res, reverse[v].size);
}

void TGraph::direct_point(uint v, uint t, uint *res)  {
        if (v>=nodes || tgraph[v].neighbors == 0) return;

        uint *timep = new uint[BLOCKSIZE*tgraph[v].csize_cedgetimesize];
	uint *changesp = new uint[BLOCKSIZE+tgraph[v].neighbors];
	uint *edgesp = new uint[BLOCKSIZE+tgraph[v].neighbors];
	uint *edgetimesizep = new uint[BLOCKSIZE+tgraph[v].neighbors];
	
	cc->Decompress(tgraph[v].cchanges, changesp, tgraph[v].neighbors);
	cc->Decompress(tgraph[v].cedges, edgesp, tgraph[v].neighbors);
	decodediff(edgesp, tgraph[v].neighbors);
	
	cc->Decompress(tgraph[v].cedgetimesize, edgetimesizep, tgraph[v].neighbors);
	
        uint i=0;
        
        for(uint j=0; j < tgraph[v].neighbors; j++) {
                decodetime(v, j, edgetimesizep, changesp, timep);
                /*
                uint c=0;
                for (uint k=0; k < changesp[j]; k++) {
                        if (timep[k] <= t) {
                                c++;
                        }
                }
                
                if (c%2 == 1) res[++i] = edgesp[j];
                */
                uint *low = lower_bound(timep, timep+changesp[j], t+1);
                if ( (low-timep)%2 == 1) res[++i] = edgesp[j];
                  
                
        }
        
        *res = i;

        delete [] timep;
	delete [] changesp;
	delete [] edgesp;
	delete [] edgetimesizep;
}

void TGraph::direct_interval(uint v, uint tstart, uint tend, uint semantic, uint *res)  {
        if (v>=nodes || tgraph[v].neighbors == 0) return;

        uint *timep = new uint[BLOCKSIZE*tgraph[v].csize_cedgetimesize];
	uint *changesp = new uint[BLOCKSIZE+tgraph[v].neighbors];
	uint *edgesp = new uint[BLOCKSIZE+tgraph[v].neighbors];
	uint *edgetimesizep = new uint[BLOCKSIZE+tgraph[v].neighbors];
	
	cc->Decompress(tgraph[v].cchanges, changesp, tgraph[v].neighbors);
	cc->Decompress(tgraph[v].cedges, edgesp, tgraph[v].neighbors);
	decodediff(edgesp, tgraph[v].neighbors);
	
	cc->Decompress(tgraph[v].cedgetimesize, edgetimesizep, tgraph[v].neighbors);
        
	
        uint i=0;
        
        for(uint j=0; j <  tgraph[v].neighbors; j++) {
                decodetime(v, j, edgetimesizep, changesp, timep);
                /*
                uint c=0;
                uint ci=0;
                for (uint k=0; k < changesp[j]; k++) {
                        if (timep[k] <= tstart) {
                                c++;
                        }
                        else if (timep[k] > tstart && timep[k] <= tend) {
                                ci++;
                        }
                
                        if (timep[k] > tend) break;
                }
                if (semantic == 0) {
                        if (c%2==1 || ci > 0) res[++i] = edgesp[j];
                }
                else if (semantic == 1) {
                        if (c%2==1 && ci == 0) res[++i] = edgesp[j];
                }
                */
                uint *low = lower_bound(timep, timep+changesp[j], tstart+1);
                uint *mid = lower_bound(timep,  timep+changesp[j], tend+1);
                
                //if ( (low-timep)%2 == 1) res[++i] = edgesp[j];
                
                if (semantic == 0) {
                        if ((low-timep)%2==1 || (mid-low) > 0) res[++i] = edgesp[j];
                }
                else if (semantic == 1) {
                        if ((low-timep)%2==1 && (mid-low) == 0) res[++i] = edgesp[j];
                }
                
        }
        
        *res = i;

        delete [] timep;
	delete [] changesp;
	delete [] edgesp;
	delete [] edgetimesizep;
        
        
}


void TGraph::direct_weak(uint v, uint tstart, uint tend, uint *res)  {
        direct_interval(v, tstart, tend, 0, res);
}

void TGraph::direct_strong(uint v, uint tstart, uint tend, uint *res)  {
        direct_interval(v, tstart, tend, 1, res);
}


uint TGraph::snapshot(uint t){
  
	/*
  uint *buffer = new uint [BUFFER];
  
  uint edges=0;
  for(uint v=0; v < nodes; v++) {
	  printf("processing... %.2f%%\r", (float)v/nodes*100);
	  *buffer = 0;
          direct_point(v, t, buffer);
          edges += *buffer;
  }
  
  delete [] buffer;
  
  return edges;
  */
  
    // This code fails... I dont know why :(
  uint *res = new uint [BUFFER];
        
  uint *timep = new uint[BUFFER];
	uint *changesp = new uint[BUFFER];
	uint *edgesp = new uint[BUFFER];
	uint *edgetimesizep = new uint[BUFFER];
        
  uint edges=0;
  for(uint v=0; v < nodes; v++) {
    if (v>=nodes || tgraph[v].neighbors == 0) continue;

  	cc->Decompress(tgraph[v].cchanges, changesp, tgraph[v].neighbors);
  	cc->Decompress(tgraph[v].cedges, edgesp, tgraph[v].neighbors);
  	decodediff(edgesp, tgraph[v].neighbors);
	
  	cc->Decompress(tgraph[v].cedgetimesize, edgetimesizep, tgraph[v].neighbors);
	
    uint i=0;
        
    for(uint j=0; j < tgraph[v].neighbors; j++) {
        decodetime(v, j, edgetimesizep, changesp, timep);
                  
                  //uint c=0;
                  //for (uint k=0; k < changesp[j]; k++) {
                  //        if (timep[k] <= t) {
                  //                c++;
                  //        }
                  //}
                
                  //if (c%2 == 1) res[++i] = edgesp[j];
                                    
        uint *low = lower_bound(timep, timep+changesp[j], t+1);
        if ( (low-timep)%2 == 1) res[++i] = edgesp[j];
    }
        
    *res = i;
    edges += *res;
  }
        
  delete [] timep;
	delete [] changesp;
	delete [] edgesp;
	delete [] edgetimesizep;
  delete [] res;
        
  return edges;
  
}

int TGraph::edge_point(uint v, uint u, uint t){
	 if (v>=nodes || tgraph[v].neighbors == 0) return 0;
        uint *timep = new uint[BLOCKSIZE*tgraph[v].csize_ctime];
	uint *changesp = new uint[BLOCKSIZE+tgraph[v].neighbors];
	uint *edgesp = new uint[BLOCKSIZE+tgraph[v].neighbors];
	uint *edgetimesizep = new uint[BLOCKSIZE+tgraph[v].neighbors];
	
	cc->Decompress(tgraph[v].cchanges, changesp, tgraph[v].neighbors);
	cc->Decompress(tgraph[v].cedges, edgesp, tgraph[v].neighbors);
	decodediff(edgesp, tgraph[v].neighbors);
	
	cc->Decompress(tgraph[v].cedgetimesize, edgetimesizep, tgraph[v].neighbors);
        
        
        uint *p;
        p = (uint *)bsearch(&u, edgesp, tgraph[v].neighbors, sizeof(uint), compare);
        

        int ok = 0;
        if (p == NULL) { return ok;}
        else {

                uint j = (uint)(p - edgesp);
                 decodetime(v, j, edgetimesizep, changesp, timep);
                
      /*          
                 	uint c=0;
                for (uint k=0; k < changesp[j]; k++) {
                       if (timep[k] <= t) {
                               c++;
                       }
                       if (timep[k] > t) {
			       break;
		       }
                }
		if ( c % 2 == 1) ok = 1;
    */
    uint *low = lower_bound(timep, timep+changesp[j], t+1);
    if ( (low-timep)%2 == 1) ok = 1;
    
    
        }
	
	
        
        delete [] timep;
	delete [] changesp;
	delete [] edgesp;
	delete [] edgetimesizep;
        
        
        return ok;
                
        
}

int TGraph::edge_interval(uint v, uint u, uint tstart, uint tend, uint semantic){
	if (v>=nodes || tgraph[v].neighbors == 0) return 0;
	uint *timep = new uint[BLOCKSIZE*tgraph[v].csize_ctime];
	uint *changesp = new uint[BLOCKSIZE+tgraph[v].neighbors];
	uint *edgesp = new uint[BLOCKSIZE+tgraph[v].neighbors];
	uint *edgetimesizep = new uint[BLOCKSIZE+tgraph[v].neighbors];

	cc->Decompress(tgraph[v].cchanges, changesp, tgraph[v].neighbors);
	cc->Decompress(tgraph[v].cedges, edgesp, tgraph[v].neighbors);
	decodediff(edgesp, tgraph[v].neighbors);

	cc->Decompress(tgraph[v].cedgetimesize, edgetimesizep, tgraph[v].neighbors);
        
        
        uint *p;
        p = (uint *)bsearch(&u, edgesp, tgraph[v].neighbors, sizeof(uint), compare);
	
        int ok = 0;
        if (p == NULL) return ok;
        else {
                uint j = (uint)(p - edgesp);
                 decodetime(v, j, edgetimesizep, changesp, timep);
                
                
                 /*
                 	uint c=0, ci=0;
                for (uint k=0; k < changesp[j]; k++) {
                       if (timep[k] <= tstart) {
                                c++;
                       }
                       else if (timep[k] > tstart && timep[k]<= tend) {
                               ci++;
                       }
               
                       if (timep[k] > tend) break;

                }
                
                if (semantic == 0) {
                        if (c%2==1 || ci > 0) ok = 1;
                }
                
                else if (semantic == 1) {
                         if (c%2==1 && ci == 0) ok = 1;
                }*/
                
                uint *low = lower_bound(timep, timep+changesp[j], tstart+1);
                uint *mid = lower_bound(timep,  timep+changesp[j], tend+1);
                
                //if ( (low-timep)%2 == 1) res[++i] = edgesp[j];
                
                if (semantic == 0) {
                        if ((low-timep)%2==1 || (mid-low) > 0) ok=1;
                }
                else if (semantic == 1) {
                        if ((low-timep)%2==1 && (mid-low) == 0) ok=1;
                }
        }
        
        delete [] timep;
	delete [] changesp;
	delete [] edgesp;
	delete [] edgetimesizep;
       
        
        return ok;

}

int TGraph::edge_strong(uint u, uint v, uint tstart, uint tend){
        return edge_interval(u, v, tstart, tend, 1);
}

int TGraph::edge_weak(uint u, uint v, uint tstart, uint tend){
        return edge_interval(u, v, tstart, tend, 0);
}


int TGraph::edge_next(uint v, uint u, uint t){
	if (v>=nodes || tgraph[v].neighbors == 0) return -1;
	uint *timep = new uint[BLOCKSIZE*tgraph[v].csize_ctime];
	
	uint *changesp = new uint[BLOCKSIZE+tgraph[v].neighbors];
	uint *edgesp = new uint[BLOCKSIZE+tgraph[v].neighbors];
	uint *edgetimesizep = new uint[BLOCKSIZE+tgraph[v].neighbors];

	cc->Decompress(tgraph[v].cchanges, changesp, tgraph[v].neighbors);
	cc->Decompress(tgraph[v].cedges, edgesp, tgraph[v].neighbors);
	decodediff(edgesp, tgraph[v].neighbors);

	cc->Decompress(tgraph[v].cedgetimesize, edgetimesizep, tgraph[v].neighbors);
        
        
        uint *p;
        p = (uint *)bsearch(&u, edgesp, tgraph[v].neighbors, sizeof(uint), compare);

        int tnext = -1;
        if (p == NULL) return tnext;
        else {
                uint j = (uint)(p - edgesp);
                 decodetime(v, j, edgetimesizep, changesp, timep);
                
		 uint r;
                uint c=0;
                for (uint k=0; k < changesp[j]; k++) {
			r=k;
                       if (timep[k] <=t) {
                               c++;
                       }
                       if (timep[k] > t) {
                               break;
                       }
                       
                }
		
                if (c%2 == 1) tnext = t;
                else tnext = timep[r];
		
		if (tnext < (int)t) tnext = -1;
        }
        
        delete [] timep;
        
        return tnext;

}

void TGraph::reverse_point(uint v, uint t, uint *res) {
	if (v>=nodes || reverse[v].size == 0) return;

	uint *nodep = new uint[reverse[v].size+BLOCKSIZE];

	decodereverse(v, nodep);

        uint i=0;

        for(uint j=0; j < reverse[v].size; j++) {

        	if (edge_point(nodep[j], v, t)) {
        		res[++i] = nodep[j];
        	}

        }
        *res = i;
        delete [] nodep;
}

void TGraph::reverse_weak(uint v, uint tstart, uint tend, uint *res) {
	if (v>=nodes || reverse[v].size == 0) return;

	uint *nodep = new uint[reverse[v].size+BLOCKSIZE];

	decodereverse(v, nodep);

        uint i=0;

        for(uint j=0; j < reverse[v].size; j++) {
        	if (edge_weak(nodep[j], v, tstart, tend)) {
        		res[++i] = nodep[j];
        	}

        }
        *res = i;
        delete [] nodep;
}

void TGraph::reverse_strong(uint v, uint tstart, uint tend, uint *res) {
	if (v>=nodes || reverse[v].size == 0) return;

	uint *nodep = new uint[reverse[v].size+BLOCKSIZE];

	decodereverse(v, nodep);

        uint i=0;

        for(uint j=0; j < reverse[v].size; j++) {
        	if (edge_strong(nodep[j], v, tstart, tend)) {
        		res[++i] = nodep[j];
        	}

        }
        *res = i;
        delete [] nodep;
}
