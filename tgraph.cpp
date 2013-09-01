#include "tgraphreader.h"
#include "tgraph.h"


#include "cppUtils.h"
#include "debug.h"

#include "arraysort.h"
#include <cstdlib>
#include <cstring>
using namespace std;



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
                cds_utils::saveValue<uint>(out, tgraph[i].neighbors, tgraph[i].size_neighbors);
                cds_utils::saveValue<uint>(out, tgraph[i].csize_neighbors, tgraph[i].size_neighbors);
                cds_utils::saveValue<uint>(out, tgraph[i].changes, tgraph[i].size_neighbors);
                cds_utils::saveValue<uint>(out, tgraph[i].ctime, tgraph[i].csize_time);
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
                tg->tgraph[i].neighbors = cds_utils::loadValue<uint>(in, tg->tgraph[i].size_neighbors);
                tg->tgraph[i].csize_neighbors = cds_utils::loadValue<uint>(in, tg->tgraph[i].size_neighbors);
                tg->tgraph[i].changes = cds_utils::loadValue<uint>(in, tg->tgraph[i].size_neighbors);
                tg->tgraph[i].ctime = cds_utils::loadValue<uint>(in, tg->tgraph[i].csize_time);
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
        
        map< uint, vector<uint> >::iterator it;
        
        uint size_neighbors;
        uint *neighbors;
        uint *csize_neighbors;
        uint csize_time;
        uint csize;
        uint *changest;

        for(uint i=0; i < nodes; i++) {
                if (i%1000==0) fprintf(stderr, "Compressing: %0.2f%%\r", (float)i*100/nodes);

                
                size_neighbors = tgr.tgraph[i].timepoints.size();
                neighbors = new uint[size_neighbors];
                csize_neighbors = new uint[size_neighbors];
                changest = new uint[size_neighbors];
                
                uint j=0;
                csize_time=0;
                for (it=tgr.tgraph[i].timepoints.begin(); it!=tgr.tgraph[i].timepoints.end(); ++it) {
                        encodediff(it->second);
                        
                        csize = cc->Compress(it->second.data(), ccedgebuffer, it->second.size());
                        
                        memcpy(&cctimebuffer[csize_time], ccedgebuffer, csize * sizeof(uint));
                        
                        neighbors[j] = it->first;
                        csize_neighbors[j] = csize_time;
                        changest[j] = it->second.size();
                        j++;
                        
                        csize_time += csize;
                }

                tgraph[i].ctime = new uint [csize_time];
                memcpy(tgraph[i].ctime, cctimebuffer, csize_time * sizeof(uint));

                tgraph[i].changes = changest;
                tgraph[i].csize_time = csize_time;
                tgraph[i].size_neighbors = size_neighbors;
                tgraph[i].neighbors = neighbors;
                tgraph[i].csize_neighbors = csize_neighbors;
        }
        fprintf(stderr, "\n");
        
        
        
}




void TGraph::decodetime(uint u, uint v, uint *res) {
        //if (tgraph[v].changes == 0) return;
         
        cc->Decompress(&tgraph[u].ctime[ tgraph[u].csize_neighbors[v] ], res, tgraph[u].changes[v]);
        decodediff(res, tgraph[u].changes[v]);
}


void TGraph::direct_point(uint v, uint t, uint *res)  {
        if (v>=nodes || tgraph[v].size_neighbors == 0) return;
        uint *timep = new uint[BLOCKSIZE*tgraph[v].size_neighbors];
        uint i=0;
        
        for(uint j=0; j < tgraph[v].size_neighbors; j++) {
                decodetime(v, j, timep);
                uint c=0;
                for (uint k=0; k < tgraph[v].changes[j]; k++) {
                        if (timep[k] <= t) {
                                c++;
                        }
                }
                
                if (c%2 == 1) res[++i] = tgraph[v].neighbors[j];
        }
        
        *res = i;

        delete [] timep;
}

void TGraph::direct_weak(uint v, uint tstart, uint tend, uint *res)  {
        if (v>=nodes || tgraph[v].size_neighbors == 0) return;
        
        uint *timep = new uint[BLOCKSIZE*tgraph[v].size_neighbors];
        uint i=0;
        
        for(uint j=0; j < tgraph[v].size_neighbors; j++) {
                decodetime(v, j, timep);
                
                uint c=0;
                uint ci=0;
                for (uint k=0; k < tgraph[v].changes[j]; k++) {
                        if (timep[k] <= tstart) {
                                c++;
                        }
                        else if (timep[k] > tstart && timep[k] <= tend) {
                                ci++;
                        }
                
                        if (timep[k] > tend) break;
                }
                
                if (c%2==1 || ci > 0) res[++i] = tgraph[v].neighbors[j];
                
        }
        
        *res = i;

        delete [] timep;
        
}


void TGraph::direct_strong(uint v, uint tstart, uint tend, uint *res)  {
        if (v>=nodes || tgraph[v].size_neighbors == 0) return;
        
        uint *timep = new uint[BLOCKSIZE*tgraph[v].size_neighbors];
        
        uint i=0;
        
        for(uint j=0; j < tgraph[v].size_neighbors; j++) {
                decodetime(v, j, timep);

                uint c=0, ci=0;
                for (uint k=0; k < tgraph[v].changes[j]; k++) {
                        for (uint k=0; k < tgraph[v].changes[j]; k++) {
                                if (timep[k] <= tstart) {
                                        c++;
                                }
                                else if (timep[k] > tstart && timep[k] <= tend) {
                                        ci++;
                                }
                
                                if (timep[k] > tend) break;
                        }
                
                        if (c%2==1 && ci == 0) res[++i] = tgraph[v].neighbors[j];

                }
                
        }
        
        *res = i;

        delete [] timep;
}


uint TGraph::snapshot(uint t){
        uint *buffer = new uint [BUFFER];
        
        uint edges=0;
        for(uint v=0; v < nodes; v++) {
                direct_point(v, t, buffer);
                edges += *buffer;
        }
        
        delete [] buffer;
        
        return edges;
}

int TGraph::edge_point(uint u, uint v, uint t){
        if (v>=nodes || tgraph[u].size_neighbors == 0) return 0;
        
        uint *timep = new uint[BLOCKSIZE*tgraph[v].size_neighbors];
        
        uint *p;
        p = (uint *)bsearch(&v, tgraph[u].neighbors, tgraph[u].size_neighbors, sizeof(uint), compare);
        
        int ok = 0;
        if (p == NULL) return ok;
        else {
                uint j = (uint)(p - tgraph[u].neighbors)/sizeof(uint);
                decodetime(v, j, timep);
                uint c=0;
                
                for (uint k=0; k < tgraph[v].changes[j]; k++) {
                       if (timep[k] <= t) {
                               c++;
                       }
                       
                       if (c%2==1) ok = 1;
                }
        }
        
        delete [] timep;
        
        return ok;
                
        
}

int TGraph::edge_weak(uint u, uint v, uint tstart, uint tend){
        if (v>=nodes || tgraph[u].size_neighbors == 0) return 0;
        
        uint *timep = new uint[BLOCKSIZE*tgraph[v].size_neighbors];
        
        uint *p;
        p = (uint *)bsearch(&v, tgraph[u].neighbors, tgraph[u].size_neighbors, sizeof(uint), compare);
        
        int ok = 0;
        if (p == NULL) return ok;
        else {
                uint j = (uint)(p - tgraph[u].neighbors)/sizeof(uint);
                decodetime(v, j, timep);
                
                uint c=0, ci=0;
                for (uint k=0; k < tgraph[v].changes[j]; k++) {
                       if (timep[k] <= tstart) {
                                c++;
                       }
                       else if (timep[k] > tstart && timep[k]<= tend) {
                               ci++;
                       }
               
                       if (timep[k] > tend) break;

                }
                
                if (c%2==1 || ci > 0) ok = 1;
        }
        
        delete [] timep;
        
        return ok;

}

int TGraph::edge_strong(uint u, uint v, uint tstart, uint tend){
        if (v>=nodes || tgraph[u].size_neighbors == 0) return 0;
        
        uint *timep = new uint[BLOCKSIZE*tgraph[v].size_neighbors];
        
        uint *p;
        p = (uint *)bsearch(&v, tgraph[u].neighbors, tgraph[u].size_neighbors, sizeof(uint), compare);
        
        int ok = 0;
        if (p == NULL) return ok;
        else {
                uint j = (uint)(p - tgraph[u].neighbors)/sizeof(uint);
                decodetime(v, j, timep);
                
                uint c=0, ci=0;
                for (uint k=0; k < tgraph[v].changes[j]; k++) {
                       if (timep[k] <= tstart) {
                                c++;
                       }
                       else if (timep[k] > tstart && timep[k]<= tend) {
                               ci++;
                       }
               
                       if (timep[k] > tend) break;

                }
                
                if (c%2==1 && ci == 0) ok = 1;

        }
        
        delete [] timep;
        
        return ok;

}

int TGraph::edge_next(uint u, uint v, uint t){
        if (v>=nodes || tgraph[u].size_neighbors == 0) return -1;
        
        uint *timep = new uint[BLOCKSIZE*tgraph[v].size_neighbors];
        
        uint *p;
        p = (uint *)bsearch(&v, tgraph[u].neighbors, tgraph[u].size_neighbors, sizeof(uint), compare);
        
        int tnext = -1;
        if (p == NULL) return tnext;
        else {
                uint j = (uint)(p - tgraph[u].neighbors)/sizeof(uint);
                decodetime(v, j, timep);
               
                uint c=0;
                for (uint k=0; k < tgraph[v].changes[j]; k++) {
                       if (timep[k] <=t) {
                               c++;
                       }
                       if (timep[k] > t) {
                               if (c%2 == 0) tnext = t;
                               else tnext = timep[k];
                               
                               break;
                       }
                       
                }
        }
        
        delete [] timep;
        
        return tnext;

}
