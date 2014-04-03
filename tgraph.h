#ifndef TGRAPH_H_
#define TGRAPH_H_ 

#include <fstream>
#include <sys/types.h>
#include "tgraphreader.h"
#include "coding_policy.h"
#include "debug.h"

#define BLOCKSIZE 128
#define BUFFER 67108864 //256 megabytes

#define CP_S9 "s9"
#define CP_S16 "s16"
#define CP_VBYTE "vbyte"
#define CP_RICE "rice"
#define CP_PFOR "pfor:32:s16:32"

using namespace std;

enum CP_FORMAT {
	S9, S16, VBYTE, RICE, PFOR,
};

class TGraphReader;

class TGraphReverse {
public:
	uint size; //out degree, number of elements
	uint csize; //compressed size
	uint *clist; //pointer to compressed list
};



class TGraphEventList {
public:

        uint neighbors; //number of edges
	
	// DELTA ENCODE
	uint *cedges; //compressed list of edges
	uint csize_cedges;
	
	uint *cchanges; //compressed number of changes by edge
	uint csize_cchanges;
	
	// DELTA ENCODE
	uint *ctime;  //compressed timepoints
	uint csize_ctime;
	
	uint *cedgetimesize; //size of each compressed list of time edges
	uint csize_cedgetimesize;
	
};

class TGraph {
public:
        uint nodes;
        uint edges;
        uint changes;
        uint maxtime;
        
	enum CP_FORMAT cp;
	
        TGraphEventList* tgraph;

        TGraphReverse *reverse;

        CodingPolicy *cc;
        
	void set_policy(enum CP_FORMAT c) {
		cp = c;
	}
	
        void loadpolicy() {
                cc = new CodingPolicy(CodingPolicy::kPosition);
		
		switch(cp) {
			case S9: cc->LoadPolicy(CP_S9); INFO("Using S9"); break;
			case S16: cc->LoadPolicy(CP_S16); INFO("Using S16"); break;
			case VBYTE: cc->LoadPolicy(CP_VBYTE); INFO("Using VBYTE"); break;
			case RICE: cc->LoadPolicy(CP_RICE); INFO("Using RICE"); break;
			case PFOR: cc->LoadPolicy(CP_PFOR); INFO("Using PFOR"); break;
			default: cc->LoadPolicy(CP_PFOR); INFO("Using S9"); break;
		}
		
                
        }

        
        
        void save(ofstream &out);
        
        static TGraph* load(ifstream &in);
        
        void create(TGraphReader &tgr);
        
        void decodetime(uint u, uint v, uint *edgetimesize, uint *changes, uint *res);

        void decodereverse(uint v, uint *res);


        int isEdgeActive(uint v, uint j, uint t, uint *timep);
        
        uint snapshot(uint t);
        
        int edge_point(uint u, uint v, uint t);
        int edge_weak(uint u, uint v, uint tstart, uint tend);
        int edge_strong(uint u, uint v, uint tstart, uint tend);
        int edge_interval(uint u, uint v, uint tstart, uint tend, uint semantic);
        
        int edge_next(uint u, uint v, uint t);
        
	void direct_point(uint node, uint t, uint *res) ;
	void direct_weak(uint node, uint tstart, uint tend, uint *res) ;
	void direct_strong(uint node, uint tstart, uint tend, uint *res) ;	
	
	void direct_interval(uint node, uint tstart, uint tend, uint semantic, uint *res);
	
        

        void reverse_point(uint node, uint t, uint *res);
	void reverse_weak(uint node, uint tstart, uint tend, uint *res);
	void reverse_strong(uint node, uint tstart, uint tend, uint *res);

};


void decodediff(uint *k, uint size);
void encodediff(vector<uint> &t);

#endif
