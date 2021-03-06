#ifndef TGRAPH_H_
#define TGRAPH_H_ 

#include <fstream>
#include <sys/types.h>
#include <cstring>
#include <vector>
#include "tgraphreader.h"
#include "coding_policy.h"
#include "debug.h"

#define BLOCKSIZE 128
#define BUFFER 4*67108864 //4*256 megabytes, 67108864 integers

#define CP_S9 "s9"
#define CP_S16 "s16"
#define CP_VBYTE "vbyte"
#define CP_RICE "rice"
#define CP_PFOR "pfor:32:s16:32"

using namespace std;

enum TypeGraph {
	kInterval, kGrowth, kPoint
};

struct opts {
	char *outfile;
	char cp[100]; //coding policy
	enum TypeGraph typegraph;
};

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
	uint *cedges; //compressed list of edges (size: neighbors)
	uint csize_cedges;
	
	uint *cchanges; //compressed number of changes by edge (size: neighbors)
	uint csize_cchanges;
	
	// DELTA ENCODE
	uint *ctime;  //compressed timepoints (size: unknown, upper bound is csize_ctime*32 integers )
	uint csize_ctime;
	
	uint *cedgetimesize; //size of each compressed list of time edges (size: neighbors)
	uint csize_cedgetimesize;
	
};

class TGraph {
public:
        uint nodes;
        uint edges;
        uint changes;
        uint maxtime;
        
        struct opts opts;
	
        TGraphEventList* tgraph;

        TGraphReverse *reverse;

        CodingPolicy *cc;
        
        void set_opts(struct opts o) {
        	opts = o;
        }
	        void loadpolicy() {
                cc = new CodingPolicy(CodingPolicy::kPosition);
				cc->LoadPolicy(opts.cp);
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

    size_t change_point(uint t);
    size_t change_interval(uint ts, uint te);
    size_t actived_point(uint t);
    size_t actived_interval(uint ts, uint te);
    size_t deactived_point(uint t);
    size_t deactived_interval(uint ts, uint te);


    // point contact graphs
    void direct_interval_pg(uint node, uint tstart, uint tend, uint *res) ;
    //void reverse_interval_pg(uint node, uint tstart, uint tend, uint *res) ; //already captured by edge_interval_pg

    int edge_interval_pg(uint u, uint v, uint tstart, uint tend) ;
    int edge_next_pg(uint u, uint v, uint t);
    size_t snapshot_pg(uint t);
    //size_t change_interval_pg(uint ts, uint te); //do not need an upgrade
    size_t actived_interval_pg(uint ts, uint te);
    size_t deactived_interval_pg(uint ts, uint te);
};


void decodediff(uint *k, uint size);
void encodediff(vector<uint> &t);

#endif
