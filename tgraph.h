#ifndef TGRAPH_H_
#define TGRAPH_H_ 

#include <fstream>
#include <sys/types.h>
#include "tgraphreader.h"
#include "coding_policy.h"

#define BLOCKSIZE 128
#define BUFFER 67108864 //256 megabytes

#define CODINGPOLICY "rice"

class TGraphEventList {
public:

        //size of each neighbor
        uint size_neighbors; //number of neighbors
        uint *neighbors; //list of edges
        uint *csize_neighbors; // compressed size of each neighbor
        
        uint *changes; //changes for each neighbor
        
        //array of compressed time points
        uint csize_time;
        uint *ctime;
};

class TGraph {
public:
        uint nodes;
        uint edges;
        uint changes;
        uint maxtime;
        
        TGraphEventList* tgraph;
        CodingPolicy *cc;
        
        void loadpolicy() {
                cc = new CodingPolicy(CodingPolicy::kPosition);
                cc->LoadPolicy(CODINGPOLICY);
        }
        
        
        void save(ofstream &out);
        
        static TGraph* load(ifstream &in);
        
        void create(TGraphReader &tgr);
        
        void decodetime(uint u, uint v, uint *res);
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
	
        
	/*
        void reverse_point(uint node, uint t, uint *res) const;
	void reverse_weak(uint node, uint tstart, uint tend, uint *res) const;
	void reverse_strong(uint node, uint tstart, uint tend, uint *res) const;
        */
};


void decodediff(uint *k, uint size);
void encodediff(vector<uint> &t);

#endif