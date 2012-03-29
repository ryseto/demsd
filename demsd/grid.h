/*
 *  grid.h
 *  DEMsd
 *
 *  Created by Ryohei Seto on 09/08/18.
 *  Copyright 2009 Ryohei Seto. All rights reserved.
 *
 */
#ifndef grid_h
#define grid_h 1
#include <vector>
#include <map>
#include "vec3d.h"
#include "particle.h"
using namespace std;
const int max_number_point =20000;
class Particle;

struct GridPoint {
	unsigned int x;
	unsigned int y;
	unsigned int z;
};

class Grid{
	double h;
    double l_center;
	vector<GridPoint> gl;
	vector<int > ***vcell;
	vector<GridPoint> ***neighbor_cell;
	unsigned int  num_of_particle;
	GridPoint gp_tmp;
	unsigned int gx_max;
	unsigned int gy_max;
	unsigned int gz_max;
    bool allocated;
    void clear_vcell();
public:
	Grid();
	~Grid();
    
	void init(unsigned int num_of_particle , double grid_size);
	void remake(vector<Particle *> &particle);
	void entry(const vec3d &p, int i);
	GridPoint p_to_grid(const vec3d &p);
	inline vector<int> *particle_in_cell(int x, 
                                         int y,
                                         int z){
        return &vcell[x][y][z];
    }
	void reset();
	inline int val_gx_max(){ return gx_max; }
	inline int val_gy_max(){ return gy_max; }
	inline int val_gz_max(){ return gz_max; }
	inline int min_gy(int gy){ return ( gy >= 0 ? gy : 0 ) ; }
	inline int max_gy(int gy){ return ( gy <= gy_max ? gy : gy_max ); }
	int gy(int i){ return gl[i].y;}
	int size(int x, int y, int z) { return vcell[x][y][z].size(); }
	void gl_resize(int n){ gl.resize(n); }
	void get_neighbor_list(const vec3d &p, vector<int> &neighbor);
	void get_neighbor_list_pointer(const vec3d &p, vector< vector<int>* > &p_neighbor);
	
};
#endif
