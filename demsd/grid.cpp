/*
 *  grid.cpp
 *  DEMsd
 *
 *  Created by Ryohei Seto on 09/08/18.
 *  Copyright 2009 Ryohei Seto. All rights reserved.
 *
 */
#include "grid.h"
#include <algorithm>
Grid::Grid(){
    // Particle should be in the box...
    // 
    l_center = 50;
}

Grid::~Grid(){
    if (allocated){
        for (int gx = 0; gx < gx_max; ++gx){
            for (int gy = 0; gy < gy_max ; ++gy){
                delete [] vcell[gx][gy];
                delete [] neighbor_cell[gx][gy];
            }
            delete [] vcell[gx];
            delete [] neighbor_cell[gx];
        }
        delete [] vcell;
    }
}

void Grid::init(unsigned int num_of_particle_,
                double grid_size)
{
    allocated= true;
	num_of_particle = num_of_particle_;
	h = grid_size;
	gx_max = (unsigned int)( 2*l_center/h );
	gy_max = (unsigned int)( 2*l_center/h );
	gz_max = (unsigned int)( 2*l_center/h );
	vcell = new vector<int> ** [gx_max];
	neighbor_cell = new	vector<GridPoint> ** [gx_max];
	for (int gx = 0; gx < gx_max; ++gx){
		vcell[gx] = new vector<int> * [gy_max];
		neighbor_cell[gx] = new	vector<GridPoint> * [gy_max];
		for (int gy = 0; gy < gy_max; ++gy){	
			vcell[gx][gy] = new vector<int> [gz_max];
			neighbor_cell[gx][gy] = new	vector<GridPoint> [gz_max];
		}
	} 
	gl.resize(num_of_particle);
	for (unsigned int gx = 1; gx < gx_max-1; ++gx){
		for (unsigned int gy = 1; gy < gy_max-1; ++gy){	
			for (unsigned int gz = 1; gz < gz_max-1; ++gz){
				GridPoint gp;
				for(gp.x = gx-1; gp.x <= gx+1; ++gp.x){	
					for(gp.y = gy-1; gp.y <= gy+1; ++gp.y){	
						for(gp.z = gz-1; gp.z <= gz+1; ++gp.z){
                            GridPoint gp_tmp = gp;
                            neighbor_cell[gx][gy][gz].push_back( gp_tmp );
						}
					}
				}
			}
		}
	}
}

void Grid::clear_vcell(){
    for (int i = 0; i < num_of_particle; ++i ){
        vcell[gl[i].x][gl[i].y][gl[i].z].clear();
    }
}

void Grid::remake(vector<Particle *> &particle){
    clear_vcell();
	for (int i = 0; i < num_of_particle; ++i ){
        gl[i] = p_to_grid( particle[i]->p ); 
		vcell[gl[i].x][gl[i].y][gl[i].z].push_back(i);
	}
}

GridPoint Grid::p_to_grid(const vec3d &p){
    gp_tmp.x = (unsigned int)((p.x+l_center)/h);
	gp_tmp.y = (unsigned int)((p.y+l_center)/h);
	gp_tmp.z = (unsigned int)((p.z+l_center)/h);
	return gp_tmp;  
}

void Grid::entry(const vec3d &p, int i){
	GridPoint gp = p_to_grid(p);
	vcell[gp.x][gp.y][gp.z].push_back(i);
	gl[i] = gp; 
}

void Grid::get_neighbor_list(const vec3d &p, vector<int> &neighbor){
    gp_tmp = p_to_grid(p);
	foreach(vector<GridPoint>, neighbor_cell[gp_tmp.x][gp_tmp.y][gp_tmp.z], gp){
		foreach(vector<int>, vcell[(*gp).x][(*gp).y][(*gp).z], iter){
			neighbor.push_back(*iter);
		}
	}
}

void Grid::get_neighbor_list_pointer(const vec3d &p, vector< vector<int>* > &p_neighbor){
	p_neighbor.clear();
    gp_tmp = p_to_grid(p);
    foreach(vector<GridPoint>, neighbor_cell[gp_tmp.x][gp_tmp.y][gp_tmp.z], gp){
		p_neighbor.push_back( &(vcell[(*gp).x][(*gp).y][(*gp).z]));
	}
}
