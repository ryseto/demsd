//
//  sdsystem.h
//  DEMsd
//
//  Created by Ryohei SETO on 12/03/20.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#ifndef stodyn_sdsystem_h
#define stodyn_sdsystem_h
#include <string>
#include <vector>
#include "stokes.h" //Ichiki-san
#include "ewald-3fts.h" //Ichiki-san
#include "ewald-3fts-matrix.h" //Ichiki-san
#include "ewald-3ft.h" //Ichiki-san
#include "libiter.h" //Ichiki-san
#include "vec3d.h"
using namespace std;

class SDsystem {
private:
    char type_of_flow; // s: shear flow, u: uniform flow
    vec3d cl_velocity; // For rigid clusters
    vec3d cl_omega; // For rigid clusters
public:
    struct stokes *sd; // composition of libstokes (Ichiki)
    int np;
    int nm;
    /* 
     * 
     */
    double * velocity;
    double * omega;
    double * strain_velocity;
    double * force;
    double * torque;
    double * stresslet;
    vec3d cl_force; // force for the rigid cluster
    vec3d cl_torque; // torque for the rigid cluster
    double cl_stresslet[5]; // stresslet for the rigid cluster
    double cl_stresslet_norm;
    
	SDsystem();
	~SDsystem();
    void setFlowType(char type_of_flow_);
    /* prepare_vectors: Whether allocate or not of memories for velocity and force vectors.
     *
     */
    void initFlowModel(int num_of_particle_,
                       int lub_correction_,
                       bool prepare_vectors);
    void setSDIterationMethod();
//    void setPositionLibStokes(vector <vec3d> &_pos);
	void setPositionSD(int, const vec3d &);
//    void setLubrication(int lub_);
    void setOutputPrecision(int);
    void setMotionRigidCluster(double vx, double vy, double vz,
                               double ox, double oy, double oz);
    void importCluster(char*, int skipline);    
    char typeOfFlow(){return type_of_flow;}
    char* infoString();
    void calcGrandMovMatrix(double* mov_mtrx);
    void solveStokesianDynamics();
    inline vec3d Position(int i){
        return vec3d (sd->pos[i*3    ],
                      sd->pos[i*3 + 1],
                      sd->pos[i*3 + 2]);
    }
    inline vec3d Velocity(int i){
        return vec3d (velocity[i*3],
                      velocity[i*3 + 1],
                      velocity[i*3 + 2]);
    }
    inline vec3d Omega(int i){
        return vec3d (omega[i*3],
                      omega[i*3 + 1],
                      omega[i*3 + 2]);
    }
    inline vec3d Force(int i){
        return vec3d (force[i*3],
                      force[i*3 + 1],
                      force[i*3 + 2]);
    }
    inline vec3d Torque(int i){
        return vec3d (torque[i*3],
                      torque[i*3 + 1],
                      torque[i*3 + 2]);
    }
};
#endif
