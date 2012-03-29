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
    char type_of_flow;
    int lubrication;

    int output_precision;
    int output_width;
    char cluster_file[64];
    vector<vec3d> init_cluster;
    vec3d center_of_mass;
    void calcCenterOfMass(vector<vec3d> &p);
    vec3d cl_velocity;
    vec3d cl_omega;
    
public:
	SDsystem();
	~SDsystem();
    void setFlowType(char type_of_flow_);
    void initFlowModel(int num_of_particle_);
    void setSDIterationMethod();
    void setPositionLibStokes();
	void setPositionSD(int, const vec3d &);
    void setLubrication(int lub_);
    void setOutputPrecision(int);
    void setMotionRigidCluster(double vx, double vy, double vz,
                               double ox, double oy, double oz);
    
    vec3d getPosition(int i){
        return vec3d (sd->pos[i*3    ],
                      sd->pos[i*3 + 1],
                      sd->pos[i*3 + 2]);
    }
    vec3d getVelocity(int i){
        return vec3d (velocity[i*3],
                      velocity[i*3 + 1],
                      velocity[i*3 + 2]);
    }
    vec3d getOmega(int i){
        return vec3d (omega[i*3],
                      omega[i*3 + 1],
                      omega[i*3 + 2]);
    }
    vec3d getForce(int i){
        return vec3d (force[i*3],
                      force[i*3 + 1],
                      force[i*3 + 2]);
    }
    vec3d getTorque(int i){
        return vec3d (torque[i*3],
                      torque[i*3 + 1],
                      torque[i*3 + 2]);
    }
    
    void importCluster(char*, int skipline);    
    char typeOfFlow(){return type_of_flow;}
    char* infoString();
    /*************************
     *   Stokesian Dynamics
     *************************/
    struct stokes *sd;
    void calcGrandMovMatrix(double* mov);
    void solveStokesianDynamics();
    /* method_hydroint
     * 0: Free-draining approximation
     * 1: Stokesian dynamics without lubrication
     * 2: Stokesian dynamics with lubrication
     */
    int method_hydroint;
    int np;
    int nm;
    /*
     *
     */
    vec3d *pos;
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
};

#endif
