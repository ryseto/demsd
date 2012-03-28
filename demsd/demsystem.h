//
//  demsystem.h
//  DEMsd
//
//  Created by SETO Ryohei on 31/07/2011.
//  Copyright 2011 Ryohei Seto. All rights reserved.
//

#ifndef stodyn_demsystem_h
#define stodyn_demsystem_h 1
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "sdsystem.h"
#include "vec3d.h"
#include "bond.h"
#include "particle.h"
#include "grid.h"
#include "quaternion.h"
#include "common.h"
using namespace std;

class Particle;
class SDsystem;
class Bond;
class Grid;

class DEMsystem {
private:
    SDsystem *sd_sys;
    double *mov; // Mobility matrix of SD
    double *vos; // velocity, omega, stresslet
    double *fte; // force, torque, rate-of-strain
//    vec3d *pos; 
    vec3d *pos_init; 
    bool mov_allocated;
    ////////////////////////////////////////////    
	double dist_generate;
//    int count_MD_steps;
    int count_SD_calc;
    vec3d total_torque;
	vec3d total_force;
    double total_stresslet[5];
    double norm_total_stresslet;
    
    int cnt_grad_method;
    /* Simulation box */
//    double lx;
//    double ly;
//    double lz;
    /* 
     * !!!BUG!!!
     * Simulation box 
     * For the moment, 
     * L_dem != L_sd makes conflict.
     * 
     */    
    double L_dem;
    double L_sd;
    int n11;
    /* DEM */
    bool firsttime_dem;
    bool restructuring;
	char parameters_file[128];
	char parameters[32];
    char init_config_file[128];
    string bond0_file;
    string bond1_file;
    string shear_process_file;
    string rootdir;
    
    double radius_of_gyration;
    double init_radius_of_gyration;

    /* Local strains */
    double r_min;  
    double r_max;  
    double bending_angle_max;
    double torsional_angle_max;
    double sliding_disp_max;
    
    double ave_bondforce;
    double ave_force;
    double max_force;
    double max_velocity;
    double max_ang_velocity; 
    double max_x;
    double min_x;
    
    double force_trac_max;
    double force_comp_max;
    double force_slid_max;
    double moment_bend_max;
    double moment_tors_max;
    vec3d f_tmp;
    vec3d t_tmp;
    vec3d v;
    vec3d o;
    double e[5];
    ofstream fout_dem;
    ofstream fout_dem_log;
    ofstream fout_conf;
    ofstream fout_conf_def;
    ofstream fout_trace;
    ofstream fout_data;
    
    void shiftCenterOfMass(vector<vec3d> &p);
    void set_FDA();
    double evaluateObjFunction(quaternion & q_, vec3d *po);
    void readParameterKey(const string &codeword, 
						  const string &value);
    void readShearProcessKey(const string &codeword,                               
                             const string &value);

public:
	DEMsystem(SDsystem &sd_sys_);
	~DEMsystem();

    int np; // number of particle
    vector<Particle *> particle;
    vector<Bond *> bond;
    ConTable *ct;
    Grid *grid;
    string shear_process;
    double shearrate; // = GammaDot (not real shear-rate)
    double shearrate_min;
    double shearrate_max;
    double incrementratio_shearrate;
    int stepnumber_shearrate;
    double simu_time; // unit of time is [F0/6 pi eta a^2] second.
    double h_stepsize; // unit of this is [1/GammaDot]

    double shearstrain;
    
    double shearstrain_step;

    double critical_deformation_SD;
    bool val_restructuring(){
        return restructuring;
    }
    double interval_strain_output;
    
    void setStepSize();
    void set_force_trac_max(double val){
        force_trac_max = val;
    }

    void setParameterFileDEM(char *parameters_file_);
    void readParameterFileDEM();
    void readParameterShearProcess();
    void readParameterBond();

    void checkFailure();
    void pos_from_DEM_to_SD();
    void setVersion(char *);
    void calcInterParticleForce();
    void getMovMatrix();
    void moveOverdampedMotionSDMov();
    void calcStresslet();
    void revertStresslet();
    void calcVelocityOmega();
    
    void initDEM();
    void importCluster(char*, int skipline);    

    void makeInitialBond(double generation_distance);
    void makeNeighbor();

    void checkState();
    void generateBond();
    void resetDeformation();
    void calcTotalDeformation();

	void initGrid();
    void setBox(double lx_, double ly_, double lz_);
	void setDEMParameters();
    void freeDrainingApproximation();
    
    void writeHeader(ofstream &out);
    void output(){
   
    }
	void openOutputFileDEM();
    void outputLogDEM();
    void outputYaplotDEM();
    void outputData();
    void outputConfiguration();
    void outputDeformationConf();

    void setSpringConstant(double &kn,
                           double &ks,
                           double &kb,
                           double &kt);

    void regeneration_onebyone();
    void rupture_onebyone();

    void setInitialMotion();
    void timeEvolution();
    void rupture();
    void shiftClusterToCenter();


   // void estimateClusterRotation();

    void estimateClusterRotation();
    double findObjMinimum(quaternion &q_try, vec3d *po, 
                          double delta_init, double conv_diff);

    void calcLocalStrains();
    void calcGyrationRadius();
    void calcTotalFTS();
    

    BondParameter bond0;
    BondParameter bond1;
    char version[3];
    bool initialprocess;
    vec3d pos_center_of_mass;
    double rot_ang;
    double step_deformation;
    double total_deformation;
    int init_continuity; // 1: continuous, 0: discontinuous
    quaternion q_rot;
//    quaternion q_rot_old;
    quaternion q_rot_total;
    quaternion q_rot_total_step;
    
    
    //double diff_q_rot_0;
//    vec3d diff_q_rot_q;

       double calc_rotation_time;

    int n_bond;
    int rup_normal;
	int rup_shear;
	int rup_bend;
	int rup_torsion;
    vector<int> regeneration_bond;
    vector<int> rupture_bond;
	int counterRegenerate;
	int counterBreak;
    double cp_f_n_max;
	double cp_f_s_max;
	double cp_m_b_max;
	double cp_m_t_max;
    double robust_bond_compression;

//    double lx0;
//    double ly0;
//    double lz0;
    double sq_dist_generate;


    };
#endif



