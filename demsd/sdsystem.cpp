//
//  sdsystem.cpp
//  DEMsd
//
//  Created by Ryohei SETO on 12/03/20.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#include "sdsystem.h"
#include <iostream>
#include <iomanip>
#include <fstream>

SDsystem::SDsystem(){
    np = -1;
	lubrication = -1;
//	pos = NULL;
	velocity = NULL;
	omega = NULL;
	strain_velocity = NULL;
	stresslet = NULL;
	force = NULL;
	torque = NULL;
}

SDsystem::~SDsystem(){
}

char* SDsystem::infoString(){
    char *sd_info;
    sd_info = new char [32];
    sprintf(sd_info, "flow_%c_lub_%d_%s", type_of_flow, lubrication, cluster_file);
    return sd_info;
}


void SDsystem::setLubrication(int lub_){
    lubrication = lub_;
}

void SDsystem::setOutputPrecision(int value){
    output_precision = value;
    output_width = value + 9;
}


/* Import the configuration (x,y,z) of particles.
 * The center of mass is set to (0,0,0)
 * INPUT 
 * importfilename: File name
 * skipline: Line number of the header
 */
void SDsystem::importCluster(char* cluster_file_, int skipline){
	sprintf( cluster_file, "%s", cluster_file_);
	string s_filename = cluster_file;
	int i_backslash = s_filename.find_last_of( "/") + 1;
	int i_extention = s_filename.find( ".dat" );	
	sprintf(cluster_file, "%s",
			(s_filename.substr(i_backslash,i_extention-i_backslash)).c_str());
	ifstream fin;
	fin.open( cluster_file_ );
	double x, y, z;
	char buf[1000];
	for (int i = 0; i< skipline; i++){
		fin.getline( buf, 1000 );
	}
    
	do{
		fin >> x >> y >> z;
        if( fin.fail() )
            break; 
        vec3d new_p(x, y, z);
		init_cluster.push_back(new_p);
	} while (!fin.eof());
    np = init_cluster.size();
    /*
     * Set center-of-mass to (0,0,0)
     */
    calcCenterOfMass(init_cluster);
    for (int i=0; i < np; i++){
        init_cluster[i] -= center_of_mass;
    }
}

void SDsystem::calcCenterOfMass(vector<vec3d> &pos_vec){
    center_of_mass.set(0,0,0);
    for (int i=0; i < np; i ++){
        center_of_mass += pos_vec[i];
    }
    center_of_mass *= 1.0/np;
}

void SDsystem::setFlowType(char type_of_flow_){
    /* s : shear flow
     * u : uniform flow
     */
    type_of_flow = type_of_flow_;
}

/* Initialize libstokes
 * Basic parameters for the libstorks
 * and pos[] are set from init_aggregate.
 */
void SDsystem::initFlowModel(int num_of_particle_){
    if (num_of_particle_ > 0)
        np = num_of_particle_;
    nm = np ;	/* nm : number of mobile particles  */
    
    if (method_hydroint != 0){
        sd = stokes_init(); // ---> libstokes
        sd->twobody_lub = lubrication;
        // For shear flows,"FTS version" is required.
        sd->version = 2; /* 0 = F, 1 = FT, 2 = FTS  */
        sd->periodic = 0; /* 0 = non periodic, 1 = periodic */
        stokes_set_np(sd, np, nm);
        if (type_of_flow == 'u'){
            /*
             * The velocity of uniform flow is alywas (1,0,0)
             */
            stokes_set_Ui(sd, 1.0, 0.0, 0.0);
        } else if (type_of_flow == 's'){
            /*
             * The shear rate of shear flow is always 1.0
             */
            //sd->Ui[0] = -1.0*lz/2 ;
            sd->Ui[0] = 0;
            sd->Oi[1] = 1.0/2;
            sd->Ei[2] = 1.0/2;
        } else {
            cerr << "Type of flow should be given. (shear or uniform)" << endl;
        }
    }
    try{
//        pos = new vec3d [np];        
        velocity = new double [ np*3 ];
        omega = new double [ np*3 ];
        strain_velocity = new double [ np*5 ];
        force = new double [ np*3 ];
        torque = new double [ np*3 ];
        stresslet = new double [ np*5 ];
    } catch (bad_alloc &){
        cerr << "bad_alloc at System::init()" << endl;
    }
}

void SDsystem::setPositionLibStokes(){
    /* For SD calculation,
     * sd->pos[*] needs to be set.
     */
    double scale = 1.0;
    for (int i = 0; i < np; i++){
        setPositionSD(i, scale*init_cluster[i]);
    }
}

void SDsystem::setPositionSD(int i, const vec3d &position){
    // Positions are set in libstokes.
	int j = 3*i;
//    pos[i] = position;
	sd->pos[j] = position.x;
	sd->pos[j+1] = position.y;
	sd->pos[j+2] = position.z;
}

void SDsystem::setMotionRigidCluster(double vx, double vy, double vz,
                                     double ox, double oy, double oz){                                   
    cl_velocity.set(vx, vy, vz);
    cl_omega.set(ox, oy, oz);
    
    for (int i = 0; i < np; i++){
        vec3d pos(sd->pos[3*i], sd->pos[3*i+1], sd->pos[3*i+2]);
		vec3d v = cross(cl_omega, pos) + cl_velocity;
        //vec3d v = cross_vec_array(cl_omega, &(sd->pos[3*i]))+ cl_velocity;
		velocity[i*3] = v.x;
		velocity[i*3+1] = v.y;
		velocity[i*3+2] = v.z;
		omega[i*3] = cl_omega.x;
		omega[i*3+1] = cl_omega.y;
		omega[i*3+2] = cl_omega.z;
	}
}

void SDsystem::setSDIterationMethod(){
    /* 
     * Set parameters for the Ichiki's codes.
     */ 
    /* setup ewald summation */
	double ewald_tr = 60.25;
	double xi = xi_by_tratio(sd, ewald_tr);
	double ewald_eps = 1.0e-12;
	stokes_set_xi(sd, xi, ewald_eps);
	//fprintf (stdout, "xi = %f\n", xi);
	//sys->lubmin = 2.0000000001;
	
	sd->lubmin2 = 4.0000000001;
	sd->lubmax = 4.00;
	
	//sd->lubmax = 2.0;
	/* set iter param
	 * INPUT
	 *   solver : string indicating the solver
	 *            sta, sta2, gpb, otmk , or gmres (default)
	 *   eps and log10_eps
	 *   max (and restart)
	 *   debug = 0 : no debug info
	 *         = 1 : iteration numbs and residue
	 *   out   : FILE * to output debug info.
	 */
	stokes_set_iter(sd, "gmres", 2000, 20, 1.0e-6, 0, stdout);
	//stokes_set_iter(sd, "sta", 2000, 20, 1.0e-6, 0, stdout);
}


// (U,O,S) = mov*(F,T,E)
void SDsystem::calcGrandMovMatrix(double* mov){
    calc_mob_3fts_matrix(sd, mov);
}

void SDsystem::solveStokesianDynamics(){
    if (sd->twobody_lub == 0){
        solve_res_3fts(sd, velocity, omega, strain_velocity,
                       force, torque, stresslet);
    } else {
        solve_res_lub_3fts(sd,
                           velocity, omega, strain_velocity,
                           force, torque, stresslet);
    }    
}










