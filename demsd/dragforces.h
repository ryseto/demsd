//
//  calcDragAndTorque.h
//  DEMsd
//
//  Created by Ryohei SETO on 12/03/20.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//

#ifndef stodyn_calcDragAndTorque_h
#define stodyn_calcDragAndTorque_h
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "sdsystem.h"

void outputVector(ofstream &fout, vec3d &vec){
    fout.setf( ios::scientific, ios::floatfield ); // set math types
    int output_precision = 6;
    int output_width = output_precision + 9;
    fout.precision(output_precision); // set precision
    fout << setw(output_width) << vec.x;
    fout << setw(output_width) << vec.y;
    fout << setw(output_width) << vec.z;
}


vec3d calcCenterOfMass(vector<vec3d> &pos_vec){
    vec3d center_of_mass(0,0,0);
    int np = pos_vec.size();
    for (int i=0; i < np; i ++){
        center_of_mass += pos_vec[i];
    }
    center_of_mass *= 1.0/np;
    return center_of_mass;
}


/* Import the configuration (x,y,z) of particles.
 * The center of mass is set to (0,0,0)
 * INPUT 
 * importfilename: File name
 * skipline: Line number of the header
 */
void importCluster(char* cluster_file_, int skipline, vector <vec3d> &pos){
    //	sprintf( cluster_file, "%s", cluster_file_);
    //	string s_filename = cluster_file_;
	//int i_backslash = s_filename.find_last_of( "/") + 1;
    //	int i_extention = s_filename.find( ".dat" );	
    //	sprintf(cluster_file, "%s",
    //			(s_filename.substr(i_backslash,i_extention-i_backslash)).c_str());

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
        vec3d new_pos(x,y,z);
        pos.push_back(new_pos);
	} while (!fin.eof());
    /*
     * Set center-of-mass to (0,0,0)
     */
    vec3d com = calcCenterOfMass(pos);
    for (int i=0; i < pos.size() ; i++){
        pos[i] -= com;
    }
}


void calcDragAndTorque(int argc, char** argv){    
    SDsystem sd_sys;
	/* Shear-rate or typical velocity is set one.
	 *
	 */
	if (argc != 5){
		/* POSITIONS indicates the path to a file including (x,y,z) of particles
		 *
		 */
		cerr << "Usage: stodyn u POSITIONS skipline lub" << endl;
		cerr << "Usage: stodyn s POSITIONS skipline lub" << endl;
        cerr << "lub=0 : Stokesian dynamics without lublication correction" << endl;
		cerr << "lub=1 : Stokesian dynamics with lublication correction " << endl;
        cerr << "This calculation will be conducted by dimensionless variables." << endl;
        cerr << "The unit of length is given by the radius of particle." << endl;
        cerr << "The unit of velocity is given by the velocity of uniform flow." << endl;
        cerr << "Or, the unit of velocity is given by the product of the shear rate and the unit of length." << endl;
        cerr << "Here, U=1.0 is used for uniform flow." << endl;
        cerr << "G=1.0 is used for shear flow." << endl;
        cerr << "F=U is expected for one-body solution." << endl;
        cerr << "===================" << endl;
        cerr << "For rescaling:" << endl;
        cerr << "* Uniform flow. " << endl;
        cerr << " a and U are given." << endl;
        cerr << "L0 = a, U0 = U, F0 = 6*M_PI*eta*L0*U0" << endl;
        cerr << "* Shear flow. " << endl;
        cerr << " a and G are given." << endl;
        cerr << "L0 = a, U0 = G*a, F0 = 6*M_PI*eta*L0*U0" << endl;
        cerr << "-------------------" << endl;
        cerr << "F ---> F*F0 " << endl;
        cerr << "T ---> T*F0*L0 " << endl;
        cerr << "-------------------" << endl;
		return;
	}
    vector <vec3d> pos;
    char type_of_flow = argv[1][0];
    char *cluster_file = argv[2];
    int skip_line =  atoi(argv[3]); 
    int lub_correction = atoi(argv[4]);
    sd_sys.setFlowType(type_of_flow);
    /* init() is called in this function*/
	importCluster(cluster_file, skip_line, pos);
    int np = pos.size();
    sd_sys.initFlowModel(np, lub_correction, true);

    for (int i = 0; i < np; i++){
        sd_sys.setPositionSD(i, pos[i]);
    }
    sd_sys.setMotionRigidCluster(0,0,0,0,0,0); // (vx,vy,vz,ox,oy,oz)
	sd_sys.setSDIterationMethod();
	sd_sys.solveStokesianDynamics();

    ofstream fout;
    char fout_name[128];
    sprintf(fout_name, "DragTorque_%s.dat", sd_sys.infoString());
    fout.open( fout_name );

    fout << "# N " << sd_sys.np << endl; 
    fout << "# x y z vx vy vz ox oy oz fx fy fz tx ty tz"<< endl;
	for (int i = 0; i < sd_sys.np; i ++){
        vec3d position = sd_sys.Position(i);
        vec3d velocity = sd_sys.Velocity(i);
        vec3d omega = sd_sys.Omega(i);
        vec3d force = -sd_sys.Force(i);
        vec3d torque = -sd_sys.Torque(i);
        outputVector(fout, position);
        outputVector(fout, velocity);
        outputVector(fout, omega);
        outputVector(fout, force);
        outputVector(fout, torque);
		fout << endl;
    }
    fout.close();
}

#endif
