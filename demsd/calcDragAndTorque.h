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
#include "SDsystem.h"


void outputVector(ofstream &fout, vec3d &vec){
    fout.setf( ios::scientific, ios::floatfield ); // set math types
    int output_precision = 6;
    int output_width = output_precision + 9;
    fout.precision(output_precision); // set precision
    fout << setw(output_width) << vec.x;
    fout << setw(output_width) << vec.y;
    fout << setw(output_width) << vec.z;

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
		cerr << "Usage: stodyn u POSITIONS skipline method" << endl;
		cerr << "Usage: stodyn s POSITIONS skipline method" << endl;
        cerr << "method=1 : Stokesian dynamics without lublication correction" << endl;
		cerr << "method=2 : Stokesian dynamics with lublication correction " << endl;
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
    sd_sys.setFlowType(argv[1][0]);
    sd_sys.method_hydroint = atoi(argv[4]); 
    switch( sd_sys.method_hydroint ){
		case 1:
            // without lubrication corrections
            sd_sys.setLubrication(0);
			break;
		case 2:
            // with lubrication corrections
			sd_sys.setLubrication(1);
			break;
        default:
            cerr << "4th argument\n";
            cerr << "1: without lubrication\n";
            cerr << "2: with lubrication\n";
            exit(1);
	} 
    /* init() is called in this function*/
    //
    sd_sys.setBox(100, 100, 100);
	sd_sys.importCluster(argv[2], atoi(argv[3]));
    sd_sys.initLibStokes();
    sd_sys.setPositionLibStokes();
    sd_sys.setMotionRigidCluster(0,0,0,0,0,0); // (vx,vy,vz,ox,oy,oz)
	sd_sys.setSDIterationMethod();

	sd_sys.solveStokesianDynamics();

    ofstream fout;
    char fout_name[128];
    sprintf(fout_name, "DragTorque_%s.dat", sd_sys.infoString());
    cerr << fout_name << endl;
    fout.open( fout_name );
    
    cerr << fout_name << endl;

    fout << "# N " << sd_sys.np << endl; 
    fout << "# x y z vx vy vz ox oy oz fx fy fz tx ty tz"<< endl;

	for (int i = 0; i < sd_sys.np; i ++){
        vec3d position = sd_sys.getPosition(i);
        vec3d velocity = sd_sys.getVelocity(i);
        vec3d omega = sd_sys.getOmega(i);
        vec3d force = -sd_sys.getForce(i);
        vec3d torque = -sd_sys.getTorque(i);
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
