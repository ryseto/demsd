//
//  shearRheology.h
//  demsd
//
//  Created by Ryohei Seto on 11/7/12.
//
//

#ifndef demsd_shearRheology_h
#define demsd_shearRheology_h
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include "stokes.h" //Ichiki-san
#include "ewald-3fts.h" //Ichiki-san
#include "ewald-3fts-matrix.h" //Ichiki-san
#include "ewald-3ft.h" //Ichiki-san
#include "libiter.h" //Ichiki-san
#include "lub-matrix.h"
#include "f.h"
#include "vec3d.h"
void shearRheology(int argc, char** argv){
	struct stokes *sd; // composition of libstokes (Ichiki)
	int np = 10;
	int nm = 10;
	sd = stokes_init(); // ---> libstokes

	sd->version = 2; /* 0 = F, 1 = FT, 2 = FTS  */
    sd->periodic = 1; /* 0 = non periodic, 1 = periodic */

    stokes_set_np(sd, np, nm);
	//	sd->twobody_nmax = 100;
	sd->twobody_lub = 1;
	/*
     * Set parameters for the Ichiki's codes.
     */
    /* setup ewald summation */
	//fprintf (stdout, "xi = %f\n", xi);
	//sys->lubmin = 2.0000000001;
	sd->lx = 20;
	sd->ly = 20;
	sd->lz = 20;

	
	double ewald_tr = 60.25;
	double xi = xi_by_tratio(sd, ewald_tr);
	double ewald_eps = 1.0e-12;
	stokes_set_xi(sd, xi, ewald_eps);
	
	//fprintf (stdout, "xi = %f\n", xi);
	//sys->lubmin = 2.0000000001;

	sd->lubmin2 = 4.0000000001;
	sd->lubmax = 4.0000000001;
	
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
	/////////////////////////

	//sd->twobody_lub = lub_correction_;
    // For shear flows,"FTS version" is required.
  

	sd->Ui[0] = 0;
	sd->Oi[1] = 1.0/2;
	sd->Ei[2] = 1.0/2;
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

	vector <vec3d> pos;
	double *a;
	a = new double [10];
	for (int i=0;i<10; i++){
		if (i < 5){
			a[i] = 1;
		} else {
			a[i] = 1.5;
		}
	}
	
	for (int i=0;i<100; i++){
		bool overlap;
		vec3d tmp;
		do{
			tmp.set(19*drand48()+1  ,
					19*drand48()+1 ,
					19*drand48()+1 );
			overlap = false;
			for (int j=0; j < pos.size(); j++){
				vec3d dr = tmp - pos[j];
				double ai_aj = a[i] + a[j];
				if ( dr.sq_norm() < ai_aj*ai_aj){
					overlap = true;
					break;
				}
			}
		}while (overlap == true);
		pos.push_back(tmp);
	}
	
	cerr << "init done" << endl;
	for (int i=0;i<10; i++){
		int j = 3*i;
		sd->pos[j] = pos[i].x -2.5;
		sd->pos[j+1] = pos[i].y-2.5;
		sd->pos[j+2] = pos[i].z-2.5;
		
	}
	stokes_set_radius(sd, a);

	double force[3*np];
	double torque[3*np];
	
	double velocity[3*np];
	double omega[3*np];
	double rate_of_strain[5*np];
	double stresslet[5*np];
	for (int i=0; i < 3*np; i++){
		force[i] = 0;
		torque[i] = 0;
	}
	for (int i=0; i < 5*np; i++){
		rate_of_strain[i] = 0;
	}
	
	//	solve_mob_3ft(sd, force, torque, veolocity, omega);
	solve_mob_3fts(sd, force, torque, rate_of_strain, velocity, omega, stresslet);

	for (int i=0;i<10; i++){
		int j = 3*i;
		cout << "r " << sd->a[i] << endl;
		cout << "c " << ' ';
		cout << sd->pos[j] << ' ';
		cout << sd->pos[j+1] << ' ';
		cout << sd->pos[j+2] << endl;
	}
	
	for (int i=0;i<10; i++){
		int j = 3*i;
		cout << "@ 3" << endl;
		cout << "l" << ' ';
		cout << sd->pos[j] << ' ';
		cout << sd->pos[j+1] << ' ';
		cout << sd->pos[j+2] << ' ';
		cout << sd->pos[j] + velocity[j]<< ' ';
		cout << sd->pos[j+1] + velocity[j+1] << ' ';
		cout << sd->pos[j+2] + velocity[j+2] << endl;
	}

	
}


#endif
