//
//  main.cpp
//  DEMsd
//
//  Created by Ryohei SETO on 12/03/21.
//  Copyright (c) 2012 Ryohei Seto. All rights reserved.
//
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "my_utilities.h"
#include "dem.h"
#include "common.h"
#include "dragforces.h"
#include "ShearRheology.h"
using namespace std;

int main (int argc, char** argv) {
	if ( argc <= 1 ){
		cerr << "Usage: stodyn TYPE ..." << endl;
        cerr << "D: DEM simulation" << endl;
		cerr << "u: calcDragAndTorque in uniform flows" << endl;
		cerr << "s: calcDragAndTorque in shear flows" << endl;
		return 0;
	}
    
    switch (argv[1][0]){
		case 'R':
			shearRheology(argc, argv);
			break;
			
        case 'D':
            /* DEM simulation for an isolated cluster in shear flow.
             * The method is explained in the publification;
             * "Restructuring of colloidal aggregates in shear flow: 
             * Coupling interparticle contact models with Stokesian dynamics".
             */
            if (argc == 2){
                cerr << "usage:" << endl;
                cerr << "stodyn D parameters cluster skip shearrate version" << endl;
                return 0;
            }
            demSimulation(argc, argv);
            break;
        case 'C':
            /*
             * Test for contact model
             */ 
            testContactModel(argc, argv);
            break;
        case 'u':
        case 's':
            /* Drag forces in uniform or shear flow.
             * $ demsd u/s POSITIONS skipline method
             * OUTPUT:
             * 
             */
            calcDragAndTorque(argc, argv);
            break;
        default:
            cerr << "D/u/s" << endl;
    }
	return 0;
}
