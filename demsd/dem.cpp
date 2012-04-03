//
//  dem.cpp
//  demsd
//
//  Created by Ryohei SETO on 4/2/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "dem.h"

void demSimulation(int argc, char** argv){
    SDsystem sd_sys;    
    DEMsystem dem(sd_sys);
	dem.setParameterFileDEM(argv[2]);   
    /* Import positions of the cluster: (file, skip lines)*/
	dem.importCluster(argv[3], atoi(argv[4]));
    /* You can add arbitary string for name of outpuf file. */
    dem.setVersion(argv[5]);
    /* Read parameter file */
    dem.readParameterFileDEM();
    dem.readParameterBond();
    dem.readParameterShearProcess();    
	/* Setup simulation */
    sd_sys.setFlowType('s'); // set shear flow
	dem.initDEM();
    sd_sys.initFlowModel(dem.np);
    /* Main simulation */
	if ( dem.shear_process == "stepwise"){
        shearStepwiseIncreaseTest(sd_sys, dem);
    } else if ( dem.shear_process == "you_can_define_test"){
        // youCanDefineTest(sd_sys, dem);
    } else {
        cerr << dem.shear_process << " is not programed." << endl;
    }
    
	return;
}

void shearStepwiseIncreaseTest(SDsystem &sd_sys,
                               DEMsystem &dem){    
    double r_exponent = 1.0 /(dem.stepnumber_shearrate - 1);
    dem.incrementratio_shearrate = pow(dem.shearrate_max/dem.shearrate_min , r_exponent);
    dem.openOutputFileDEM();
	dem.outputYaplotDEM();
    dem.outputDeformationConf();
    dem.setInitialMotion();
    dem.shiftClusterToCenter();
    if (sd_sys.method_hydroint != 0){
        dem.getSDMovMatrix();
    }
    dem.initialprocess = false;
    dem.outputConfiguration();
    
    dem.shearrate = dem.shearrate_min;
    for(int m = 0; m < dem.stepnumber_shearrate; m++){
        /*
         * Change shear rate.
         */ 
        if (m > 0){
            dem.shearrate *= dem.incrementratio_shearrate;
        }
        cerr << "shear rate = "  << dem.shearrate << endl;
        dem.setStepSize();
        double shearstrain_end = dem.shearstrain + dem.shearstrain_step;
        int step_counter =0;
        do{
            double nexttime_output = dem.shearstrain + dem.interval_strain_output;
            dem.makeNeighbor();
            
            while (dem.shearstrain < nexttime_output){
                dem.calcInterParticleForce();
                if (sd_sys.method_hydroint == 0){
                    dem.moveOverdampedMotionFDA();
                } else {
                    dem.moveOverdampedMotionSDMov();
                }
                dem.simu_time += (dem.h_stepsize / dem.shearrate);
                dem.shearstrain += dem.h_stepsize;
                dem.checkFailure();
                if (!dem.regeneration_bond.empty()){
                    dem.regeneration_onebyone();
                    dem.regeneration_bond.clear();
                }
                dem.generateBond();
                if (sd_sys.method_hydroint != 0){
                    dem.estimateClusterRotation();
                    if ( dem.step_deformation > dem.critical_deformation_SD ){
                        cerr << "renew matrix" << ' ' << dem.step_deformation << endl;
                        dem.shiftClusterToCenter();                    
                        dem.getSDMovMatrix();
                        dem.resetDeformation();
                        if (dem.version[0] == 'o'){
                            dem.outputYaplotDEM();
                            dem.outputDeformationConf();
                        }
                        if (dem.step_deformation > 100){
                            cerr << "the simulation failed." << endl;
                            exit(1);
                        }                        
                    }
                }
                if ( step_counter++ % 200 == 0){
                    dem.shiftClusterToCenter(); 
                    dem.q_rot.normalize();
                    for (int i=0; i < sd_sys.np; i++){
                        dem.particle[i]->setNorm_u(); /*********** IMPORTANT ************/
                        dem.particle[i]->orientation.normalize(); /*********** IMPORTANT ************/
                    }
                }
            }
            dem.shiftClusterToCenter();
            dem.calcGyrationRadius();
            dem.checkState();
            dem.calcTotalFTS();
            dem.calcTotalDeformation();
            /*
             * Data output
             */
            dem.outputDeformationConf();
            dem.outputYaplotDEM();
            dem.outputLogDEM();
            dem.outputData();
            dem.outputConfiguration();
        } while (dem.shearstrain < shearstrain_end);    
    }
}

