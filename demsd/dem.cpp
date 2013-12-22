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
    /* lubrication correction is not used.
     */
    int lub_correction = 0;
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
	dem.initDEM();
    if (dem.method_hydroint > 0){
        sd_sys.setFlowType('s'); // Set shear flow in SD.
        sd_sys.initFlowModel(dem.np, lub_correction, false);
    }
    /* Main simulation */
	if ( dem.shear_process == "stepwise"){
        shearStepwiseIncreaseTest(sd_sys, dem);
    } else if ( dem.shear_process == "you_can_define_test"){
        // youCanDefineTest(sd_sys, dem);
    } else {
        cerr << dem.shear_process << " is not programmed." << endl;
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
    if (dem.method_hydroint == 1){
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
        dem.setStepSize();
        double shearstrain_end = dem.shearstrain + dem.shearstrain_step;
        int step_counter =0;
        do{
            double nexttime_output = dem.shearstrain + dem.interval_strain_output;
            dem.makeNeighbor();
            cerr << "ss=" << dem.shearstrain << endl;
            while (dem.shearstrain <= nexttime_output){
                dem.calcInterParticleForce();
                if (dem.method_hydroint == 0){
                    dem.moveOverdampedMotionFDAShear();
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
                if (dem.method_hydroint == 1){
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
                    for (int i=0; i < dem.np; i++){
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

void testContactModel(int argc, char** argv){
    //
    // UNFINISHED...
    //
    DEMsystem dem;
    dem.setBondFile(argv[2],argv[3]);
    dem.readParameterBond();
    dem.setTwoParticleSystem();
	dem.initDEM();
    dem.openOutputFileDEM();
    dem.outputYaplotDEM();
    dem.particle[1]->p.z += 0.3;
    dem.h_stepsize = 0.0001;
    dem.outputYaplotDEM();
    while(1){
        dem.calcInterParticleForce();
        dem.moveOverdampedMotionFDA();
        dem.outputYaplotDEM();
    }
	
}
