/*
 *  dem.h
 *  DEMsd
 *
 *  Created by Ryohei Seto on 10/07/07.
 *  Copyright 2010 Ryohei Seto. All rights reserved.
 *
 */
#include "demsystem.h"
#include "common.h"
#include "bond.h"
#include "particle.h"
#include "SDsystem.h"
#include <vector>
using namespace std;

void shearStepwiseIncreaseTest(SDsystem &sd_sys, DEMsystem &dem);

void demSimulation(int argc, char** argv){
    SDsystem sd_sys;    
    DEMsystem dem(sd_sys);
    dem.initialprocess = true;
    /* parameter file */
	dem.setParameterFileDEM(argv[2]);   
    /* Cluster configuration : num. of skip lines */
	dem.importCluster(argv[3], atoi(argv[4]));
    dem.setVersion(argv[5]);
    dem.readParameterFileDEM();
    dem.readParameterBond();
    dem.readParameterShearProcess();
    sd_sys.setFlowType('s');
	/* DEM */
	dem.setDEMParameters();
	dem.initDEM();
    /* Stokesian Dynamics */
    sd_sys.initLibStokes(dem.np);

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
    dem.openOutputFileDEM();
	dem.outputYaplotDEM();
    dem.outputDeformationConf();

	ofstream out_povray ;
    dem.setInitialMotion();
    dem.shiftClusterToCenter();
    dem.getMovMatrix();
    dem.initialprocess = false;
    dem.outputConfiguration();
    
    dem.shearrate = dem.shearrate_min;
    for(int m = 0; m < dem.stepnumber_shearrate; m++){
        if (m > 0){
            dem.shearrate *= dem.incrementratio_shearrate;
        }
        cerr << "shear rate = "  << dem.shearrate << endl;
        dem.setStepSize();
        double shearstrain_end = dem.shearstrain + dem.shearstrain_step;
        do{
            double nexttime_output = dem.shearstrain + dem.interval_strain_output;
            dem.makeNeighbor();
            int i =0;
            while (dem.shearstrain < nexttime_output){

                dem.calcInterParticleForce();
                if (sd_sys.method_hydroint == 0){
                    dem.freeDrainingApproximation();
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
//                    dem.outputYaplotDEM();
                    dem.estimateClusterRotation();
                    if ( dem.step_deformation > dem.critical_deformation_SD ){
                        cerr << "renew matrix" << ' ' << dem.step_deformation << endl;
                        dem.shiftClusterToCenter();                    
                        dem.getMovMatrix();
                        dem.resetDeformation();
                        //sd_sys.set_dr_from_sdpos();
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

                
                i ++; 
                if ( i % 200 == 0){
                    dem.shiftClusterToCenter(); 
                    dem.q_rot.normalize();
                    for (int i=0; i < sd_sys.np; i++){
                        dem.particle[i]->setNorm_u(); /*********** IMPORTANT ************/
                        dem.particle[i]->orientation.normalize(); /*********** IMPORTANT ************/
                    }
                }
                //   dem.outputYaplotDEM();
                
            }
            dem.shiftClusterToCenter();
            dem.calcGyrationRadius();
            dem.checkState();
            dem.calcTotalFTS();
            dem.calcTotalDeformation();
            dem.outputDeformationConf();
            dem.outputYaplotDEM();
            dem.outputLogDEM();
            dem.outputData();
            dem.outputConfiguration();
        }while (dem.shearstrain < shearstrain_end);    
    }

}