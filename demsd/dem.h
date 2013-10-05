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
#include "sdsystem.h"
#include <vector>
using namespace std;

void demSimulation(int argc, char** argv);
void shearStepwiseIncreaseTest(SDsystem &sd_sys, DEMsystem &dem);


void testContactModel(int argc, char** argv);