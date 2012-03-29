/*
 *  my_utilities.h
 *  DEMsd
 *
 *  Created by Ryohei Seto on 10/05/03.
 *  Copyright (c) 2010 Ryohei Seto. All rights reserved.
 *
 */

#ifndef my_utilities_h
#define my_utilities_h 1
#include <fstream>
#include <iomanip>
#include <sys/time.h>
using namespace std;

/* math tools */
inline double pow_10(int index){
	double val = 1.0;
	if ( index > 0 ){
		while ((index --) > 0){
			val *= 10;
		}
	} else {
		while ((index ++) < 0){
			val /= 10;
		}
	}
	return val ;
}

#endif





