/*
 *  contable.cpp
 *  DEMsd
 *
 *  Created by Ryohei Seto on 07/01/18.
 *  Copyright 2007 Ryohei Seto. All rights reserved.
 *
 */

#include "contable.h"
#include <iostream>
#include <iomanip>

ConTable::~ConTable(){
	if (allocate){
		for (int i=0; i<n+2; ++i){
			delete [] contact_table[i];
		}
		delete [] contact_table;
	}
}

void ConTable::set(int num_of_particle){
	allocate = true;
	n = num_of_particle;
	contact_table = new bool* [n+2];
	for (int i=0; i < n+2; ++i){
		contact_table[i] = new bool[n+2];
	}
	for (int i=0; i < n+2; ++i){
		for (int j=0; j < i; ++j){
			contact_table[i][j] = false;
			contact_table[j][i] = false;
		}
		contact_table[i][i] = true;
	}
}

void ConTable::on_connect(int i, int j){
	contact_table[i][j] = true;
	contact_table[j][i] = true;
}

void ConTable::off_connect(int i, int j){
	contact_table[i][j] = false;
	contact_table[j][i] = false;
}

void ConTable::reset(){
	for (int i=0; i<n+2; ++i){
		for (int j=0; j<n+2; ++j){
			contact_table[i][j] = false;
		}
		contact_table[i][i] = true;
	}
}
