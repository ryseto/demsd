/*
 *  contable.h
 *  DEMsd
 *
 *  Created by SETO Ryohei on 07/01/18.
 *  Copyright 2007 Ryohei Seto. All rights reserved.
 *
 */

#ifndef contable_h
#define contable_h 1

/*
 * Contact table:
 * When particle i and j are contacting each other,
 * contact_table[i][j] = 1.
 */

class ConTable{
	bool allocate;
	bool **contact_table;
	int n;
public:
	ConTable():allocate(false) {}
	~ConTable();
	inline bool connect(int i, int j){
		return contact_table[i][j];
	}

	void set(int particleNumber);
	void reset();
	void on_connect(int i, int j);
	void off_connect(int i, int j);	
};
#endif
