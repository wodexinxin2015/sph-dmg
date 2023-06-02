/************************************************************************************
SPH Discrete Model Generator
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/
#include <stdio.h>
#include "SPHD_Class_Data.h"

#pragma once
class Particle_Gen2D {
public:
	int parti_gen2d(void);
private:
	int read_para(Model_Control *mdc, Model_Para *mdp);
	void outtofile(Particle *parti, int num, FILE *fop);
	void boundary_box(Particle *bdparti, Model_Para *mdp, double dr);
	void particle_shape(Particle *parti, Model_Para *mdp, double dr,
		int type, int matype, int permtype, int etype);
	void structure_2d(Particle *parti, Model_Para *mdp, double dr);
	void mixture1_2d(Particle *parti, Model_Para mdp, double base, double dr);
	void mixture2_2d(Particle *parti, Model_Para mdp, double base, double dr);
	void exca_back_2d(Particle *parti, Model_Para mdp, int type, int etype);
	void rain_shape(Particle* parti, Model_Para* mdp, double dr,
		int type, int matype, int permtype, int etype);
};

class Particle_Gen3D {
public:
	int parti_gen3d(void);
private:
	int read_para(Model_Control *mdc, Model_Para *mdp);
	void outtofile(Particle *parti, int num, FILE *fop);
	void boundary_box(Particle *bdparti, Model_Para *mdp, double dr);
	void particle_shape(Particle *parti, Model_Para *mdp, double dr,
		int type, int matype, int permtype, int etype);
	void structure_2d(Particle *parti, Model_Para *mdp, double dr);
	void structure_3d(Particle *parti, Model_Para mdp, double dr, int str3_flag);
	void mixture1_2d(Particle *parti, Model_Para mdp, double base, double dr);
	void mixture2_2d(Particle *parti, Model_Para mdp, double base, double dr);
	void exca_back_3d(Particle *parti, Model_Para mdp, int type, int etype);
	void extruding_3d(Particle *parti, Model_Para *mdp, int type, double ext_dist, double dr);
	void mixture_3d(Particle *parti, Model_Para mdp, int num_2d);
	void rain_shape(Particle* parti, Model_Para* mdp, double dr,
		int type, int matype, int permtype, int etype);
};