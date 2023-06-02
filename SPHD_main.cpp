/************************************************************************************
SPH Discrete Model Generator
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/
/************************************************************************************
2019-06-13: V1.1.1 update the SPH Simulating Tool version to 1.0.2
2020-01-01: V1.2.0 number of turning points in bottom can be different from that of top
2020-01-06: V1.2.1 version update
2020-06-29: V1.2.2 fix the array index overflow
2020-07-21: V1.3.1 update the product version and file version
2020-07-21: V1.3.1 the generation of 3D boundary box and 2D/3D rain particles
2020-07-21: V1.3.2 revise a bug in the particle generation of model bottom
2021-01-20: V1.3.3 revise a bug in 3D model generation
2021-04-07: V1.4.1 revise the generation process of 3D model
************************************************************************************/
# define _CRT_SECURE_NO_WARNINGS
#include "SPHD_Class_Functions.h"

int err_t;

void statement(void) {
	printf("-------------------------------------------------------------------\n");
	printf("SPH Discrete Model Generator V1.5.1\n");
	printf("--Copyright (c) 2016-2021 Weijie ZHANG, GeoHohai, Hohai University.\n");
	printf("-------------------------------------------------------------------\n");
	printf("2020-07-21: V1.3.1 update the product version and file version.\n");
	printf("2020-07-21: V1.3.1 the generation of 3D boundary box and 2D/3D rain particles.\n");
	printf("2020-07-21: V1.3.2 revise a bug in the particle generation of model bottom.\n");
	printf("2020-11-30: V1.3.3 revise a bug in 3D model generation.\n");
	printf("2021-01-20: V1.3.3 revise a bug in 3D model generation\n");
	printf("2021-04-07: V1.4.1 revise the generation process of 3D model\n");
	printf("2021-04-14: V1.5.1 add the generation of several piles (less than 10)\n");
	printf("-------------------------------------------------------------------\n");
}

int main(void) {
	int err_t = 0;
	int dim;

	//statement
	statement();

	//instance for the class
	Particle_Gen2D inst_gen2d;
	Particle_Gen3D inst_gen3d;

	//selecting dimension
	printf("Please input the dimension of problem:\n");
	err_t = scanf("%d", &dim);

	//running module
	if (dim == 2) {
		err_t = inst_gen2d.parti_gen2d();
	}
	else if (dim == 3) {
		err_t = inst_gen3d.parti_gen3d();
	}
	else {
		printf("Error: Dimension is not 2 or 3.\n");
		return 1;
	}

	//error dealing
	if (err_t != 0) {
		printf("Error: during the particle generation.\n");
		return err_t;
	}

	return 0;
}