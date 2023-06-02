/************************************************************************************
SPH Discrete Model Generator
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/

#pragma once

class Particle {
public:

	//variables
	double x[3];  /*position vector */
	double vx[3];  /*previous position vector */
	int type; /*type of particle */
	int matype; /*sub type of material */
	int permtype; /*permeability type*/
	int etype; /*existing type: 0--not exist; 1--exist*/
	
	Particle();
	~Particle();
};

class Model_Control {
public:

	int dim; /*dimension of problem*/
	int b_flag; /*type of boundary generation: 1-2D boundary box; 2--irregular shape; 3-3D boundary box*/
	int w_flag; /*type of water particles: 0--not exist; 1--exist*/
	int s_flag; /*type of soil particles: 0--not exist; 1--exist*/
	int a_flag; /*type of air particles: 0--not exist; 1--exist*/
	int r_flag; /*type of rain particles: 0--not exist; 1--exist*/

	int mix_flag; /*mixture type: 1--layered; 2--particle by particle;
				  3--3D particle by particle*/
	int str_flag; /*structure type: 0--not exist; 1--exist*/
	int str3_flag; /*3D structure type: 0--not exist; 1--exist*/
	int eob_flag; /*excavation type: 
				  1--watre and soil; 2--water only; 3--soil only*/
	double dr; /*initial particle spacing*/
	double ext_dst; /*3D extruding distance*/

	Model_Control();
	~Model_Control();
};

class Model_Para {
public:
	int n_boundary; /*number of boundary turning points*/
	int n_water_top; /*number of water turning points: top surface*/
	int n_soil_top; /*number of soil turning points: top surface*/
	int n_air_top; /*number of air turning points: top surface*/
	int n_water_bot; /*number of water turning points: bottom surface*/
	int n_soil_bot; /*number of soil turning points: bottom surface*/
	int n_air_bot; /*number of air turning points: bottom surface*/

	int num_boundary;  /*number of boundary pariticles*/
	int num_water;  /*number of water pariticles*/
	int num_soil;  /*number of soil pariticles*/
	int num_air;  /*number of air pariticles*/
	int num_rain;  /*number of rain pariticles*/

	double b_lb[3]; /*left bottom of boundary*/
	double b_rt[3]; /*right top of boundary*/
	double str_lb[2]; /*left bottom of structure*/
	double str_rt[2]; /*right top of structure*/
	double str3_lb[10][3]; /*left bottom of structure for 3D*/
	double str3_rt[10][3]; /*right top of structure for 3D*/
	double eob_lb[2]; /*left bottom of excavation and backfill*/
	double eob_rt[2]; /*right top of sexcavation and backfill*/
	double eob3_lb[3]; /*left bottom of excavation and backfill for 3D*/
	double eob3_rt[3]; /*right top of excavation and backfill for 3D*/

	double (*b_bottom)[2]; /*bottom of the boundary turing points*/
	double (*b_top)[2]; /*top of the boundary turing points*/
	double (*w_bottom)[2]; /*bottom of the water turing points*/
	double (*w_top)[2]; /*top of the water turing points*/
	double (*s_bottom)[2]; /*bottom of the soil turing points*/
	double (*s_top)[2]; /*top of the soil turing points*/
	double (*a_bottom)[2]; /*bottom of the air turing points*/
	double (*a_top)[2]; /*top of the air turing points*/
	double left_down[3]; /*coordinate of left and down corner*/
	double right_up[3]; /*coordinate of right and up corner*/

	Model_Para();
	~Model_Para();

	void allocate_bdp(); //allocate memeory for boundary turning points
	void allocate_wtp(); //allocate memeory for water turning points
	void allocate_slp(); //allocate memeory for soil turning points
	void allocate_arp(); //allocate memeory for air turning points
};