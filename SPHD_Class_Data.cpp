/************************************************************************************
SPH Discrete Model Generator
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/
#include <stdlib.h>
#include "SPHD_Class_Data.h"

Particle::Particle() {
	for (int i = 0; i < 3; i++) {
		x[i] = 0.0;
		vx[i] = 0.0;
	}

	type = -1;
	matype = -1;
	permtype = -1;
	etype = -1;
}

Particle::~Particle() {

}

Model_Control::Model_Control() {
	dim = 2;
	
	b_flag = 0;
	w_flag = 0;
	s_flag = 0;
	a_flag = 0;
	r_flag = 0;

	mix_flag = 0;
	str_flag = 0;
	str3_flag = 0;
	eob_flag = 0;

	dr = 0.0;
	ext_dst = 0.0;
}

Model_Control::~Model_Control() {

}

Model_Para::Model_Para() {
	n_boundary = 0;
	n_water_top = 0;
	n_soil_top = 0;
	n_air_top = 0;
	n_water_bot = 0;
	n_soil_bot = 0;
	n_air_bot = 0;

	num_boundary = 0;
	num_water = 0;
	num_soil = 0;
	num_air = 0;
	num_rain = 0;

	b_lb[0] = 0.0;
	b_lb[1] = 0.0;
	b_lb[2] = 0.0;
	b_rt[0] = 0.0;
	b_rt[1] = 0.0;
	b_rt[2] = 0.0;

	str_lb[0] = 0.0;
	str_lb[1] = 0.0;
	str_rt[0] = 0.0;
	str_rt[1] = 0.0;

	for (int i = 0; i < 10; i++) {
		str3_lb[i][0] = 0.0;
		str3_lb[i][1] = 0.0;
		str3_lb[i][2] = 0.0;
		str3_rt[i][0] = 0.0;
		str3_rt[i][1] = 0.0;
		str3_rt[i][2] = 0.0;
	}

	eob_lb[0] = 0.0;
	eob_lb[1] = 0.0;
	eob_rt[0] = 0.0;
	eob_rt[1] = 0.0;

	eob3_lb[0] = 0.0;
	eob3_lb[1] = 0.0;
	eob3_lb[2] = 0.0;
	eob3_rt[0] = 0.0;
	eob3_rt[1] = 0.0;

	left_down[0] = 0.0;
	left_down[1] = 0.0;
	left_down[2] = 0.0;

	right_up[0] = 0.0;
	right_up[1] = 0.0;
	right_up[2] = 0.0;

	b_bottom = NULL;
	b_top = NULL;

	w_bottom = NULL;
	w_top = NULL;

	s_bottom = NULL;
	s_top = NULL;

	a_bottom = NULL;
	a_top = NULL;

}

Model_Para::~Model_Para() {
	if (b_bottom != NULL) delete[]b_bottom;
	if (b_top != NULL) delete[]b_top;
	if (w_bottom != NULL) delete[]w_bottom;
	if (w_top != NULL) delete[]w_top;
	if (s_bottom != NULL) delete[]s_bottom;
	if (s_top != NULL) delete[]s_top;
	if (a_bottom != NULL) delete[]a_bottom;
	if (a_top != NULL) delete[]a_top;
}

void Model_Para::allocate_bdp() {
	if (n_boundary >= 1) {
		b_bottom = new double[n_boundary][2];
		b_top = new double[n_boundary][2];
	}
	else {
		b_bottom = NULL;
		b_top = NULL;
	}
}

void Model_Para::allocate_wtp() {
	if (n_water_bot >= 1) {
		w_bottom = new double[n_water_bot][2];
	}
	if (n_water_top >= 1) {
		w_top = new double[n_water_top][2];
	}
}

void Model_Para::allocate_slp() {
	if (n_soil_bot >= 1) {
		s_bottom = new double[n_soil_bot][2];
	}
	if (n_soil_top >= 1) {
		s_top = new double[n_soil_top][2];
	}
}

void Model_Para::allocate_arp() {
	if (n_air_bot >= 1) {
		a_bottom = new double[n_air_bot][2];
	}
	if (n_air_top >= 1) {
		a_top = new double[n_air_top][2];
	}
}