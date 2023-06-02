/************************************************************************************
SPH Discrete Model Generator
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/
#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include "SPHD_Class_Functions.h"
#include "SPHD_option.h"

extern int err_t;

int Particle_Gen3D::parti_gen3d(void) {
	//variables
	int err_t, num_2d;
	FILE *fop;

	//model variables
	Model_Control inst_mdc;
	Model_Para inst_mdp;

	//particle variables
	Particle *bd_particle = new Particle[size / 4];
	Particle *wt_particle = new Particle[size / 4];
	Particle *sl_particle = new Particle[size / 4];
	Particle *ar_particle = new Particle[size / 4];
	Particle* rn_particle = new Particle[size / 4];

	//read file
	inst_mdc.dim = 3;
	err_t = read_para(&inst_mdc, &inst_mdp);
	if (err_t > 0)  return err_t;
	
	//generate particles
	num_2d = 0;
	if (inst_mdc.b_flag == 1) {
		boundary_box(bd_particle, &inst_mdp, inst_mdc.dr); //3D boxing boundary
	}
	else if (inst_mdc.b_flag == 2) {
		particle_shape(bd_particle, &inst_mdp, inst_mdc.dr, 0, 0, 0, 0);
		//extruding boundary particles
		extruding_3d(bd_particle, &inst_mdp, 0, inst_mdc.ext_dst, inst_mdc.dr);
	}
	//water
	if (inst_mdc.w_flag == 1) {
		particle_shape(wt_particle, &inst_mdp, inst_mdc.dr, 1, 0, 0, 1);
		num_2d = inst_mdp.num_water;

		//mixture
		if (inst_mdc.mix_flag == 1) {
			mixture1_2d(wt_particle, inst_mdp, inst_mdp.w_bottom[0][1], inst_mdc.dr);
		}
		else if (inst_mdc.mix_flag == 2 || inst_mdc.mix_flag == 3) {
			mixture2_2d(wt_particle, inst_mdp, inst_mdp.w_bottom[0][0], inst_mdc.dr);
		}

		//extruding water particles
		extruding_3d(wt_particle, &inst_mdp, 1, inst_mdc.ext_dst, inst_mdc.dr);

		//3D particle by particle mixture
		if (inst_mdc.mix_flag == 3) {
			mixture_3d(wt_particle, inst_mdp, num_2d);
		}
	}
	//soil
	if (inst_mdc.s_flag == 1) {
		particle_shape(sl_particle, &inst_mdp, inst_mdc.dr, 2, 0, 0, 2);

		//structure
		if (inst_mdc.str_flag == 1) {
			structure_2d(sl_particle, &inst_mdp, inst_mdc.dr);
		}

		//extruding soil particles
		extruding_3d(sl_particle, &inst_mdp, 2, inst_mdc.ext_dst, inst_mdc.dr);
	}
	//air
	if (inst_mdc.a_flag == 1) {
		particle_shape(ar_particle, &inst_mdp, inst_mdc.dr, 3, 0, 0, 3);

		//extruding air particles
		extruding_3d(ar_particle, &inst_mdp, 3, inst_mdc.ext_dst, inst_mdc.dr);
	}
	//rain
	if (inst_mdc.r_flag == 1) {
		rain_shape(rn_particle, &inst_mdp, inst_mdc.dr, 7, 0, 0, 1);
	}
	//3D structure
	if (inst_mdc.str3_flag >= 1) {
		structure_3d(sl_particle, inst_mdp, inst_mdc.dr, inst_mdc.str3_flag);
	}

	//3D exacavation and backfill
	if (inst_mdc.eob_flag == 1) { //both water and soil
		exca_back_3d(wt_particle, inst_mdp, 1, 7);
		exca_back_3d(sl_particle, inst_mdp, 2, 7);
	}
	else if (inst_mdc.eob_flag == 2) { //for water only
		exca_back_3d(wt_particle, inst_mdp, 1, 7);
	}
	else if (inst_mdc.eob_flag == 3) { //for soil only
		exca_back_3d(sl_particle, inst_mdp, 2, 7);
	}

	//output to the file
	fop = fopen("model.out", "w");
	outtofile(bd_particle, inst_mdp.num_boundary, fop);
	outtofile(wt_particle, inst_mdp.num_water, fop);
	outtofile(sl_particle, inst_mdp.num_soil, fop);
	outtofile(ar_particle, inst_mdp.num_air, fop);
	outtofile(rn_particle, inst_mdp.num_rain, fop);
	fclose(fop);

	//dellocate memory
	delete[]bd_particle;
	delete[]wt_particle;
	delete[]sl_particle;
	delete[]ar_particle;
	delete[]rn_particle;

	return 0;
}

void Particle_Gen3D::mixture_3d(Particle *parti, Model_Para mdp, int num_2d) {
	for (int i = num_2d; i < mdp.num_water; i++) {
		int temp = i - num_2d;
		if (parti[temp].matype == 1) parti[i].matype = 0;
		else if (parti[temp].matype == 0) parti[i].matype = 1;
	}
}

void Particle_Gen3D::extruding_3d(Particle *parti, Model_Para *mdp, int type, double ext_dist, double dr) {
	int ntotal, k, i;
	double z;
	Particle *temp_parti = new Particle[size];

	ntotal = 0;
	if (type == 0) ntotal = mdp->num_boundary;
	if (type == 1) ntotal = mdp->num_water;
	if (type == 2) ntotal = mdp->num_soil;
	if (type == 3) ntotal = mdp->num_air;	

	k = 0;
	z = 0.0;
	do {
		for (i = 0; i < ntotal; i++) {
			temp_parti[k].x[0] = z;
			temp_parti[k].x[1] = parti[i].x[0];
			temp_parti[k].x[2] = parti[i].x[1];
			temp_parti[k].vx[0] = 0.0;
			temp_parti[k].vx[1] = parti[i].vx[0];
			temp_parti[k].vx[2] = parti[i].vx[1];
			temp_parti[k].type = parti[i].type;
			temp_parti[k].matype = parti[i].matype;
			temp_parti[k].permtype = parti[i].permtype;
			temp_parti[k].etype = parti[i].etype;
			k = k + 1;
		}
		z = z + dr;
	} while (z <= ext_dist);

	for (i = 0; i < k; i++) {
		parti[i].x[0] = temp_parti[i].x[0];
		parti[i].x[1] = temp_parti[i].x[1];
		parti[i].x[2] = temp_parti[i].x[2];
		parti[i].vx[0] = temp_parti[i].vx[0];
		parti[i].vx[1] = temp_parti[i].vx[1];
		parti[i].vx[2] = temp_parti[i].vx[2];
		parti[i].type = temp_parti[i].type;
		parti[i].matype = temp_parti[i].matype;
		parti[i].permtype = temp_parti[i].permtype;
		parti[i].etype = temp_parti[i].etype;
	}

	if (type == 0) mdp->num_boundary = k;
	if (type == 1) mdp->num_water = k;
	if (type == 2) mdp->num_soil = k;
	if (type == 3) mdp->num_air = k;

	delete[]temp_parti;
}

void Particle_Gen3D::exca_back_3d(Particle *parti, Model_Para mdp, int type, int etype) {
	int ntotal = 0;
	if (type == 1) ntotal = mdp.num_water;
	if (type == 2) ntotal = mdp.num_soil;

	for (int i = 0; i < ntotal; i++) {
		if ((parti[i].x[0] >= mdp.eob3_lb[0]) && (parti[i].x[0] <= mdp.eob3_rt[0]))
		{
			if ((parti[i].x[1] >= mdp.eob3_lb[1]) && (parti[i].x[1] <= mdp.eob3_rt[1])) {
				if ((parti[i].x[2] >= mdp.eob3_lb[2]) && (parti[i].x[2] <= mdp.eob3_rt[2])) {
					parti[i].etype = etype;
				}
			}
		}
	}
}

void Particle_Gen3D::mixture2_2d(Particle *parti, Model_Para mdp, double base, double dr) {
	int flag, tflag;
	double corx;

	for (int i = 0; i < mdp.num_water; i++) {
		double dst = parti[i].x[0] - base;
		int temp = (int)(dst / dr);

		if (i != 0) {
			if ((parti[i].x[0] != parti[i - 1].x[0]) && (temp % 2 == 0)) {
				tflag = 1;
				flag = tflag;
				corx = parti[i].x[1];
			}
			else if ((parti[i].x[0] != parti[i - 1].x[0]) && (temp % 2 == 1)) {
				tflag = -1;
				flag = tflag;
				corx = parti[i].x[1];
			}
			else if (parti[i].x[0] == parti[i - 1].x[0]) {
				double dst1 = parti[i].x[1] - corx;
				int temp1 = (int)(dst1 / dr);
				if (temp1 % 2 == 0) flag = tflag;
				else flag = -tflag;
			}
		}
		else {
			flag = 1;
			tflag = 1;
			corx = parti[i].x[1];
		}

		if (flag == 1) {
			parti[i].matype = 1;
		}
		else if (flag == -1) {
			parti[i].matype = 0;
		}
	}
}

void Particle_Gen3D::mixture1_2d(Particle *parti, Model_Para mdp, double base, double dr) {
	for (int i = 0; i < mdp.num_water; i++) {
		double dst = parti[i].x[1] - base;
		int temp = (int)(dst / dr + 0.1);

		if (fabs(temp % 2) < 0.00001) {
			parti[i].matype = 1;
		}
		else {
			parti[i].matype = 0;
		}
	}
}

void Particle_Gen3D::structure_2d(Particle* parti, Model_Para* mdp, double dr) {
	int id;
	double x_low, x_up;
	double y_low, y_up;
	double x, y;

	//initialization
	x_low = mdp->str_lb[0];
	y_low = mdp->str_lb[1];
	x_up = mdp->str_rt[0] + 0.01 * dr;
	y_up = mdp->str_rt[1] + 0.01 * dr;

	//x loop
	id = mdp->num_soil;
	x = x_low;
	do {
		y = y_low;
		do {
			parti[id].x[0] = x;
			parti[id].x[1] = y;
			parti[id].vx[0] = 0.0;
			parti[id].vx[1] = 0.0;
			parti[id].type = 4;
			parti[id].matype = 6;
			parti[id].permtype = 0;
			parti[id].etype = 4;
			id += 1;
			y += dr;
		} while (y <= y_up);
		x += dr;
	} while (x <= x_up);
	//update the number of soil particles
	mdp->num_soil = id;
}

void Particle_Gen3D::structure_3d(Particle *parti, Model_Para mdp, double dr, int str3_flag) {
	double temp = 0.01*dr;
	for (int kt = 0; kt < str3_flag; kt++) {
		for (int i = 0; i < mdp.num_soil; i++) {
			if ((parti[i].x[0] >= (mdp.str3_lb[kt][0] - temp)) && (parti[i].x[0] <= (mdp.str3_rt[kt][0] + temp))) {
				if ((parti[i].x[1] >= (mdp.str3_lb[kt][1] - temp)) && (parti[i].x[1] <= (mdp.str3_rt[kt][1] + temp))) {
					if ((parti[i].x[2] >= (mdp.str3_lb[kt][2] - temp)) && (parti[i].x[2] <= (mdp.str3_rt[kt][2] + temp))) {
						parti[i].type = 4;
						parti[i].matype = 6;
						parti[i].permtype = 0;
						parti[i].etype = 4;
					}
				}
			}
		}
	}
}

void Particle_Gen3D::boundary_box(Particle *bdparti, Model_Para *mdp, double dr) {
	//variables
	double x, y, z;
	int b_num;
	double x_low, x_up;
	double y_low, y_up;
	double z_low, z_up;
	//initialization
	b_num = 0;
	//bottom layers
	x_low = mdp->b_lb[0] - 3.00 * dr;
	x_up = mdp->b_rt[0] + 3.00 * dr;
	y_low = mdp->b_lb[1] - 3.00 * dr;
	y_up = mdp->b_rt[1] + 3.00 * dr;
	z_low = mdp->b_lb[2] - 3.00 * dr;
	z_up = mdp->b_lb[2] - 0.99 * dr;
	//z loop
	z = z_low;
	do {
		//y loop
		y = y_low;
		do {
			//x loop
			x = x_low;
			do {
				if (b_num < size / 4) {
					bdparti[b_num].x[0] = x;
					bdparti[b_num].x[1] = y;
					bdparti[b_num].x[2] = z;
					bdparti[b_num].vx[0] = 0.0;
					bdparti[b_num].vx[1] = 0.0;
					bdparti[b_num].vx[2] = 0.0;
					bdparti[b_num].type = 0;
					bdparti[b_num].matype = 0;
					bdparti[b_num].permtype = 0;
					bdparti[b_num].etype = 0;
				}
				else break;
				x = x + dr;
				b_num += 1;
			} while (x <= x_up);
			//x loop
			y = y + dr;
		} while (y <= y_up);
		//y loop
		z = z + dr;
	} while (z <= z_up);
	//z loop
	//front layers
	x_low = mdp->b_rt[0] + 1.00 * dr;
	x_up = mdp->b_rt[0] + 3.01 * dr;
	y_low = mdp->b_lb[1] - 3.00 * dr;
	y_up = mdp->b_rt[1] + 3.00 * dr;
	z_low = mdp->b_lb[2];
	z_up = mdp->b_rt[2];
	//z loop
	z = z_low;
	do {
		//y loop
		y = y_low;
		do {
			//x loop
			x = x_low;
			do {
				if (b_num < size / 4) {
					bdparti[b_num].x[0] = x;
					bdparti[b_num].x[1] = y;
					bdparti[b_num].x[2] = z;
					bdparti[b_num].vx[0] = 0.0;
					bdparti[b_num].vx[1] = 0.0;
					bdparti[b_num].vx[2] = 0.0;
					bdparti[b_num].type = 0;
					bdparti[b_num].matype = 0;
					bdparti[b_num].permtype = 0;
					bdparti[b_num].etype = 0;
				}
				else break;
				x = x + dr;
				b_num += 1;
			} while (x <= x_up);
			//x loop
			y = y + dr;
		} while (y <= y_up);
		//y loop
		z = z + dr;
	} while (z <= z_up);
	//z loop
	//back layers
	x_low = mdp->b_lb[0] - 3.00 * dr;
	x_up = mdp->b_lb[0] - 0.99 * dr;
	y_low = mdp->b_lb[1] - 3.00 * dr;
	y_up = mdp->b_rt[1] + 3.00 * dr;
	z_low = mdp->b_lb[2];
	z_up = mdp->b_rt[2];
	//z loop
	z = z_low;
	do {
		//y loop
		y = y_low;
		do {
			//x loop
			x = x_low;
			do {
				if (b_num < size / 4) {
					bdparti[b_num].x[0] = x;
					bdparti[b_num].x[1] = y;
					bdparti[b_num].x[2] = z;
					bdparti[b_num].vx[0] = 0.0;
					bdparti[b_num].vx[1] = 0.0;
					bdparti[b_num].vx[2] = 0.0;
					bdparti[b_num].type = 0;
					bdparti[b_num].matype = 0;
					bdparti[b_num].permtype = 0;
					bdparti[b_num].etype = 0;
				}
				else break;
				x = x + dr;
				b_num += 1;
			} while (x <= x_up);
			//x loop
			y = y + dr;
		} while (y <= y_up);
		//y loop
		z = z + dr;
	} while (z <= z_up);
	//z loop
	//left layers
	x_low = mdp->b_lb[0];
	x_up = mdp->b_rt[0] + 0.01 * dr;
	y_low = mdp->b_lb[1] - 3.00 * dr;
	y_up = mdp->b_lb[1] - 0.99 * dr;
	z_low = mdp->b_lb[2];
	z_up = mdp->b_rt[2];
	//z loop
	z = z_low;
	do {
		//y loop
		y = y_low;
		do {
			//x loop
			x = x_low;
			do {
				if (b_num < size / 4) {
					bdparti[b_num].x[0] = x;
					bdparti[b_num].x[1] = y;
					bdparti[b_num].x[2] = z;
					bdparti[b_num].vx[0] = 0.0;
					bdparti[b_num].vx[1] = 0.0;
					bdparti[b_num].vx[2] = 0.0;
					bdparti[b_num].type = 0;
					bdparti[b_num].matype = 0;
					bdparti[b_num].permtype = 0;
					bdparti[b_num].etype = 0;
				}
				else break;
				x = x + dr;
				b_num += 1;
			} while (x <= x_up);
			//x loop
			y = y + dr;
		} while (y <= y_up);
		//y loop
		z = z + dr;
	} while (z <= z_up);
	//z loop
	//right layers
	x_low = mdp->b_lb[0];
	x_up = mdp->b_rt[0] + 0.01 * dr;
	y_low = mdp->b_rt[1] + 1.00 * dr;
	y_up = mdp->b_rt[1] + 3.01 * dr;
	z_low = mdp->b_lb[2];
	z_up = mdp->b_rt[2];
	//z loop
	z = z_low;
	do {
		//y loop
		y = y_low;
		do {
			//x loop
			x = x_low;
			do {
				if (b_num < size / 4) {
					bdparti[b_num].x[0] = x;
					bdparti[b_num].x[1] = y;
					bdparti[b_num].x[2] = z;
					bdparti[b_num].vx[0] = 0.0;
					bdparti[b_num].vx[1] = 0.0;
					bdparti[b_num].vx[2] = 0.0;
					bdparti[b_num].type = 0;
					bdparti[b_num].matype = 0;
					bdparti[b_num].permtype = 0;
					bdparti[b_num].etype = 0;
				}
				else break;
				x = x + dr;
				b_num += 1;
			} while (x <= x_up);
			//x loop
			y = y + dr;
		} while (y <= y_up);
		//y loop
		z = z + dr;
	} while (z <= z_up);
	//z loop
	//update the number of boundary particles
	mdp->num_boundary = b_num;
}

void Particle_Gen3D::particle_shape(Particle *parti, Model_Para *mdp, double dr,
	int type, int matype, int permtype, int etype) {
	double min_y;
	double x, y;
	int k, l;
	Particle* out_parti = new Particle[size];
	Particle* iner_parti = new Particle[size];
	double temp = 0.01 * dr;
	double slope, inter;
	int ntotal_bot, ntotal_top;
	double(*p_b)[2], (*p_t)[2];
	int i, j;

	if (type == 0) {
		ntotal_top = mdp->n_boundary;
		ntotal_bot = ntotal_top;
		p_b = mdp->b_bottom;
		p_t = mdp->b_top;
	}
	else if (type == 1) {
		ntotal_bot = mdp->n_water_bot;
		ntotal_top = mdp->n_water_top;
		p_b = mdp->w_bottom;
		p_t = mdp->w_top;
	}
	else if (type == 2) {
		ntotal_bot = mdp->n_soil_bot;
		ntotal_top = mdp->n_soil_top;
		p_b = mdp->s_bottom;
		p_t = mdp->s_top;
	}
	else if (type == 3) {
		ntotal_bot = mdp->n_air_bot;
		ntotal_top = mdp->n_air_top;
		p_b = mdp->a_bottom;
		p_t = mdp->a_top;
	}
	else
	{
		ntotal_bot = 0;
		ntotal_top = 0;
		p_b = mdp->b_bottom;
		p_t = mdp->b_top;
	}

	//find the minimum of y coordinate
	min_y = 1000000.0;
	for (i = 0; i < ntotal_bot; i++) {
		min_y = fmin(min_y, p_b[i][1]);
	}
	min_y -= 2.0 * dr;

	//particles from baseline to top
	k = 0;
	for (i = 0; i < ntotal_top - 1; i++) {
		if (fabs(p_t[i + 1][0] - p_t[i][0]) > temp) {
			slope = (p_t[i + 1][1] - p_t[i][1])
				/ (p_t[i + 1][0] - p_t[i][0]);
			inter = p_t[i][1] - p_t[i][0] * slope;

			//generating
			x = p_t[i][0];
			do {
				y = min_y;

				do {
					if (k < size) {
						out_parti[k].x[0] = x;
						out_parti[k].x[1] = y;
						out_parti[k].vx[0] = 0.0;
						out_parti[k].vx[1] = 0.0;
						out_parti[k].type = type;
						out_parti[k].matype = matype;
						out_parti[k].permtype = permtype;
						out_parti[k].etype = etype;
					}
					else break;
					k = k + 1;
					y = y + dr;
				} while (y <= (inter + slope * x + 0.01 * dr));
				x = x + dr;
			} while (x <= p_t[i + 1][0] - 0.01 * dr);
		}
	}

	//particles from baseline to bottom
	l = 0;
	for (i = 0; i < ntotal_bot - 1; i++) {
		if (fabs(p_b[i + 1][0] - p_b[i][0]) > temp) {
			slope = (p_b[i + 1][1] - p_b[i][1])
				/ (p_b[i + 1][0] - p_b[i][0]);
			inter = p_b[i][1] - p_b[i][0] * slope;

			//generating
			x = p_b[i][0];
			do {
				y = min_y;
				do {
					if (l < size) {
						iner_parti[l].x[0] = x;
						iner_parti[l].x[1] = y;
						iner_parti[l].vx[0] = 0.0;
						iner_parti[l].vx[1] = 0.0;
					}
					else break;

					l = l + 1;
					y = y + dr;
				} while (y <= (inter + slope * x - 0.01 * dr));
				x = x + dr;
			} while (x <= p_b[i + 1][0] + 0.01 * dr);
		}
	}

	//delete the redundant particles
	temp = 0.01 * dr;
	for (i = 0; i < k; i++) {
		for (j = 0; j < l; j++) {
			if ((fabs(out_parti[i].x[0] - iner_parti[j].x[0]) < temp)
				&& (fabs(out_parti[i].x[1] - iner_parti[j].x[1]) < temp)) {
				out_parti[i].type = 100;
			}
		}
	}

	//save information to variables
	l = 0;
	for (i = 0; i < k; i++) {
		if (out_parti[i].type == type) {
			parti[l].x[0] = out_parti[i].x[0];
			parti[l].x[1] = out_parti[i].x[1];
			parti[l].vx[0] = out_parti[i].vx[0];
			parti[l].vx[1] = out_parti[i].vx[1];
			parti[l].type = out_parti[i].type;
			parti[l].matype = out_parti[i].matype;
			parti[l].permtype = out_parti[i].permtype;
			parti[l].etype = out_parti[i].etype;
			l++;
		}
	}

	if (type == 0) mdp->num_boundary = l;
	if (type == 1) mdp->num_water = l;
	if (type == 2) mdp->num_soil = l;
	if (type == 3) mdp->num_air = l;

	delete[]out_parti;
	delete[]iner_parti;
}

int Particle_Gen3D::read_para(Model_Control *mdc, Model_Para *mdp) {
	char c[3];
	FILE* pFile = NULL;

	pFile = fopen("input.txt", "r");
	if (pFile == NULL) return 3;

	//inputing of controlling flags
	rewind(pFile);
	while (!feof(pFile)) {
		err_t = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'c' && c[1] == 'o') { //flags
			err_t = fscanf(pFile, "%d%*[^\n]%*c", &mdc->b_flag);
			err_t = fscanf(pFile, "%d%*[^\n]%*c", &mdc->w_flag);
			err_t = fscanf(pFile, "%d%*[^\n]%*c", &mdc->s_flag);
			err_t = fscanf(pFile, "%d%*[^\n]%*c", &mdc->a_flag);
			err_t = fscanf(pFile, "%d%*[^\n]%*c", &mdc->r_flag);
			err_t = fscanf(pFile, "%d%*[^\n]%*c", &mdc->mix_flag);
			err_t = fscanf(pFile, "%d%*[^\n]%*c", &mdc->str_flag);
			err_t = fscanf(pFile, "%d%*[^\n]%*c", &mdc->str3_flag);
			err_t = fscanf(pFile, "%d%*[^\n]%*c", &mdc->eob_flag);
			break;
		}
	}

	//inputing of numbers of turning points
	rewind(pFile);
	while (!feof(pFile)) {
		err_t = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'p' && c[1] == 'a') { //flags
			err_t = fscanf(pFile, "%d%*[^\n]%*c", &mdp->n_boundary);
			err_t = fscanf(pFile, "%d%d%*[^\n]%*c", &mdp->n_water_bot, &mdp->n_water_top);
			err_t = fscanf(pFile, "%d%d%*[^\n]%*c", &mdp->n_soil_bot, &mdp->n_soil_top);
			err_t = fscanf(pFile, "%d%d%*[^\n]%*c", &mdp->n_air_bot, &mdp->n_air_top);
			err_t = fscanf(pFile, "%lf%*[^\n]%*c", &mdc->dr);
			err_t = fscanf(pFile, "%lf%*[^\n]%*c", &mdc->ext_dst);
			break;
		}
	}

	//allocate memory for turning points
	mdp->allocate_bdp();
	mdp->allocate_wtp();
	mdp->allocate_slp();
	mdp->allocate_arp();

	//inputing of turning points
	if (mdc->b_flag == 1) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == 'b' && c[1] == 'o') { //flags
				err_t = fscanf(pFile, "%lf%lf%lf%*[^\n]%*c", &mdp->b_lb[0], &mdp->b_lb[1], &mdp->b_lb[2]);
				err_t = fscanf(pFile, "%lf%lf%lf%*[^\n]%*c", &mdp->b_rt[0], &mdp->b_rt[1], &mdp->b_rt[2]);
				break;
			}
		}
	}
	else if (mdc->b_flag == 2) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == 'b' && c[1] == 'o') { //flags
				for (int i = 0; i < mdp->n_boundary; i++) {
					err_t = fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->b_bottom[i][0], &mdp->b_bottom[i][1]);
				}
				for (int i = 0; i < mdp->n_boundary; i++) {
					err_t = fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->b_top[i][0], &mdp->b_top[i][1]);
				}
				break;
			}
		}
	}

	if (mdc->w_flag > 0) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == 'w' && c[1] == 'a') { //flags
				for (int i = 0; i < mdp->n_water_bot; i++) {
					err_t = fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->w_bottom[i][0], &mdp->w_bottom[i][1]);
				}
				for (int i = 0; i < mdp->n_water_top; i++) {
					err_t = fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->w_top[i][0], &mdp->w_top[i][1]);
				}
				break;
			}
		}
	}

	if (mdc->s_flag > 0) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == 's' && c[1] == 'o') { //flags
				for (int i = 0; i < mdp->n_soil_bot; i++) {
					err_t = fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->s_bottom[i][0], &mdp->s_bottom[i][1]);
				}
				for (int i = 0; i < mdp->n_soil_top; i++) {
					err_t = fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->s_top[i][0], &mdp->s_top[i][1]);
				}
				break;
			}
		}
	}

	if (mdc->a_flag > 0) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == 'a' && c[1] == 'i') { //flags
				for (int i = 0; i < mdp->n_air_bot; i++) {
					err_t = fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->a_bottom[i][0], &mdp->a_bottom[i][1]);
				}
				for (int i = 0; i < mdp->n_air_top; i++) {
					err_t = fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->a_top[i][0], &mdp->a_top[i][1]);
				}
				break;
			}
		}
	}

	//inputing of other parameters
	if (mdc->str_flag > 0) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == '2' && c[1] == 'd') { //flags
				err_t = fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->str_lb[0], &mdp->str_lb[1]);
				err_t = fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->str_rt[0], &mdp->str_rt[1]);
				break;
			}
		}
	}

	if (mdc->str3_flag > 0) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == '3' && c[1] == 'd') { //flags
				for (int kt = 0; kt < mdc->str3_flag; kt++) {
					err_t = fscanf(pFile, "%lf%lf%lf%*[^\n]%*c", &mdp->str3_lb[kt][0], &mdp->str3_lb[kt][1], &mdp->str3_lb[kt][2]);
					err_t = fscanf(pFile, "%lf%lf%lf%*[^\n]%*c", &mdp->str3_rt[kt][0], &mdp->str3_rt[kt][1], &mdp->str3_rt[kt][2]);
				}
				break;
			}
		}
	}

	if (mdc->eob_flag > 0) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == 'e' && c[1] == 'o') { //flags
				err_t = fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->eob_lb[0], &mdp->eob_lb[1]);
				err_t = fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->eob_rt[0], &mdp->eob_rt[1]);
				break;
			}
		}
	}

	if (mdc->r_flag > 0) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t = fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == 'r' && c[1] == 'a') { //flags
				err_t = fscanf(pFile, "%lf%lf%lf%*[^\n]%*c", &mdp->left_down[0], &mdp->left_down[1], &mdp->left_down[2]);
				err_t = fscanf(pFile, "%lf%lf%lf%*[^\n]%*c", &mdp->right_up[0], &mdp->right_up[1], &mdp->right_up[2]);
				break;
			}
		}
	}

	fclose(pFile);
	return 0;
}

void Particle_Gen3D::outtofile(Particle *parti, int num, FILE *fop) {
	for (int i = 0; i < num; i++) {
		fprintf(fop, "%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %4d %4d %4d %4d \n",
			parti[i].x[0], parti[i].x[1], parti[i].x[2], parti[i].vx[0], parti[i].vx[1], parti[i].vx[2],
			parti[i].type, parti[i].matype, parti[i].permtype, parti[i].etype);
	}
}

//generating rain particles
void Particle_Gen3D::rain_shape(Particle* parti, Model_Para* mdp, double dr,
	int type, int matype, int permtype, int etype) {
	double p1[3], p2[3];
	double x, y, z;
	int r_num;

	//initialization of varibales
	p1[0] = mdp->left_down[0];
	p1[1] = mdp->left_down[1];
	p1[2] = mdp->left_down[2];

	p2[0] = mdp->right_up[0];
	p2[1] = mdp->right_up[1];
	p2[2] = mdp->right_up[2];

	//loop for rain particles
	r_num = 0;
	//z loop
	z = p1[2];
	do {
		//y loop
		y = p1[1];
		do {
			//x loop
			x = p1[0];
			do {
				if (r_num < size / 4) {
					parti[r_num].x[0] = x;
					parti[r_num].x[1] = y;
					parti[r_num].x[2] = z;
					parti[r_num].vx[0] = 0.0;
					parti[r_num].vx[1] = 0.0;
					parti[r_num].vx[2] = 0.0;
					parti[r_num].type = type;
					parti[r_num].matype = matype;
					parti[r_num].permtype = permtype;
					parti[r_num].etype = etype;
				}
				else break;
				x = x + dr;
				r_num += 1;
			} while (x <= p2[0]);
			//x loop
			y = y + dr;
		} while (y <= p2[1]);
		//y loop
		z = z + dr / 4.0;
	} while (z <= p2[2]);
	//z loop
	//update the rain particle number
	mdp->num_rain = r_num;
}