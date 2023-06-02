/************************************************************************************
SPH Discrete Model Generator
--Copyright (c) 2018-2023, Weijie ZHANG, GeoHohai, Hohai University.
************************************************************************************/
#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include "SPHD_Class_Functions.h"
#include "SPHD_option.h"

extern int err_t;

int Particle_Gen2D::parti_gen2d(void) {
	//variables
	int err_t;
	FILE *fop;

	//model variables
	Model_Control inst_mdc;
	Model_Para inst_mdp;

	//particle variables
	Particle* bd_particle = new Particle[size / 4];
	Particle* wt_particle = new Particle[size / 4];
	Particle* sl_particle = new Particle[size / 4];
	Particle* ar_particle = new Particle[size / 4];
	Particle *rn_particle = new Particle[size / 4];

	//read file
	inst_mdc.dim = 2;
	err_t = read_para(&inst_mdc, &inst_mdp);
	if (err_t > 0)  return err_t;

	//generate particles
	if (inst_mdc.b_flag == 1) {
		boundary_box(bd_particle, &inst_mdp, inst_mdc.dr);
	}
	else if (inst_mdc.b_flag == 2) {
		particle_shape(bd_particle, &inst_mdp, inst_mdc.dr, 0, 0, 0, 0);
	}
	//water
	if (inst_mdc.w_flag == 1) {
		particle_shape(wt_particle, &inst_mdp, inst_mdc.dr, 1, 0, 0, 1);

		//mixture
		if (inst_mdc.mix_flag == 1) {
			mixture1_2d(wt_particle, inst_mdp, inst_mdp.w_bottom[0][1], inst_mdc.dr);
		}
		else if (inst_mdc.mix_flag == 2) {
			mixture2_2d(wt_particle, inst_mdp, inst_mdp.w_bottom[0][0], inst_mdc.dr);
		}
	}
	//soil
	if (inst_mdc.s_flag == 1) {
		particle_shape(sl_particle, &inst_mdp, inst_mdc.dr, 2, 0, 1, 2);

		//structure
		if (inst_mdc.str_flag == 1) {
			structure_2d(sl_particle, &inst_mdp, inst_mdc.dr);
		}
	}
	//air
	if (inst_mdc.a_flag == 1) {
		particle_shape(ar_particle, &inst_mdp, inst_mdc.dr, 3, 0, 0, 3);
	}
	//rain
	if (inst_mdc.r_flag == 1) {
		rain_shape(rn_particle, &inst_mdp, inst_mdc.dr, 7, 0, 0, 1);
	}
	//eob
	if (inst_mdc.eob_flag == 1) { //both water and soil
		exca_back_2d(wt_particle, inst_mdp, 1, 7);
		exca_back_2d(sl_particle, inst_mdp, 2, 7);
	}
	else if (inst_mdc.eob_flag == 2) { //for water only
		exca_back_2d(wt_particle, inst_mdp, 1, 7);
	}
	else if (inst_mdc.eob_flag == 3) { //for soil only
		exca_back_2d(sl_particle, inst_mdp, 2, 7);
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
	delete []bd_particle;
	delete []wt_particle;
	delete []sl_particle;
	delete []ar_particle;
	delete []rn_particle;

	return 0;
}

void Particle_Gen2D::exca_back_2d(Particle *parti, Model_Para mdp, int type, int etype) {
	int ntotal = 0;
	if (type == 1) ntotal = mdp.num_water;
	if (type == 2) ntotal = mdp.num_soil;

	for (int i = 0; i < ntotal; i++) {
		if ((parti[i].x[0] >= mdp.eob_lb[0]) && (parti[i].x[0] <= mdp.eob_rt[0]))
		{
			if ((parti[i].x[1] >= mdp.eob_lb[1]) && (parti[i].x[1] <= mdp.eob_rt[1])) {
				parti[i].etype = etype;
			}
		}
	}
}

void Particle_Gen2D::mixture2_2d(Particle *parti, Model_Para mdp, double base, double dr) {
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
				if (temp1%2 == 0) flag = tflag;
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

void Particle_Gen2D::mixture1_2d(Particle *parti, Model_Para mdp, double base, double dr) {
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

void Particle_Gen2D::structure_2d(Particle *parti, Model_Para *mdp, double dr) {
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

void Particle_Gen2D::boundary_box(Particle *bdparti, Model_Para *mdp, double dr) {
	int k, l;
	Particle *out_parti = new Particle[size];
	Particle *iner_parti = new Particle[size];
	double c[2], d[2], x;
	double e[2], f[2], y;
	double temp;

	//outer scope
	c[0] = mdp->b_lb[0] - 3.0*dr;
	c[1] = mdp->b_lb[1] - 3.0*dr;

	d[0] = mdp->b_rt[0] + 3.0*dr;
	d[1] = mdp->b_rt[1] + 3.0*dr;

	//inner scope
	e[0] = mdp->b_lb[0];
	e[1] = mdp->b_lb[1];

	f[0] = mdp->b_rt[0];
	f[1] = mdp->b_rt[1];

	//outer particles
	k = 0;
	x = c[0];
	do {
		y = c[1];
		do {
			if (k < size) {
				out_parti[k].x[0] = x;
				out_parti[k].x[1] = y;
				out_parti[k].vx[0] = 0.0;
				out_parti[k].vx[1] = 0.0;
				out_parti[k].type = 0;
				out_parti[k].matype = 0;
				out_parti[k].permtype = 0;
				out_parti[k].etype = 0;
			}
			k++;
			y = y + dr;
		} while (y <= d[1]);
		x = x + dr;
	} while (x <= d[0]);

	//inner particles
	l = 0;
	x = e[0];
	do {
		y = e[1];
		do {
			if (l < size) {
				iner_parti[l].x[0] = x;
				iner_parti[l].x[1] = y;
				iner_parti[l].vx[0] = 0.0;
				iner_parti[l].vx[1] = 0.0;
			}
			l++;
			y = y + dr;
		} while (y <= f[1]);
		x = x + dr;
	} while (x <= f[0]);

	//delete the redundant particles
	temp = 0.0001*dr;
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < l; j++) {
			if ((fabs(out_parti[i].x[0] - iner_parti[j].x[0]) < temp)
				&& (fabs(out_parti[i].x[1] - iner_parti[j].x[1]) < temp)) {
				out_parti[i].type = 11;
			}
		}
	}

	//save information to variables
	l = 0;
	for (int i = 0; i < k; i++) {
		if (out_parti[i].type == 0) {
			bdparti[l].x[0] = out_parti[i].x[0];
			bdparti[l].x[1] = out_parti[i].x[1];
			bdparti[l].vx[0] = out_parti[i].vx[0];
			bdparti[l].vx[1] = out_parti[i].vx[1];
			bdparti[l].type = 0;
			bdparti[l].matype = 0;
			bdparti[l].permtype = 0;
			bdparti[l].etype = 0;
			l++;
		}
	}
	mdp->num_boundary = l;
	delete []out_parti;
	delete []iner_parti;
}

void Particle_Gen2D::particle_shape(Particle *parti, Model_Para *mdp, double dr,
	int type, int matype, int permtype, int etype) {
	double min_y;
	double x, y;
	int k, l;
	Particle *out_parti = new Particle[size];
	Particle *iner_parti = new Particle[size];
	double temp = 0.01*dr;
	double slope, inter;
	int ntotal_bot, ntotal_top;
	double (*p_b)[2], (*p_t)[2];
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
	min_y -= 2.0*dr;

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
				} while (y <= (inter + slope * x - 0.01* dr));
				x = x + dr;
			} while (x <= p_b[i + 1][0] + 0.01 * dr);
		}
	}

	//delete the redundant particles
	temp = 0.01*dr;
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

int Particle_Gen2D::read_para(Model_Control *mdc, Model_Para *mdp) {
	char c[3];
	FILE *pFile = NULL;

	pFile = fopen("input.txt", "r");
	if (pFile == NULL) return 3;

	//inputing of controlling flags
	rewind(pFile);
	while (!feof(pFile)) {
		err_t=fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
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
		err_t=fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
		if (c[0] == 'p' && c[1] == 'a') { //flags
			err_t=fscanf(pFile, "%d%*[^\n]%*c", &mdp->n_boundary);
			err_t=fscanf(pFile, "%d%d%*[^\n]%*c", &mdp->n_water_bot, &mdp->n_water_top);
			err_t=fscanf(pFile, "%d%d%*[^\n]%*c", &mdp->n_soil_bot, &mdp->n_soil_top);
			err_t=fscanf(pFile, "%d%d%*[^\n]%*c", &mdp->n_air_bot, &mdp->n_air_top);
			err_t=fscanf(pFile, "%lf%*[^\n]%*c", &mdc->dr);
			err_t=fscanf(pFile, "%lf%*[^\n]%*c", &mdc->ext_dst);
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
			err_t=fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == 'b' && c[1] == 'o') { //flags
				err_t=fscanf(pFile, "%lf%lf%lf%*[^\n]%*c", &mdp->b_lb[0], &mdp->b_lb[1], &mdp->b_lb[2]);
				err_t=fscanf(pFile, "%lf%lf%lf%*[^\n]%*c", &mdp->b_rt[0], &mdp->b_rt[1], &mdp->b_rt[2]);
				break;
			}
		}
	}
	else if (mdc->b_flag == 2) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t=fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == 'b' && c[1] == 'o') { //flags
				for (int i = 0; i < mdp->n_boundary; i++) {
					err_t=fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->b_bottom[i][0], &mdp->b_bottom[i][1]);
				}
				for (int i = 0; i < mdp->n_boundary; i++) {
					err_t=fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->b_top[i][0], &mdp->b_top[i][1]);
				}
				break;
			}
		}
	}

	if (mdc->w_flag > 0) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t=fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == 'w' && c[1] == 'a') { //flags
				for (int i = 0; i < mdp->n_water_bot; i++) {
					err_t=fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->w_bottom[i][0], &mdp->w_bottom[i][1]);
				}
				for (int i = 0; i < mdp->n_water_top; i++) {
					err_t=fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->w_top[i][0], &mdp->w_top[i][1]);
				}
				break;
			}
		}
	}

	if (mdc->s_flag > 0) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t=fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == 's' && c[1] == 'o') { //flags
				for (int i = 0; i < mdp->n_soil_bot; i++) {
					err_t=fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->s_bottom[i][0], &mdp->s_bottom[i][1]);
				}
				for (int i = 0; i < mdp->n_soil_top; i++) {
					err_t=fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->s_top[i][0], &mdp->s_top[i][1]);
				}
				break;
			}
		}
	}

	if (mdc->a_flag > 0) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t=fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == 'a' && c[1] == 'i') { //flags
				for (int i = 0; i < mdp->n_air_bot; i++) {
					err_t=fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->a_bottom[i][0], &mdp->a_bottom[i][1]);
				}
				for (int i = 0; i < mdp->n_air_top; i++) {
					err_t=fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->a_top[i][0], &mdp->a_top[i][1]);
				}
				break;
			}
		}
	}

	//inputing of other parameters
	if (mdc->str_flag > 0) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t=fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == '2' && c[1] == 'd') { //flags
				err_t=fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->str_lb[0], &mdp->str_lb[1]);
				err_t=fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->str_rt[0], &mdp->str_rt[1]);
				break;
			}
		}
	}

	if (mdc->eob_flag > 0) {
		rewind(pFile);
		while (!feof(pFile)) {
			err_t=fscanf(pFile, "%c%c%*[^\n]%*c", &c[0], &c[1]);
			if (c[0] == 'e' && c[1] == 'o') { //flags
				err_t=fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->eob_lb[0], &mdp->eob_lb[1]);
				err_t=fscanf(pFile, "%lf%lf%*[^\n]%*c", &mdp->eob_rt[0], &mdp->eob_rt[1]);
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
//output to file
void Particle_Gen2D::outtofile(Particle *parti, int num, FILE *fop) {
	for (int i = 0; i < num; i++) {
		fprintf(fop, "%7.3f %7.3f %7.3f %7.3f %4d %4d %4d %4d \n",
			parti[i].x[0], parti[i].x[1], parti[i].vx[0], parti[i].vx[1],
			parti[i].type, parti[i].matype, parti[i].permtype, parti[i].etype);
	}
}

//generating rain particles
void Particle_Gen2D::rain_shape(Particle* parti, Model_Para* mdp, double dr,
	int type, int matype, int permtype, int etype) {
	double p1[3], p2[3];
	double x, y;
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
	y = p1[1];
	do {
		x = p1[0];
		do {
			if (r_num < size / 4) {
				parti[r_num].x[0] = x;
				parti[r_num].x[1] = y;
				parti[r_num].vx[0] = 0.0;
				parti[r_num].vx[1] = 0.0;
				parti[r_num].type = type;
				parti[r_num].matype = matype;
				parti[r_num].permtype = permtype;
				parti[r_num].etype = etype;
			}
			else break;
			x = x + dr;
			r_num += 1;
		} while (x <= p2[0]);
		y = y + dr / 4.0;
	} while (y <= p2[1]);

	//update the rain particle number
	mdp->num_rain = r_num;
}