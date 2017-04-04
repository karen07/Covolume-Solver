#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double dl, vl, pl, cl, dr, vr, pr, cr;
double p_s, v_s;
double gama;
double b;
double eps;

double solve_star_reg(double start, void(*fun)(double, double *, double *))
{
	double f, fd, x;

	x = start;
	for (int i = 1; i <= 20; i++) {
		fun(x, &f, &fd);
		x = x - f / fd;
		if (fabs(f) <= eps)
			break;
		if (x < 0.0)
			x = eps;
	}
	return x;
}

double adiabat_val(double p) {
	double v;

	if (p > pr) {
		v = vr + (p - pr) * sqrt(2 * (1 / dr - b) / (p * (gama + 1) + pr * (gama - 1)));
	}
	else {
		v = vr + 2.0 / (gama - 1.0) * cr * (1 - b * dr) * (pow(p / pr, (gama - 1.0) / (2.0 * gama)) - 1.0);
	}

	return v;
}

void adiabats_sub(double p, double *v, double *vd) {

	if (p > pr) {
		*v = vr + (p - pr) * sqrt(2 * (1 / dr - b) / (p * (gama + 1) + pr * (gama - 1)));
		*vd = sqrt(2 * (1 / dr - b) / (p * (gama + 1) + pr * (gama - 1))) - (gama + 1) * (p - pr) * sqrt(1 / dr - b) / (sqrt(2) * pow(p * (gama + 1) + pr * (gama - 1), 3.0 / 2.0));
	}
	else {
		*v = vr + 2.0 / (gama - 1.0) * cr * (1 - b * dr) * (pow(p / pr, (gama - 1.0) / (2.0 * gama)) - 1.0);
		*vd = cr / (gama * pr) * (1 - b * dr) * (pow(p / pr, (-gama - 1.0) / (2.0 * gama)));
	}

	if (p > pl) {
		*v -= vl - (p - pl) * sqrt(2 * (1 / dl - b) / (p * (gama + 1) + pl * (gama - 1)));
		*vd -= -sqrt(2 * (1 / dl - b) / (p * (gama + 1) + pl * (gama - 1))) + (gama + 1) * (p - pl) * sqrt(1 / dl - b) / (sqrt(2) * pow(p * (gama + 1) + pl * (gama - 1), 3.0 / 2.0));
	}
	else {
		*v -= vl - 2.0 / (gama - 1.0) * cl * (1 - b * dl) * (pow(p / pl, (gama - 1.0) / (2.0 * gama)) - 1.0);
		*vd -= -cl / (gama * pl) * (1 - b * dl) * (pow(p / pl, (-gama - 1.0) / (2.0 * gama)));
	}
}

double d_equation_l(double ro, double s) {
	double i = vl + 2 / (gama - 1) * cl * (1 - b * dl);

	double bet = pow((gama - 1) * (i - s), 2) / (gama * pl * pow((1 / dl - b), gama));

	return pow(ro, gama - 1) * pow(gama + 1 - 2 * b * ro, 2)
		- pow(1 - b * ro, gama + 1) * bet;
}

double d_equation_r(double ro, double s) {
	double i = vr - 2 / (gama - 1) * cr * (1 - b * dr);

	double bet = pow((gama - 1) * (i - s), 2) / (gama * pr * pow((1 / dr - b), gama));

	return pow(ro, gama - 1) * pow(gama + 1 - 2 * b * ro, 2)
		- pow(1 - b * ro, gama + 1) * bet;
}

double solve_d_equation(double ds, double d, double s, double(*ro_val)(double, double)) {
	double left_bord = d;
	double right_bord = ds;
	double middle;
	while (fabs(ro_val(left_bord, s)) > eps) {
		middle = left_bord - (left_bord - right_bord) / 2;
		if (signbit(ro_val(left_bord, s)) == signbit(ro_val(middle, s))) {
			left_bord = middle;
		}
		else {
			right_bord = middle;
		}
	}
	return left_bord;
}

void solution_profile(const double pm, const double um, const double s, double *d_r, double *v_r, double *p_r)
{
	double cs, flow_l, flow_r, shl, shr, sl, sr, stl, str;
	double d, v, p;

	if (s <= um) {
		// sampling point lies to the left of the contact discontinuity
		if (pm <= pl) {
			// left rarefaction
			shl = vl - cl;
			if (s <= shl) {
				d = dl;
				v = vl;
				p = pl;
			}
			else {
				d = 1 / (pow(pl / pm, 1 / gama) * (1 / dl - b) + b);
				cs = sqrt(pm * gama / (d * (1 - b * d)));
				stl = um - cs;
				if (s > stl) {
					v = um;
					p = pm;
				}
				else {
					d = solve_d_equation(d, dl, s, d_equation_l);
					p = pl * pow((1 / dl - b) / (1 / d - b), gama);
					v = vl - 2.0 / (gama - 1.0) * cl * (1 - b * dl) * (pow(p / pl, (gama - 1.0) / (2.0 * gama)) - 1.0);
				}
			}
		}
		else {
			// left shock
			flow_l = sqrt((pm * (gama + 1) + pl * (gama - 1)) / (2 / dl - 2 * b));
			d = (pm * (gama + 1) + pl * (gama - 1)) / (pm * (gama - 1) / dl + pl * (gama + 1) / dl + 2 * b * (pm - pl));
			sl = vl - flow_l / dl;
			if (s <= sl) {
				// sampled point is left data state
				d = dl;
				v = vl;
				p = pl;
			}
			else {
				// sampled point is star left state
				v = um;
				p = pm;
			}
		}
	}
	else {
		// sampling point lies to the right of the contact discontinuity
		if (pm > pr) {
			// right shock
			flow_r = sqrt((pm * (gama + 1) + pr * (gama - 1)) / (2 / dr - 2 * b));
			d = (pm * (gama + 1) + pr * (gama - 1)) / (pm * (gama - 1) / dr + pr * (gama + 1) / dr + 2 * b * (pm - pr));
			sr = vr + flow_r / dr;
			if (s >= sr) {
				// sampled point is right data state
				d = dr;
				v = vr;
				p = pr;
			}
			else {
				// sampled point is star right state
				v = um;
				p = pm;
			}
		}
		else {
			// right rarefaction
			shr = vr + cr;
			if (s >= shr) {
				// sampled point is right data state
				d = dr;
				v = vr;
				p = pr;
			}
			else {
				d = 1 / (pow(pr / pm, 1 / gama) * (1 / dr - b) + b);
				cs = sqrt(pm * gama / (d * (1 - b * d)));
				str = um + cs;
				if (s <= str) {
					// sampled point is star right state
					v = um;
					p = pm;
				}
				else {
					d = solve_d_equation(d, dr, s, d_equation_r);
					p = pr * pow((1 / dr - b) / (1 / d - b), gama);
					v = vr + 2.0 / (gama - 1.0) * cr * (1 - b * dr) * (pow(p / pr, (gama - 1.0) / (2.0 * gama)) - 1.0);
				}
			}
		}
	}

	*d_r = d;
	*v_r = v;
	*p_r = p;
}

void riman_solver_covolume(double dl_t, double vl_t, double pl_t, double dr_t, double vr_t, double pr_t, double b_t, double gama_t, double len, double dia_pos, double time, double cells, double eps_t, double **d_arr, double **v_arr, double **p_arr, double **x) {
	dl = dl_t; vl = vl_t; pl = pl_t; dr = dr_t; vr = vr_t; pr = pr_t;
	b = b_t;
	eps = eps_t;
	gama = gama_t;

	cl = sqrt(gama * pl / (dl * (1 - b * dl)));
	cr = sqrt(gama * pr / (dr * (1 - b * dr)));

	if (2.0 / (gama - 1.0) * (cl + cr) <= (vr - vl)) {
		*d_arr = 0;
		*v_arr = 0;
		*p_arr = 0;
		*x = 0;
		return;
	}

	p_s = solve_star_reg(pr, adiabats_sub);
	v_s = adiabat_val(p_s);

	double dx = len / double(cells);

	*d_arr = (double *)malloc(cells * sizeof(double));
	*v_arr = (double *)malloc(cells * sizeof(double));
	*p_arr = (double *)malloc(cells * sizeof(double));
	*x = (double *)malloc(cells * sizeof(double));

	double xpos, s, d_res, v_res, p_res;

	for (int i = 0; i < cells; i++) {
		xpos = (i - 0.5)*dx;
		s = (xpos - dia_pos) / time;
		solution_profile(p_s, v_s, s, &d_res, &v_res, &p_res);

		(*d_arr)[i] = d_res;
		(*v_arr)[i] = v_res;
		(*p_arr)[i] = p_res;
		(*x)[i] = xpos;
	}
}