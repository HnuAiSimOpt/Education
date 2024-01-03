#include "mathfem.h"
#include <cmath>
#include "DoubleArray.h"
#define CUBIC_ZERO 1.0e-100
//Õª×ÔOOFEM mathfem.cÎÄ¼þ
void cubic(double a, double b, double c, double d, double *r1, double *r2, double *r3, int *num)
{
	double aa, p, q, D, u, v, phi;
	double help;

	double norm = 1e-6 * (fabs(a) + fabs(b) + fabs(c)) + CUBIC_ZERO;
	if (fabs(a) <= norm) {
		if (fabs(b) <= norm) {
			if (fabs(c) <= norm) {
				if (fabs(d) <= norm) {
					*r1 = 0.0;
					*num = 1;
					return;
				}
				else {
					*num = 0;
					return;
				}
			}
			else {
				*r1 = -d / c;
				*num = 1;
				return;
			}
		}
		else {
			if ((D = c * c - 4.0 * b * d) < 0.0) {
				*num = 0;
				return;
			}
			else {
				//*r1 = (-c + sqrt(D)) / 2.0 / b;
				//*r2 = (-c - sqrt(D)) / 2.0 / b;
				if (fabs(c) < norm) {
					help = -d / b;
					if (help > 0.) {
						*r1 = sqrt(help);
						*r2 = -sqrt(help);
						*num = 2;
						return;
					}
					else {
						*num = 0;
						return;
					}
				}
				else {
					help = -(c + sgn(c) * sqrt(D)) / 2.0;
					*r1 = help / b;
					*r2 = d / help;
					*num = 2;
				}
			}
		}
	}
	else {
		aa = a;
		a = b / (aa * 3.0);
		b = c / aa;
		c = d / aa;
		p = b - a * a * 3.0;
		q = 2.0 * a * a * a - a * b + c;
		D = q * q / 4.0 + p * p * p / 27.0;
		if (fabs(D) < CUBIC_ZERO) {
			if (fabs(p * q) < CUBIC_ZERO) {
				*r1 = 0.0 - a;
				*r2 = *r1;
				*r3 = *r1;
				*num = 3;
			}
			else {
				*r2 = cbrt(q / 2.0) - a;
				*r1 = -2.0 * * r2 - a;
				*num = 2;
			}
		}
		else {
			if (D > 0.0) {
				u = -q / 2.0 + sqrt(D);
				v = -q - u;

				*r1 = u + v - a;
				*r1 = cbrt(u) + cbrt(v) - a;
				*num = 1;
			}
			else {
				p = sqrt(fabs(p) / 3.0);
				help = (-q / (2.0 * p * p * p));
				if (fabs(help) > 1.0) {
					help = sgn(help);            // prevent rounding errors
				}

				phi = acos(help) / 3.0;
				double cp = cos(phi);
				double sp = sqrt(3.0) * sin(phi);
				*r1 = 2 * p * cp - a;
				*r2 = -p * (cp + sp) - a;
				*r3 = -p * (cp - sp) - a;

				// I'm getting some pretty bad accuracy, a single iteration like this would help alot
				//* r1 -= (d + c*(*r1) + b*(*r1)*(*r1) + a*(*r1)*(*r1)*(*r1))/(c + 2*b*(*r1) + 3*a*(*r1)*(*r1));
				//* r2 -= (d + c*(*r2) + b*(*r2)*(*r2) + a*(*r2)*(*r2)*(*r2))/(c + 2*b*(*r2) + 3*a*(*r2)*(*r2));
				//* r3 -= (d + c*(*r3) + b*(*r3)*(*r3) + a*(*r3)*(*r3)*(*r3))/(c + 2*b*(*r3) + 3*a*(*r3)*(*r3));
				* num = 3;
			}
		}
	}
}

void eigOfSymmtricMatrix(DoubleMatrix& vec, DoubleArray& val, DoubleMatrix src)
{
	if (!src.isSquare())
	{
		ERROR("cannot calculate the eigvals of !");
	}
	int n;
	double so, thresh, g, h, t, theta, c, s, z;
	n = src.GetColumn();
	val.resize(n);
	vec.resize(n, n);
	vec.beUnitMatrix();
	for (int i = 1; i <= n; i++)
	{
		val(i) = src(i, i);
	}
	for (int nite = 1; nite <= 50; nite++)
	{
		so = 0.;
		for (int p = 1; p <= n - 1; p++)
		{
			for (int q = p + 1; q <= n; q++)
			{
				so += abs(src(p, q));
			}
		}
		// convergence
		if (so <= 1e-12)
		{
			return;
		}
		if (nite < 4)
		{
			thresh = 0.2 * so / (n * n);
		}
		else
		{
			thresh = 0.;
		}
		// do sweep
		for (int p = 1; p <= n - 1; p++)
		{
			for (int q = p + 1; q <= n; q++)
			{
				g = 100.0 * abs(src(p, q));
				if ((nite > 4) && (abs(val(p)) + g == abs(val(p))) && (abs(val(q)) + g == abs(val(q))))
				{
					src(p, q) = 0.;
				}
				else if (abs(src(p,q)) > thresh)
				{
					h = val(q) - val(p);
					if (abs(h) + g == abs(h))
					{
						t = src(p, q) / h;
					}
					else
					{
						theta = 0.5 * h / src(p, q);
						if (theta < 0.)
						{
							t = -1 / (sqrt(1. + theta * theta) - theta);
						}
						else
						{
							t = 1 / (sqrt(1. + theta * theta) + theta);
						}
					}
					c = 1. / sqrt(1. + t * t);
					s = t * c;
					z = t * src(p, q);

					src(p, q) = 0.;
					val(p) = val(p) - z;
					val(q) = val(q) + z;
					for (int r = 1; r <= p - 1; r++)
					{
						t = src(r, p);
						src(r, p) = c * t - s * src(r, q);
						src(r, q) = s * t + c * src(r, q);
					}
					for (int r = p + 1; r <= q - 1; r++)
					{
						t = src(p, r);
						src(p, r) = c * t - s * src(r, q);
						src(r, q) = s * t + c * src(r, q);
					}
					for (int r = q + 1; r <= n; r++)
					{
						t = src(p, r);
						src(p, r) = c * t - s * src(q, r);
						src(q, r) = s * t + c * src(q, r);
					}
					for (int r = 1; r <= n; r++)
					{
						t = vec(r, p);
						vec(r, p) = c * t - s * vec(r, q);
						vec(r, q) = s * t + c * vec(r, q);
					}
				}
			}
		}
	}
	ERROR("not convergented!");
}
