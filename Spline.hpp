#pragma once
#include <cstdio>

static const int SPLINE_MAX_INPUT = 1024;

// Thomasアルゴリズム(https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm)で3重対角行列を解く
// a:下, b:対角, c:上
void triDiagonalMatrixSolve(double a[], double b[], double c[], double d[], double x[], int n)
{
	static double cdash[SPLINE_MAX_INPUT]{ 0 };
	static double ddash[SPLINE_MAX_INPUT]{ 0 };

	cdash[0] = c[0] / b[0];
	ddash[0] = d[0] / b[0];

	for (int i = 1; i < n; i++)
	{
		cdash[i] = c[i] / (b[i] - (a[i] * cdash[i - 1]));
		ddash[i] = (d[i] - (a[i] * ddash[i - 1])) / (b[i] - (a[i] * cdash[i - 1]));
	}

	x[n - 1] = ddash[n - 1];

	for (int i = n - 2; i > -1; i--)
	{
		x[i] = ddash[i] - (cdash[i] * x[i + 1]);
	}
}

void splineInterpolation(const double x[], const double y[], int n, const double inpx[], double inpy[], int inpn)
{
	static double h[SPLINE_MAX_INPUT]{ 0 };
	static double d[SPLINE_MAX_INPUT]{ 0 };
	static double l[SPLINE_MAX_INPUT]{ 0 };
	static double u[SPLINE_MAX_INPUT]{ 0 };
	static double w[SPLINE_MAX_INPUT]{ 0 };
	static double v[SPLINE_MAX_INPUT]{ 0 };

	for (int i = 0; i <= n - 1; i++)
	{
		h[i] = x[i + 1] - x[i];
	}

	// 左辺作成
	for (int i = 0; i < n - 1; i++)
	{
		d[i] = 2 * (h[i] + h[i + 1]);

		if (i == 0)
			l[i] = 0;
		else
			l[i] = h[i];

		if (i == n - 2)
			u[i] = 0;
		else
			u[i] = h[i + 1];
	}

	// 右辺作成
	for (int i = 1; i < n - 1; i++)
	{
		w[i - 1] = 6 * (((y[i + 1] - y[i]) / h[i]) - ((y[i] - y[i - 1]) / h[i - 1]));
	}

	// 連立方程式解く
	triDiagonalMatrixSolve(l, d, u, w, v, n - 2);

	// 補間値算出
	int j = 0;
	for (int i = 0; i < n - 1; i++)
	{
		double vi = (i == 0) ? 0 : v[i - 1];
		double vip1 = (i == n - 2) ? 0 : v[i];
		double bi = vi / 2.0;
		double ai = (vip1 - vi) / (6.0 * (x[i + 1] - x[i]));
		double di = y[i];
		double ci = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (x[i + 1] - x[i]) * (2.0 * vi + vip1) / 6.0;

		printf("a, b, c, d = %.3lf, %.3lf, %.3lf, %.3lf\n", di, ci, bi, ai);

		while (inpx[j] < x[i + 1])
		{
			double xdiff = inpx[j] - x[i];
			inpy[j] = ((ai * xdiff + bi) * xdiff + ci) * xdiff + di;
			j++;
		}
	}
}