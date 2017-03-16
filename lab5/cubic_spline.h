#pragma once

#include <limits>

class cubic_spline
{
private:
	struct spline_tuple
	{
		double a, b, c, d, x;
	};

	spline_tuple *splines;
	int n;

	void free_mem();

public:
	cubic_spline();
	~cubic_spline();

	void build_spline(const double *x, const double *y, int n);

	double f(double x) const;
};



