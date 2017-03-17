#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "cubic_spline.h"

#define FIELD_WIDTH 20
#define PRESICION 6

using namespace std;

double origFunc(double x);
double sixDer(double x);

double newtonPolinom(double testX, int count, double ** data);
double lagrangePolinom(double testX, int count, double ** data);
double magoranta(double testX, int count, double ** data);

int factorial(int n)
{
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


int main(int argc, const char * argv[])
{
	char * inFile = "input.txt";
	char * outFile = "output.txt";

	ifstream in(inFile);
	ofstream out(outFile);

	int count;
	in >> count;

	double ** data = new double*[2];
	data[0] = new double[count];
	data[1] = new double[count];

	double ** dataRev = new double*[2];
	dataRev[0] = new double[count];
	dataRev[1] = new double[count];

	for (int i = 0; i < count; i++)
	{
		in >> data[0][i] >> data[1][i];
		dataRev[0][count - 1 - i] = data[0][i];
		dataRev[1][count - 1 - i] = data[1][i];
	}		

	double delta = abs((data[0][1] - data[0][0]) / 10);

	cubic_spline spline;

	double *x = new double[count];
	double *y = new double[count];
	for (int i = 0; i < count; i++)
	{
		x[i] = data[0][i];
		y[i] = data[1][i];
	}

	spline.build_spline(x, y, count);

	out << setw(FIELD_WIDTH) << "x";
	out << setw(FIELD_WIDTH) << "Original func";
	out << setw(FIELD_WIDTH) << "Newton Pol -> fault";
	out << setw(FIELD_WIDTH) << "Newton Pol <- fault";
	out << setw(FIELD_WIDTH) << "Lagrange fault";
	out << setw(FIELD_WIDTH) << "Splines fault";
	out << setw(FIELD_WIDTH) << "Magoranta" << endl;
	for (int i = 0; i < FIELD_WIDTH * 7; i++)
		out << "-";
	out << endl;

	for (double i = data[0][0]; i <= data[0][count-1]+delta; i += delta)
	{
		out << setw(FIELD_WIDTH) << setprecision(PRESICION);
		out << i <<"\t";
		out << setw(FIELD_WIDTH) << setprecision(PRESICION);
		out << origFunc(i) << "\t";
		out << setw(FIELD_WIDTH) << setprecision(PRESICION);
		out << abs(newtonPolinom(i, count, data) - origFunc(i) ) << "\t";
		out << setw(FIELD_WIDTH) << setprecision(PRESICION);
		out << abs(newtonPolinom(i, count, dataRev) - origFunc(i) ) << "\t";
		out << setw(FIELD_WIDTH) << setprecision(PRESICION);
		out << abs(lagrangePolinom(i, count, data) - origFunc(i) ) << "\t";
		out << setw(FIELD_WIDTH) << setprecision(PRESICION);
		out << abs(spline.f(i) - origFunc(i)) << "\t";
		out << setw(FIELD_WIDTH);
		out << magoranta(i,count,data) << endl;
	}

	/*for (double i = data[0][0]; i <= data[0][count - 1]; i += delta)
	{
		out << i << "\t";
		out << origFunc(i) << "\t";
		out << newtonPolinom(i, count, data)<< "\t";
		out << newtonPolinom(i, count, dataRev) << "\t";
		out << lagrangePolinom(i, count, data) << "\t";
		out << spline.f(i) << endl;
	}*/

	delete[] data[0];
	delete[] data[1];
	delete[] data;

	delete[] dataRev[0];
	delete[] dataRev[1];
	delete[] dataRev;

	delete[] x;
	delete[] y;

	in.close();
	out.close();

	return 0;
}

double origFunc(double x)
{
	return log(sin(0.5)) - ((pow(cos(12*x),2.0))/(sin(24*x)*24));
}

double sixDer(double x)
{
	return (-124416.0*(302.0*cos(12.0*x)+57.0*cos(36.0*x)+cos(60.*x)))/(pow(sin(12.*x),7.));
}

double newtonPolinom(double testX, int count, double ** data)
{
	double res = data[1][0], divDiff, den;
	int i, j, k;
	for (i = 1; i < count; i++)
	{
		divDiff = 0;
		for (j = 0; j <= i; j++)
		{
			den = 1;
			for (k = 0; k <= i; k++)
			{
				if (k != j) den *= (data[0][j] - data[0][k]);
			}
			divDiff += data[1][j] / den;
		}
		for (k = 0; k < i; k++)
			divDiff *= (testX - data[0][k]);
		res += divDiff;
	}

	return res;
}

double lagrangePolinom(double testX, int count, double ** data)
{
	double res = 0, L;
	int i, j, k;
	for (int i = 0; i < count; i++)
	{
		L = 1;
		for (int j = 0; j < count; j++)
		{
			if (j != i)
				L *= (testX - data[0][j]) / (data[0][i] - data[0][j]);
		}
		res += data[1][i] * L;
	}

	return res;
}

double magoranta(double testX, int count, double ** data)
{
	double sup = abs(sixDer(data[0][0]));
	double der = 0;
	double mult = 1;
	for (int i = 0; i < count; i++)
	{
		der = abs(sixDer(data[0][i]));
		if (sup < der) sup = der;
		mult *= abs(testX - data[0][i]);
	}

	return (sup*mult)/factorial(count+1);
}

