#pragma once

#include <limits>

class cubic_spline
{
private:
	struct spline_tuple
	{
		double a, b, c, d, x;
	};

	spline_tuple *splines; // Сплайн
	int n; // Количество узлов сетки

	void free_mem(); // Освобождение памяти

public:
	cubic_spline();
	~cubic_spline();

					 // Построение сплайна
					 // x - узлы сетки, должны быть упорядочены по возрастанию, кратные узлы запрещены
					 // y - значения функции в узлах сетки
					 // n - количество узлов сетки
	void build_spline(const double *x, const double *y, int n);

	// Вычисление значения интерполированной функции в произвольной точке
	double f(double x) const;
};



