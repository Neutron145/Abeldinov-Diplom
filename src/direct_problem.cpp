#define _USE_MATH_DEFINES
#include "math.h"
#include "stdio.h"
#include "time.h"
#include "direct_problem.h"

bool direct_problem::is_print = false;
bool direct_problem::is_save = false;
const char* direct_problem::path = 0;

double direct_problem::step_lambda = 0.1;
double direct_problem::epsilon = 1e-5;
double direct_problem::epsilon_integral = 1e-5;

direct_problem::direct_problem()
{
	q = nullptr;
}

direct_problem::~direct_problem()
{

}

double direct_problem::f1(double x, double u1, double u2, double lambda) { return u2; }
double direct_problem::f2(double x, double u1, double u2, double lambda) { return (q(x) - constant) * u1 - lambda * u1; }

double direct_problem::runge_kutta(double lambda, double* u1, double* u2, int size) {
	u1[0] = 0;
	u2[0] = 1;
	double h = M_PI / (double)(size - 1);
	double x = 0;
	for (int i = 1; i < size; i++) {
		double k11 = h * f1(x, u1[i - 1], u2[i - 1], lambda);
		double k12 = h * f2(x, u1[i - 1], u2[i - 1], lambda);

		double k21 = h * f1(x + (h / 2), u1[i - 1] + (k11 / 2), u2[i - 1] + (k12 / 2), lambda);
		double k22 = h * f2(x + (h / 2), u1[i - 1] + (k11 / 2), u2[i - 1] + (k12 / 2), lambda);

		double k31 = h * f1(x + (h / 2), u1[i - 1] + (k21 / 2), u2[i - 1] + (k22 / 2), lambda);
		double k32 = h * f2(x + (h / 2), u1[i - 1] + (k21 / 2), u2[i - 1] + (k22 / 2), lambda);

		double k41 = h * f1(x + h, u1[i - 1] + k31, u2[i - 1] + k32, lambda);
		double k42 = h * f2(x + h, u1[i - 1] + k31, u2[i - 1] + k32, lambda);

		u1[i] = u1[i - 1] + (k11 + 2 * k21 + 2 * k31 + k41) / 6;
		u2[i] = u2[i - 1] + (k12 + 2 * k22 + 2 * k32 + k42) / 6;
		x += h;
	}
	return u1[size - 1];
}

double direct_problem::newton_integral(double* values, int length) {
	double h = M_PI / (float(length) - 1.);
	double I = values[0] + values[length - 1];
	double s1 = 0, s2 = 0;
	for (int i = 1; i < length - 1; i++) {
		if (i % 3 == 0) s1 += values[i];
		else s2 += values[i];
	}
	I += 2 * s1 + 3 * s2;
	I *= (3 * h) / 8;
	return I;
}

double direct_problem::integral_runge(double lambda) {
	int length = 4;
	double h = M_PI / (double)(length - 1);
	double* y1 = new double[length];
	double* y1_ = new double[length];
	runge_kutta(lambda, y1, y1_, length);
	for(int i = 0; i < length; i++) {
		y1[i] = y1[i] * y1[i];
	}
	double I2 = newton_integral(y1, length);
	double I1 = 0;
	int i = 0;
	do {
		I1 = I2;
		length = length * 2 - 1;
		h /= 2;
		double* y2 = new double[length];
		double* y2_ = new double[length];
		runge_kutta(lambda, y2, y2_, length);
		for(int i = 0; i < length; i++) {
			y2[i] = y2[i] * y2[i];
		}
		I2 = newton_integral(y2, length);
		i++;
	} while (fabs(I1 - I2) > epsilon_integral && i < 15);
	return I2;
}

int direct_problem::rk_runge(double lambda) { 
	//!< Количество узлов
	int nodes_count = 2;
	double h = M_PI / (2. - 1);
	double* y1 = new double[nodes_count];
	double* y1_ = new double[nodes_count];
	double* errors = new double[20];
	double* hs = new double[20];
	double ans1 = 0, ans2 = 0;
	//!< Вычисляем S при максимально возможном шаге
	ans1 = runge_kutta(lambda, y1, y1_, nodes_count);
	double* y2 = new double[2 * nodes_count - 1];
	double* y2_ = new double[2 * nodes_count - 1];
	double max_err = -1;
	int iter = 0;
	do {
		//!< Вычисляем S на сетке с h/2
		ans2 = runge_kutta(lambda, y2, y2_, 2 * nodes_count - 1);
		//!< Определяем максимальную ошибку на двух сетках
		max_err = -1;
		for (int i = 0; i < nodes_count; i++) {
			double err = fabs(y1[i] - y2[i * 2]);
			if (err > max_err) max_err = err;
		}
		//!< Если ошибка большая, то h = h/2
		if (max_err > 15 * epsilon_runge) {
			ans1 = ans2;
			nodes_count = 2 * nodes_count - 1;
			h /= 2;
			delete[] y1;
			delete[] y1_;
			y1 = new double[nodes_count];
			y1_ = new double[nodes_count];
			for (int i = 0; i < nodes_count; i++) {
				y1[i] = y2[i];
			}
			delete[] y2;
			delete[] y2_;
			y2 = new double[2 * nodes_count - 1];
			y2_ = new double[2 * nodes_count - 1];
			errors[iter] = max_err;
			hs[iter] = h;
			iter++;
		}
	} while (max_err > 15 * epsilon_runge && iter < 20);
	
	delete[] y1;
	delete[] y1_;
	delete[] y2;
	delete[] y2_;
	return nodes_count;
}

void direct_problem::bisection(double* lambdas, double* alphas, int size) {
	double l_l = 0;
	double l_r = 0;
	//!< Количество найденных корней
	int n = 0;	
	double* y = new double[runge_grid_size];
	double* y_ = new double[runge_grid_size];
	double lambda = -1;

	clock_t start, end;
	start = clock();
	while(n < size) {
		l_l = l_r;
		l_r = runge_kutta(lambda, y, y_, runge_grid_size);
		if (l_l * l_r < 0) {
			double l = lambda - step_lambda, r = lambda;
			double mid = (l + r) / 2;
			do {
				double l_m = runge_kutta(mid, y, y_, runge_grid_size);
				if (l_m * l_l < 0) {
					r = mid;
					l_r = l_m;
				}
				else {
					l = mid;
					l_l = l_m;
				}
				mid = (l + r) / 2;
			} while(fabs(r - l) > epsilon); 
			end = clock();
			start = end;	
			lambdas[n] = mid + constant;
			constant = 0;
			alphas[n] = integral_runge(lambdas[n]);
			constant = constant_;
			if(n % 15 == 0 && n != 0) {
				runge_grid_size = runge_grid_size * 2 - 1;
				y = new double[runge_grid_size];
				y_ = new double[runge_grid_size];
			}
			n++;
		}
		lambda += step_lambda;
	}
}

void direct_problem::print(double* lambdas, double* alphas, int size) {
	printf("lambda\talpha\n");
	for (int i = 0; i < size; i++) {
		printf("%0.4f\t%0.4f\n", lambdas[i], alphas[i]);
	}
}

void direct_problem::save(double* lambdas, double* alphas, int size) {
	FILE* fptr;
	int code = fopen_s(&fptr, path, "w");
	if (fptr != NULL && code == 0) {
		fprintf(fptr, "%d ", size);
		for (int i = 0; i < size; i++) {
			fprintf(fptr, "%e %e ", lambdas[i], alphas[i]);
		}
	}
	if (fptr != NULL) { fclose(fptr); }
}

void direct_problem::solve_dslp(double* lambdas, double* alphas, int size) {
	//!< Вычисление погрешности для метода Рунге-Кутта
	double eps_coef = fabs(M_PI * 15 * cos(M_PI * 15) - sin(M_PI * 15)) / (2 * pow(15, 3));
	epsilon_runge = epsilon * eps_coef;

	//!< Размер сетки для метода Рунге-Кутта 
	runge_grid_size = rk_runge(pow(15, 2));

	//!< Поиск минимального значения
	for(double x = 0; x < M_PI; x += M_PI/(runge_grid_size - 1)) {
		if (q(x) < constant) constant = q(x);
	}
	constant_ = constant;

	bisection(lambdas, alphas, size);
	if (path != 0 && is_save) { save(lambdas, alphas, size); }
	if (is_print) { 
		print(lambdas, alphas, size);
	}
}