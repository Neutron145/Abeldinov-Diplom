#define _USE_MATH_DEFINES
#include "math.h"
#include <iostream>
#include <fstream>
#include "time.h"
#include "inverse_problem.h"

inverse_problem::inverse_problem(){

}

inverse_problem::~inverse_problem(){

}

//!< Точность вычисления потенциала
double inverse_problem::epsilon = 1e-6;

//!< Система уравнений первого порядка
double inverse_problem::f1(double u) { return u; }
double inverse_problem::f2(double x, double y, double lambda) { return x * y - lambda * y; }

//!< Функция загрузки спектральных данных
int inverse_problem::load_spectral_data(double** lambdas, double** alphas, const char* path) {
	std::ifstream fin(path);
	//!< Количество спектральных параметров
	int p_size;
	fin >> p_size;
	*lambdas = new double[p_size];
	*alphas = new double[p_size];
	//!< Считываем спектральные данные из файла
	for (int i = 0; i < p_size; i++)
		fin >> (*lambdas)[i] >> (*alphas)[i];
	return p_size;
}

//!< Функция вычисления решения методом Рунге-Кутта 4-го порядка 
void inverse_problem::runge_kutta(double lambda, double* q, double* y_arr, double* u_arr, int size) {
	double y = 0, u = 1;
	double h = M_PI / (double)(size - 1);
	int i = 0;
	for (double x = 0; x < M_PI; x += h) {
		y_arr[i] = y;
		u_arr[i] = u;

		double k1 = h * f1(u);
		double l1 = h * f2(q[i], y, lambda);
		double k2 = h * f1(u + l1 / 2);
		double l2 = h * f2((q[i] + q[i + 1]) / 2, y + k1 / 2, lambda);
		double k3 = h * f1(u + l2 / 2);
		double l3 = h * f2((q[i] + q[i + 1]) / 2, y + k2 / 2, lambda);
		double k4 = h * f1(u + l3);
		double l4 = h * f2(q[i + 1], y + k3, lambda);

		y += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		u += (l1 + 2 * l2 + 2 * l3 + l4) / 6;
		i++;
	}
}

void inverse_problem::q_model_discrete(double* q, int size) {
	double x = 0;
	for (int i = 0; i < size; x += M_PI / (size - 1)) { q[i++] =  q_model(x);  }
}

void inverse_problem::solve_islp(double* lambdas, double* alphas, int param_count, double* q, int nodes_count) {
	spectral_mappings(lambdas, alphas, param_count, q, nodes_count);
}

void inverse_problem::solve_islp(const char* path, int param_count, double* q, int nodes_count) {
	double* lambdas;
	double* alphas;
	int size = load_spectral_data(&lambdas, &alphas, path);
	spectral_mappings(lambdas, alphas, param_count, q, nodes_count);
}

//!< Метод спектральных отображений
void inverse_problem::spectral_mappings(double* lambdas, double* alphas, int param_count, double* q, int nodes_count) {
	//!< Получение спектральных данных для модельного потенциала 
	double* tq_lambdas = new double[param_count];
	double* tq_alphas = new double[param_count];
	direct_solver.q = q_model;
	direct_solver.solve_dslp(tq_lambdas, tq_alphas, param_count);


	//!< Первое приближение потенциала
	double* q0 = new double[nodes_count];
	q_model_discrete(q0, nodes_count);
	//!< Значения потенциала на текущем и предыдущем итерациях
	double* q1 = new double[nodes_count];
	double* q2 = new double[nodes_count];
	q_model_discrete(q2, nodes_count);

	//!< Массивы для хранения всех решений S и их производных
	double* y1 = new double[nodes_count];
	double* y2 = new double[nodes_count];
	double* y3 = new double[nodes_count];
	double* y4 = new double[nodes_count];
	double* u1 = new double[nodes_count];
	double* u2 = new double[nodes_count];
	double* u3 = new double[nodes_count];
	double* u4 = new double[nodes_count];

	bool is_err = true;
	int r = 1;
	printf("\nepsilon: %f\n", epsilon);
	double max_err = -1;
	while (fabs(max_err) > 1e-6) {
		double* eps = new double[nodes_count];
		for (int i = 0; i < nodes_count; i++) {
			q1[i] = q2[i];
			eps[i] = 0;
		}
		//!< Вычисление eps' во всех узлах x_i
		for (int k = 0; k < param_count; k++) {
			runge_kutta(lambdas[k], q1, y1, u1, nodes_count);
			runge_kutta(lambdas[k], q0, y2, u2, nodes_count);
			runge_kutta(tq_lambdas[k], q1, y3, u3, nodes_count);
			runge_kutta(tq_lambdas[k], q0, y4, u4, nodes_count);
			for (int i = 0; i < nodes_count; i++) {
				double s1 = u1[i] * y2[i] + y1[i] * u2[i];
				double s2 = u3[i] * y4[i] + y3[i] * u4[i];
				eps[i] += (s1 / alphas[k]) - (s2 / tq_alphas[k]);
			}
		}
		//!< Вычисление q_r во всех узлах x_i
		for (int i = 0; i < nodes_count; i++) {
			q2[i] = q0[i] - 2 * eps[i];
		}

		//!< Вычисление ошибки для условия остановки
		max_err = -1;
		for (int i = 0; i < nodes_count; i++) {
			double err = fabs(q1[i] - q2[i]);
			if (err > max_err) {
				max_err = err;
			}
		}

		printf("r: %d max error: %f\n", r++, max_err);
	}
	
	for (int i = 0; i < nodes_count; i++) {
		q[i] = q2[i];
	}
}