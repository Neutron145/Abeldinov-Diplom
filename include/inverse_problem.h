#include "direct_problem.h"

class inverse_problem
{
public:
	inverse_problem();
	~inverse_problem();

	//!< Точность вычислений
	static double epsilon;

	//!< Модельный потенциал (непрерывный)
	double (*q_model)(double);

	void solve_islp(double*, double*, int, double*, int);
	void solve_islp(const char*, int, double*, int);
private:
	direct_problem direct_solver;

	//!< Исходная система дифференциальных уравнений 
	double f1(double);
	double f2(double, double, double);

	//!< Загрузка спектральных данных из файла
	int load_spectral_data(double**, double**, const char*);

	//!< Классический метод Рунге-Кутта 4-го порядка
	void runge_kutta(double, double*, double*, double*, int);

	//!< Дискретизированный модельный потенциал
	void q_model_discrete(double*, int);

	//!< Функция решения обратной задачи
	void spectral_mappings(double* lambdas, double* alphas, int param_count, double* q, int size);

};