#pragma once

class direct_problem
{
public:
	direct_problem();
	~direct_problem();

    //!< Вспомогательные переменные
    static bool is_print;
    static bool is_save;
    static const char* path;

    //!< Параметры решения 
    static double step_lambda;
    static double epsilon;
    static double epsilon_integral;
    
    //!< Потенциал, при котором решается задача
    double (*q)(double);

    //!< Функция решения прямой задачи
    void solve_dslp(double* lambdas, double* alphas, int size);

private:
    double epsilon_runge;
	int runge_grid_size = 2;
    double constant = 0;
    double constant_ = 0;

    //!< Исходная система дифференциальных уравнений
    double f1(double, double, double, double);
    double f2(double, double, double, double);

    //!< Классический метод Рунге-Кутта 4-го порядка
    double runge_kutta(double, double*, double*, int);
    //!< Метод Ньютона
    double newton_integral(double*, int);

    //!< Правила Рунге для интеграла и метода Р-К
    double integral_runge(double);
    int rk_runge(double);

    //!< Метод бисекции для вычисления корней уравнения
    void bisection(double*, double*, int);
    
    //!< Вспомогательные функции
    void print(double*, double*, int);
    void save(double*, double*, int);
};