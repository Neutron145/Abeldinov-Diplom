# Решатель прямой и обратной задач Штурма-Лиувилля

Программа позволяет численно решать прямые и обратные задачи Штурма-Лиувилля с краевыми условиями Дирихле.
Данная программа также распространяется в виде динамической библиотеки для удобного подключения к различным проектам.  
Отсутсвуют зависимости от различных библиотек, собственная реализация всех математических методов.  

## Возможности
- Поиск заданного количества спектральных данных для произвольного потенциала.
- Построение произвольного потенциала по заданным спектральным данным и модельному потенциалу.
- Задание контроля точности и сохранение результатов в файл.


## Пример использования
```cpp
#define _USE_MATH_DEFINES
#include "math.h"
#include "stdio.h"
#include "time.h"
#include "direct_problem.h"
#include "inverse_problem.h"

void direct_solve_test() {
  	direct_problem::is_print = true;
  	// Объявляем экзмеляр уравнения, задаем ему потенциал
  	direct_problem standart;
	standart.q = [](double x) -> double {
		return 0; 
	};

	// Получаем для данного уравнения спектральные данные
	int size = 15;
	double* lambdas = new double[size];
	double* alphas = new double[size];
	standart.solve_dslp(lambdas, alphas, size);
}

void inverse_solve_test() {
  	direct_problem::is_print = false;
  	// Задаем уравнения, для которого будем решать обратную задачу
  	direct_problem first;
	first.q = [](double x) -> double {
		return cos(x); 
	};

  	// Получаем для него спектральные данные
	int size = 15;
	double* lambdas = new double[size];
	double* alphas = new double[size];
	first.solve_dslp(lambdas, alphas, size);

  	// Объявляем экземпляр дифференциального уравнения и задаем модельный потенциал
	inverse_problem second;
	second.q_model = [](double x) -> double {
		return 1 - ((2 * x) / M_PI);
	};

  	// Решаем обратную задачу на основе модельного потенциала и спектральных данных
	double* q_ = new double[128];
	second.solve_islp(lambdas, alphas, size, q_, 128);

	for(int i = 0; i < 128; i++) {
		printf("%f\n", q_[i]);
	}
}

int main() {
	direct_solve_test();
	inverse_solve_test();
}
```

---
Абельдинов Рафаэль  
2025
