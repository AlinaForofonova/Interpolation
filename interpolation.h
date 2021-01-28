#include <iostream>
#include <Eigen>

typedef double (*functions)(double);

using namespace std;
using namespace Eigen;

//табуляция функции, cheb = 0 - равномерная, cheb = 1 - чебышевская
void tabulationFunction(functions &f, VectorXd &x, VectorXd &func, bool cheb = 1, int n = 0, double a = 0, double b = 0);

//вычисление полинома Лагранжа в точке
double lagrangePolynomial(VectorXd &x, VectorXd &func, double x0);

//Вычисление коэффициентов сплайна
void splineCoeffCalculation(VectorXd &x, VectorXd &func, MatrixXd &coeff, bool cheb = 1);

//Вычисление сплайна в точке
double spline(VectorXd &x, MatrixXd &coeff, double x0);

//Табуляция полинома Лагранжа
void tabulationLagrange(VectorXd &xTrue, VectorXd &funcTrue, VectorXd &x, VectorXd &func, int n = 0, double a = 0, double b = 0);

//Табуляция сплайна
void tabulationSpline(VectorXd &xTrue, MatrixXd &coeff, VectorXd &x, VectorXd &func, int n = 0, double a = 0, double b = 0);