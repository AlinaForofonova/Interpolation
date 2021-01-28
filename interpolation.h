#include <iostream>
#include <Eigen>

typedef double (*functions)(double);

using namespace std;
using namespace Eigen;

//��������� �������, cheb = 0 - �����������, cheb = 1 - �����������
void tabulationFunction(functions &f, VectorXd &x, VectorXd &func, bool cheb = 1, int n = 0, double a = 0, double b = 0);

//���������� �������� �������� � �����
double lagrangePolynomial(VectorXd &x, VectorXd &func, double x0);

//���������� ������������� �������
void splineCoeffCalculation(VectorXd &x, VectorXd &func, MatrixXd &coeff, bool cheb = 1);

//���������� ������� � �����
double spline(VectorXd &x, MatrixXd &coeff, double x0);

//��������� �������� ��������
void tabulationLagrange(VectorXd &xTrue, VectorXd &funcTrue, VectorXd &x, VectorXd &func, int n = 0, double a = 0, double b = 0);

//��������� �������
void tabulationSpline(VectorXd &xTrue, MatrixXd &coeff, VectorXd &x, VectorXd &func, int n = 0, double a = 0, double b = 0);