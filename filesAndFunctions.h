#include <iostream>
#include <fstream>
#include <Eigen>

typedef double (*functions)(double);

using namespace std;
using namespace Eigen;

double f1(double x);
double f2(double x);
double f3(double x);
double f4(double x);
double f5(double x);

//����� ������� �� ������
functions choiceFunction();

//������ �������� � ����
void printToFile(string fileName, VectorXd &x, VectorXd &func);