#include <iostream>
#include <Eigen>
#include <cmath> 
#include "filesAndFunctions.h"
#include "interpolation.h"

typedef double (*functions)(double);

using namespace std;
using namespace Eigen;

int main()
{
    setlocale(LC_ALL, "RUSSIAN");

    //functions f = choiceFunction();
    functions f;
    VectorXd x, func, tabX, tabFunc;
    MatrixXd coeff;

    
    f = f1;
    cout << "y = x^2" << endl;
    tabulationFunction(f, x, func, 0, 3, -1, 1);
    cout << endl << x << endl;
    tabulationLagrange(x, func, tabX, tabFunc, 100, -1, 1);
    printToFile("tests//f1_Lagrange.csv", tabX, tabFunc);
    splineCoeffCalculation(x, func, coeff, 1);
    tabulationSpline(x, coeff, tabX, tabFunc, 100, -1, 1);
    printToFile("tests//f1_Spline.csv", tabX, tabFunc);
    
    f = f2;
    cout << "y = 1/(1+x^2)" << endl;
    tabulationFunction(f, x, func, 1, 11, -5, 5);
    tabulationLagrange(x, func, tabX, tabFunc, 500, -5, 5);
    printToFile("tests//f2_Lagrange.csv", tabX, tabFunc);
    splineCoeffCalculation(x, func, coeff, 1);
    tabulationSpline(x, coeff, tabX, tabFunc, 100, -5, 5);
    printToFile("tests//f2_Spline.csv", tabX, tabFunc);

    f = f3;
    cout << "y = 1/acrtg(1+10x^2)" << endl;
    tabulationFunction(f, x, func, 1, 5, -1, 1);
    tabulationLagrange(x, func, tabX, tabFunc, 300, -1, 1);
    printToFile("tests//f3_Lagrange.csv", tabX, tabFunc);
    splineCoeffCalculation(x, func, coeff, 1);
    tabulationSpline(x, coeff, tabX, tabFunc, 300, -3, 3);
    printToFile("tests//f3_Spline.csv", tabX, tabFunc);

    f = f4;
    cout << "страшная функция" << endl;
    tabulationFunction(f, x, func, 1, 11, -1, 1);
    tabulationLagrange(x, func, tabX, tabFunc, 100, -1, 1);
    printToFile("tests//f4_Lagrange.csv", tabX, tabFunc);
    splineCoeffCalculation(x, func, coeff, 1);
    tabulationSpline(x, coeff, tabX, tabFunc, 100, -1, 1);
    printToFile("tests//f4_Spline.csv", tabX, tabFunc);

    f = f5;
    cout << "y = exp(x)" << endl;
    tabulationFunction(f, x, func, 0, 11, 0, 2);
    tabulationLagrange(x, func, tabX, tabFunc, 100, 0, 2.5);
    printToFile("tests//f5_Lagrange.csv", tabX, tabFunc);
    splineCoeffCalculation(x, func, coeff, 0);
    tabulationSpline(x, coeff, tabX, tabFunc, 100, 0, 2.5);
    printToFile("tests//f5_Spline.csv", tabX, tabFunc);

    cout << lagrangePolynomial(x, func, 2.2) - exp(2.2) << endl;

    cout << "Выполнено!";

    cin.get();
    cin.get();

    return 0;
}