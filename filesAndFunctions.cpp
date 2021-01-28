#include <iostream>
#include <fstream>
#include <math.h>
#include <Eigen>

typedef double (*functions)(double);

using namespace std;
using namespace Eigen;

double f1(double x) { return x*x; }
double f2(double x) { return 1/(1+x*x); }
double f3(double x) { return 1/atan(1+10*x*x); }
double f4(double x) { return pow(4*x*x*x+2*x*x-4*x+2, sqrt(2)) + asin(1/(5+x-x*x)) - 5; }
double f5(double x) { return exp(x); }

//выбор функции из списка
functions choiceFunction()
{
    cout << "Выберите функцию:" << endl
         << "1. y = x^2" << endl
         << "2. y = 1/(1+x^2)" << endl
         << "3. y = 1/arctg(1+10x^2)" << endl
         << "4. страшная функция" << endl;

    int num;
    cin >> num;

    switch(num)
    {
        case 1: return f1;
        case 2: return f2;
        case 3: return f3;
        case 4: return f4;
    }

    return f1;
}

//печать значений в файл
void printToFile(string fileName, VectorXd &x, VectorXd &func)
{
    ofstream fout(fileName);

    for(int i = 0; i < x.size(); i++)
    {
        fout << x(i) << ",";
        fout << func(i) << endl;
    }

    fout.close();
}

