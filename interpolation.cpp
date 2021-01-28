#include <iostream>
#include <Eigen>
#include <cmath>

typedef double (*functions)(double);

using namespace std;
using namespace Eigen;

//табуляция функции, cheb = 0 - равномерная, cheb = 1 - чебышевская
void tabulationFunction(functions &f, VectorXd &x, VectorXd &func, bool cheb = 1, int n = 0, double a = 0, double b = 0)
{
    if (n == 0)
    {
        cout << "Введите количество узлов: ";
        cin >> n;
        cout << "Введите начало и конец: ";
        cin >> a >> b;
    }

    x.resize(n);
    func.resize(n);

    if(cheb)
    {
        const double pi = 3.14159265358979323846;
        double hp = 0.5*(a+b), hm = 0.5*(b-a), N = 2*n;
        for(int i = n; i > 0; i--)
        {
            x(n-i) = hp + hm*cos( pi*(2*i-1)/N );
            func(n-i) = f(x(n-i));
        }
    }
    else
    {
        double h = (b-a)/(n-1);
        for(int i = 0; i < n; i++)
        {
            x(i) = a+i*h;
            func(i) = f(x(i));
        }
    }
}



//вычисление полинома Лагранжа в точке
double lagrangePolynomial(VectorXd &x, VectorXd &func, double x0)
{
    double result = 0.0;
    double product;

    for(int i = 0; i < x.size(); i++)
    {
        product = 1.0;
        for(int j = 0; j < x.size(); j++)
            if(i != j) 
                product *= ( x0-x(j) ) / ( x(i)-x(j) );
        result += product * func(i);
    }

    return result;
}

//ПРОМЕЖУТОЧНЫЕ ФУНКЦИИ ТРЕБУЮЩИЕСЯ ДЛЯ ВЫЧИСЛЕНИЯ КОЭФФИЦИЕНТОВ СПЛАЙНА
//h(i)
double h(VectorXd &x, double i) { return x(i)-x(i-1); }
//2*(h(i-1) + h(i))
double H(VectorXd &x, double i) { return 2*(x(i)-x(i-2)); }
//3*(g(i+1) - g(i))
double G(VectorXd &func, VectorXd &x, double i) { return 3*( h(func,i)/h(x,i) - h(func,i-1)/h(x,i-1) ); }
//3h*(g(i) - g(i-1)) (для равномерной сетки)
double Gh(VectorXd &func, double i) { return 3*( func(i)-2*func(i-1)+func(i-2) ); }

//Вычисление коэффициентов сплайна, cheb = 0 - равномерная, cheb = 1 - чебышевская
void splineCoeffCalculation(VectorXd &x, VectorXd &func, MatrixXd &coeff, bool cheb = 1)
{
    int n = x.size();
    VectorXd alpha(n), beta(n);

    
    alpha(1) = beta(1) = 0;
    if (cheb) //чебышевская сетка
    {
        double HH; // HH(i) = H(i) - h(i-1)*beta(i-1)
        for(int i = 2; i < n; i++)
        {
            HH = H(x,i) - h(x,i-1)*beta(i-1);
            alpha(i) = ( G(func,x,i) - h(x,i-1)*alpha(i-1) ) / HH;
            beta(i) = h(x,i) / HH;
        }
    }
    else //равномерная сетка
    {
        double h = pow(x(1)-x(0), 2); //h^2
        for(int i = 2; i < n; i++)
        {
            beta(i)  = 1/(4-beta(i-1));
            alpha(i) = (Gh(func,i)/h - alpha(i-1)) * beta(i);
        }
    }
    

    //0 - a, 1 - b, 2 - c, 3 - d
    coeff.resize(4, n);
    
    for(int i = 1; i < n; i++)
        coeff(0,i) = func(i-1); 

    coeff(2,n-1) = alpha(n-1);
    coeff(1,n-1) = h(func,n-1)/h(x,n-1) - 2*coeff(2,n-1)*h(x,n-1)/3;
    coeff(3,n-1) = -coeff(2,n-1) / (3*h(x,n-1));
    for(int i = n-2; i > 1; i--)
    {
        coeff(2,i) = alpha(i) - beta(i)*coeff(2,i+1);
        coeff(1,i) = h(func,i)/h(x,i) - (coeff(2,i+1)+2*coeff(2,i)) * h(x,i)/3;
        coeff(3,i) = (coeff(2,i+1)-coeff(2,i)) / (3*h(x,i));
    }
    coeff(2,1) = 0;
    coeff(1,1) = h(func,1)/h(x,1) - coeff(2,2)*h(x,1)/3;
    coeff(3,1) = coeff(2,2) / (3*h(x,1));
}

//Вычисление сплайна в точке
double spline(VectorXd &x, MatrixXd &coeff, double x0)
{
    int i = 0;
    while(x(i) < x0 && i < x.size()-1)
        i++;
    
    double h = x0;
    if (i == 0)
        {
            h -= x(i);
            i++;
        }
    else
        h -= x(i-1);

    return coeff(0,i) + coeff(1,i)*h + coeff(2,i)*pow(h,2) + coeff(3,i)*pow(h,3);
}


//Табуляция полинома Лагранжа
void tabulationLagrange(VectorXd &xTrue, VectorXd &funcTrue, VectorXd &x, VectorXd &func, int n = 0, double a = 0, double b = 0)
{
    if(n == 0)
    {
        cout << "Введите количество точек: ";
        cin >> n;
        cout << "Введите начало и конец: ";
        cin >> a >> b;
    }

    x.resize(n);
    func.resize(n);

    double h = (b-a)/(n-1);
    for(int i = 0; i < n; i++)
    {
        x(i) = a + i*h;
        func(i) = lagrangePolynomial(xTrue, funcTrue, x(i));
    }
}

//Табуляция сплайна
void tabulationSpline(VectorXd &xTrue, MatrixXd &coeff, VectorXd &x, VectorXd &func, int n = 0, double a = 0, double b = 0)
{
    if (n == 0)
    {
        cout << "Введите количество точек: ";
        cin >> n;
        cout << "Введите начало и конец: ";
        cin >> a >> b;
    }

    x.resize(n);
    func.resize(n);

    double h = (b-a)/(n-1);
    for(int i = 0; i < n; i++)
    {
        x(i) = a + i*h;
        func(i) = spline(xTrue, coeff, x(i));
    }
}