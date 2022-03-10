#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <math.h> 

double Equation(double x)
{
    double result;

    result = powl(10, x) - 5 * x - 2;

    return result;
}

double Derivative(double x)
{
    double result;

    result = log(10) * powl(10, x) - 5;

    return result;
}

double EqEquation(double x)
{
    double result;

    result = log10(5 * x + 2);

    return result;
}

double DerEqEquation(double x)
{
    double result;

    result = 5 / ((5 * x + 2) * log(10));

    return result;
}

double NewtonForLinEq(double start_point, double eps)
{
    double solution = eps, temp;

    for(;;)
    {
        solution = start_point - Equation(start_point) / Derivative(start_point);
        if(abs(start_point - solution) <= eps)
        {
            break;
        }
        start_point = solution;
    }

    return solution;
}

double IterationMethod(double start_point, double eps)
{
    double solution;
    double q = DerEqEquation(start_point);

    for(;;)
    {
        solution = EqEquation(start_point);
        if(q / (1 - q) * abs(start_point - solution) <= eps)
        {
            break;
        }
        start_point = solution;
    }

    return solution;
}

int main()
{
    double x, eps;
    std::cout << "\nThe equation is 10^x - 5x - 2\n";
    // Read starting values
    //-------------
    std::cout << "Please, enter the starting value\n";
    std::cin >> x;
    // Read epsilon
    //-------------
    std::cout << "Please, enter the epsilon value\n";
    std::cin >> eps;

    printf("\nThe solution using Newton's Method is: %.4f", NewtonForLinEq(x, eps));
    printf("\nThe solution using Iteration Method is: %.4f", IterationMethod(x, eps));
}