#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <math.h> 

double f1(double x1, double x2)
{
    double result;

    result = x1 * x1 - x1 + x2 * x2 - 1;

    return result;
}

double f2(double x1, double x2)
{
    double result;

    result = x2 - tan(x1);

    return result;
}

double df1dx1(double x1)
{
    double result;

    result = 2 * x1 - 1.0;

    return result;
}

double df1dx2(double x2)
{
    double result;

    result = 2 * x2;

    return result;
}

double df2dx1(double x1)
{
    double result;

    result = -1 / (cos(x1) * cos(x1));

    return result;
}

double df2dx2(double x2)
{
    double result;

    result = 1.0;

    return result;
}

double detA1(double x1, double x2)
{
    double det;

    det = f1(x1, x2) - df1dx2(x2) * f2(x1, x2);

    return det;
}

double detA2(double x1, double x2)
{
    double det;

    det = df1dx1(x1) * f2(x1, x2) - f1(x1, x2) * df2dx1(x1);

    return det;
}

double detJ(double x1, double x2)
{
    double det;

    det = df1dx1(x1) * df2dx2(x2) - df2dx1(x1);

    return det;
}

double phi1(double x2)
{
    double result;

    result = atan(x2);

    return result;
}

double phi2(double x1)
{
    double result;

    result = sqrt(1 - x1 * x1 + x1);
    
    return result;
}

void PrintRoots(std::vector<double> roots)
{
    for(int i = 0; i < roots.size(); i++)
    {  
        printf("\nX%d = %.4f", i + 1, roots[i]);
    }
}

double MaxNorm(std::vector<double> error)
{
    double max = error[0];

    if(error[1] > max)
    {
        max = error[1];
    }

    return max;
}

std::vector<double> NewtonMethod(double x1, double x2, double eps)
{
    std::vector<double> roots(2, 0);
    std::vector<double> error(2, 0);
    int i = 0;

    for(;;)
    {
        roots[0] = x1 - detA1(x1, x2) / detJ(x1, x2);
        roots[1] = x2 - detA2(x1, x2) / detJ(x1, x2);

        error[0] = abs(x1 - roots[0]);
        error[1] = abs(x2 - roots[1]);

        if(MaxNorm(error) <= eps)
        {
            printf("\n Number of iterations: %d", i + 1);
            break;
        }
        x1 = roots[0];
        x2 = roots[1];
        i++;
    }

    return roots;
}

std::vector<double> IterationMethod(double x1, double x2, double eps)
{
    std::vector<double> roots(2, 0);
    std::vector<double> error(2, 0);
    int i = 0;

    for(;;)
    {
        double q = 0.5;
        roots[0] = phi1(x2);
        roots[1] = phi2(x1);

        error[0] = abs(x1 - roots[0]);
        error[1] = abs(x2 - roots[1]);

        if(q / (1 - q) * MaxNorm(error) < eps)
        {
            printf("\n Number of iterations: %d", i + 1);
            break;
        }

        x1 = roots[0];
        x2 = roots[1];
        i++;
    }

    return roots;
}

int main()
{
    double x1, x2, eps;
    std::cout << "\nThe system is x^2 - x + y^2 - 1 = 0\n              y - tg(x) = 0\n\n";

    // Read starting values
    //-------------
    std::cout << "Please, enter the starting point\n";
    std::cin >> x1;
    std::cin >> x2;

    // Read epsilon
    //-------------
    std::cout << "Please, enter the epsilon value\n";
    std::cin >> eps;

    // Find the roots of the equation
    //-------------------------------
    std::cout << "\nThe solution using Newton's Method is: ";
    PrintRoots(NewtonMethod(x1, x2, eps));
    std::cout << "\nThe solution using Iteration Method is: ";
    PrintRoots(IterationMethod(x1, x2, eps));
}

// Start: 1.5 1.5