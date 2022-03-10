#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

void NiceVectorPrint(std::vector<double> vector)
{
    printf("\n");
    for(int i = 0; i < vector.size(); i++)
    {
        printf("%.5f ", vector[i]);
    }
    printf("\n");
}

void NiceMatrixPrint(std::vector<std::vector<double>> matrix)
{
    std::vector<std::vector<double>>::iterator row;
    std::vector<double>::iterator col;
    printf("\n");
    for(row = matrix.begin(); row != matrix.end(); row++)
    {
        for (col = row->begin(); col != row->end(); col++)
        {
            printf("%.5f ", *col);
        }
        printf("\n");
    }
    printf("\n");
}

double Function(double x, double y, double z)
{
    double result;

    result = ((2 * x + 4) * z - 2 * y) / (x * (x + 4));

    return result;
}

double ExactSolution(double x)
{
    double result;
    
    result = x * x + x + 2;

    return result;
}

std::vector<double> ExactResult(double a, double b, double h)
{
    std::vector<double> result;

    for(double i = a; i <= b + h; i += h)
    {
        result.push_back(ExactSolution(i));
    }

    return result;
}

std::vector<double> MakeSequence(double a, double b, double h)
{
    std::vector<double> result;

    for(double i = a; i <= b + h; i += h)
    {
        result.push_back(i);
    }

    return result;
}

std::vector<double> RungeKutta(double y, double z, double a, double b, double h)
{
    std::vector<double> y_vec;
    double delta1, delta2;
    double k1, k2, k3, k4, l1, l2, l3, l4;
    std::vector<double> x_vec = MakeSequence(a, b, h);
    std::vector<double> z_vec;
    y_vec.push_back(y);
    z_vec.push_back(z);

    for(int i = 1; i < x_vec.size(); i++)
    {
        l1 = h * Function(x_vec[i - 1], y_vec[i - 1], z_vec[i - 1]);
        k1 = h * z_vec[i - 1];
        l2 = h * Function(x_vec[i - 1] + h / 2, y_vec[i - 1] + k1 / 2, z_vec[i - 1] + l1 / 2);
        k2 = h * (z_vec[i - 1] + l1 / 2);
        l3 = h * Function(x_vec[i - 1] + h / 2, y_vec[i - 1] + k2 / 2, z_vec[i - 1] + l2 / 2);
        k3 = h * (z_vec[i - 1] + l2 / 2);
        l4 = h * Function(x_vec[i - 1] + h, y_vec[i - 1] + k3, z_vec[i - 1] + l3);
        k4 = h * (z_vec[i - 1] + l3);

        // printf("\n %.4f %.4f %.4f %.4f", k1, k2, k3, k4);
        delta1 = (k1  + 2 * k2 + 2 * k3 + k4) / 6;
        delta2 = (l1  + 2 * l2 + 2 * l3 + l4) / 6;

        y_vec.push_back(y_vec[i - 1] + delta1);
        z_vec.push_back(z_vec[i - 1] + delta2);
    }
    
    return y_vec;
}

std::vector<double> ShootingMethod(double a, double b, double h, double y, double y_right, double eps)
{
    std::vector<double> result;
    double eta0 = 1, eta1 = 0.8, temp;
    double f0, f1;
    
    result = RungeKutta(y, eta0, a, b, h);
    f0 = result.back() - y_right;
    result = RungeKutta(y, eta1, a, b, h);
    f1 = result.back() - y_right;

    while (abs(f1) > eps)
    {
        temp = eta1;
        eta1 = eta1 - (eta1 - eta0) / (f1 - f0) * f1;
        eta0 = temp;
        f0 = f1;
        result = RungeKutta(y, eta1, a, b, h);
        f1 = result.back() - y_right;
    }

    return result;
}

std::vector<double> Error(std::vector<double> result, std::vector<double> exact_result)
{
    std::vector<double> error;
    for(int i = 0; i < result.size(); i++)
    {
        error.push_back(abs(result[i] - exact_result[i]));
    }

    return error;
}

void RungeRombergRichardson(double a, double b, double h, double y, double y_right, double eps, std::vector<double> exact_result)
{
    std::vector<double> step = ShootingMethod(a, b, h, y, y_right, eps);
    std::vector<double> halfstep = ShootingMethod(a, b, h / 2, y, y_right, eps);
    std::vector<double> result(step.size(), 0);
    for(int i = 0; i < step.size(); i++)
    {
        result[i] = step[i] + (step[i] - halfstep[i * 2]) / (0.5 * 0.5 - 1);
    }
    printf("\nShooting method with h=%.1f using Runge-Romberg-Richardson method for accuracy: ", h);
    NiceVectorPrint(result);
    printf(" Error:");
    NiceVectorPrint(Error(result, exact_result));
}

int main()
{
    double a = 1, b = 2, h = 0.1;
    double y = ExactSolution(a), eps = 0.0001;
    double x_right = 2;
    double y_right = ExactSolution(x_right);
    
    std::cout << "\nResults of Exact Solution:";
    std::vector<double> result = MakeSequence(a, b, h);
    NiceVectorPrint(result);
    std::vector<double> exact_result = ExactResult(a, b, h);      
    NiceVectorPrint(exact_result);

    result = ShootingMethod(a, b, h, y, y_right, eps);
    std::cout << "\nResults of the Shooting method:";
    NiceVectorPrint(result);
    RungeRombergRichardson(a, b, h, y, y_right, eps, exact_result);
}