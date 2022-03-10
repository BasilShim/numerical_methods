#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>

struct Condition
{
    double x;
    double y;
};

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

double Function(double x, double y)
{   
    double result;

    result = - y * (x - 2) / (x * x);
    //result = (x + y) * (x + y);

    return result;
}

double ExactSolution(double x)
{
    double result;

    result = sin(x - 1) + cos(x - 1) / x;
    //result = tan(x) - x;

    return result;
}

std::vector<std::vector<double>> ExplicitEuler(Condition con, double a, double b, double h)
{
    double eps;
    int i = 1;
    std::vector<std::vector<double>> result(2, std::vector<double>());
    result[0].push_back(con.x);
    result[1].push_back(con.y);

    for(;;)
    {
        result[1].push_back(result[1][i - 1] + h * Function(result[0][i - 1], result[1][i - 1]));
        result[0].push_back(result[0][i - 1] + h);
        eps = abs(ExactSolution(result[0][i]) - result[1][i]);
        if(result[0][i] >= b)
        {
            break;
        }
        i++;
    }
    
    return result;
}

std::vector<std::vector<double>> RungeKutta(Condition con, double a, double b, double h)
{
    double eps, delta, k1, k2, k3, k4;
    int i = 1;
    std::vector<std::vector<double>> result(2, std::vector<double>());
    result[0].push_back(con.x);
    result[1].push_back(con.y);

    for(;;)
    {
        k1 = h * Function(result[0][i - 1], result[1][i - 1]);
        k2 = h * Function(result[0][i - 1] + h / 2, result[1][i - 1] + k1 / 2);
        k3 = h * Function(result[0][i - 1] + h / 2, result[1][i - 1] + k2 / 2);
        k4 = h * Function(result[0][i - 1] + h, result[1][i - 1] + k3);
        delta = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        result[1].push_back(result[1][i - 1] + delta);
        result[0].push_back(result[0][i - 1] + h);
        if(result[0][i] >= b)
        {
            break;
        }
        i++;
    }

    return result;
}

std::vector<std::vector<double>> Adams(Condition con, double a, double b, double h, std::vector<std::vector<double>> result)
{
    int i = result[0].size() - 2;
    for(;;)
    {
        result[1][i] = result[1][i - 1] + (55 * Function(result[0][i - 1], result[1][i - 1]) - 59 * Function(result[0][i - 2], result[1][i - 2]) 
        + 37 * Function(result[0][i - 3], result[1][i - 3]) - 9 * Function(result[0][i - 4], result[1][i - 4])) * h / 24;
        if(result[0][i] >= b)
        {
            break;
        }
        i++;
    }

    return result;
}

std::vector<std::vector<double>> ExactResult(double a, double b, double h)
{
    std::vector<std::vector<double>> result(2, std::vector<double>());

    for(double i = a; i <= b + h; i+=h)
    {
        result[0].push_back(i);
        result[1].push_back(ExactSolution(i));
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

void RungeRombergRichardson(Condition con, double a, double b, double h, std::vector<double> exact_result)
{
    std::vector<double> step = ExplicitEuler(con, a, b, h)[1];
    std::vector<double> halfstep = ExplicitEuler(con, a, b, h / 2)[1];
    std::vector<double> result(step.size(), 0);
    for(int i = 0; i < step.size(); i++)
    {
        result[i] = step[i] + (step[i] - halfstep[i]) / (0.5 * 0.5 - 1);
    }
    printf("\nExplicit Euler's method with h=%.1f using Runge-Romberg-Richardson method for accuracy: ", h);
    NiceVectorPrint(result);
    printf(" Error:");
    NiceVectorPrint(Error(result, exact_result));

    step = RungeKutta(con, a, b, h)[1];
    halfstep = RungeKutta(con, a, b, h / 2)[1];
    for(int i = 0; i < step.size(); i++)
    {
        result[i] = step[i] + (step[i] - halfstep[i]) / (0.5 * 0.5 - 1);
    }
    printf("\nRunge-Kutta method with h=%.1f using Runge-Romberg-Richardson method for accuracy: ", h);
    NiceVectorPrint(result);
    printf(" Error:");
    NiceVectorPrint(Error(result, exact_result));

    step = Adams(con, a, b, h, RungeKutta(con, a, b, h))[1];
    halfstep = Adams(con, a, b, h / 2, RungeKutta(con, a, b, h / 2))[1];
    for(int i = 0; i < step.size(); i++)
    {
        result[i] = step[i] + (step[i] - halfstep[i]) / (0.5 * 0.5 - 1);
    }
    printf("\nAdam's method with h=%.1f using Runge-Romberg-Richardson method for accuracy: ", h);
    NiceVectorPrint(result);
    printf(" Error:");
    NiceVectorPrint(Error(result, exact_result));
}

int main()
{
    Condition con;
    con.x = 1;
    con.y = 1;

    double a = 1, b = 2, h = 0.1;
    //double a = 0, b = 0.5, h = 0.1;

    std::vector<std::vector<double>> result = ExactResult(a, b, h);
    std::vector<double> exact_result = result[1];
    std::cout << "\nResults of Exact Solution:";
    NiceMatrixPrint(result);
    result = ExplicitEuler(con, a, b, h);
    std::cout << "\nResults of Explicit Euler's method:";
    NiceMatrixPrint(result);
    result = RungeKutta(con, a, b, h);
    std::cout << "\nResults of 4th Order Runge-Kutta method:";
    NiceMatrixPrint(result);
    result = Adams(con, a, b, h, result);
    std::cout << "\nResults of 4th Order Adams' method:";
    NiceMatrixPrint(result);
    RungeRombergRichardson(con, a, b, h, exact_result);
}