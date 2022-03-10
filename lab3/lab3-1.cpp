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
        printf("%.4f ", vector[i]);
    }
    printf("\n");
}

double Func(double x)
{
    double result;

    result = M_PI / 2 - atan(x) + x;
    //result = sin(M_PI * x / 6);

    return result;
}

double OmegaFunc(std::vector<double> X, int pos)
{
    double result = 1;

    for(int i = 0; i < X.size(); i++)
    {
        if(i != pos)
        {
            result *= (X[pos] - X[i]);
            //printf("\nX[%d] - X[%d] = %.4f - %.4f = %.4f", pos, i, X[pos], X[i], X[pos] - X[i]);
        }
    }

    return result;
}

std::vector<double> LagrangeInterpolation(std::vector<double> x)
{
    std::vector<double> result(4, 0);

    for(int i = 0; i < x.size(); i++)
    {
        //printf("\n F(x%d) = %.4f", i, Func(X[i]));
        //printf("\n W(x%d) = %.4f", i, OmegaFunc(X, i));
        result[i] = Func(x[i]) / OmegaFunc(x, i);
    }

    return result;
}

double LagrangePolynomial(std::vector<double> coef, std::vector<double> x, double point)
{
    double result = 0, temp = 1;
    for(int i = 0; i < coef.size(); i++)
    {
        for(int j = 0; j < x.size(); j++)
        {
            if(j != i)
            {
                temp *= (point - x[j]);
            }
        }
        result += coef[i] * temp;
        temp = 1;
    }

    return result;
}

std::vector<double> NewtonInterpolation(std::vector<double> x)
{
    int n = x.size(), k = 0;
    std::vector<double> result(n, 0);
    std::vector<double> temp(n, 0);
    std::vector<double> new_temp(n, 0);

    for(int i = k; i < n; i++)
    {
        temp[i] = Func(x[i]);
    }
    result[k] = temp[k];
    k++;

    for(int j = 0; j < n - 1; j++)
    {
        for(int i = k; i < n; i++)
        {
            new_temp[i] = (temp[i - 1] - temp[i]) / (x[0] - x[k]);
            //printf("\n(temp[%d] - temp[%d]) / (x[%d] - x[%d]) = (%.4f - %.4f) / (%.4f - %.4f) = %.4f", i - 1, i, 0, k, temp[i-1], temp[i], x[0], x[k], new_temp[i]);
        }
        temp = new_temp;
        result[k] = temp[k];
        k++;
    }

    return result;
}

double NewtonPolynomial(std::vector<double> coef, std::vector<double> x, double point)
{
    double result = coef[0], temp = 1;

    for(int i = 1; i < coef.size(); i++)
    {
        for(int j = 0; j < i; j++)
        {
            temp *= (point - x[j]);
        }
        result += coef[i] * temp;
        temp = 1;
    }

    return result;
}

int main()
{
    std::vector<double> Xa, Xb;
    double point;

    // Read points from file
    //-------------------------------------
    std::ifstream inFile("lab3-1_points.txt");

    if (inFile.good())
    {
        std::string str;
        getline(inFile, str);
        std::istringstream string_stream(str);
        double num;
        while(string_stream >> num)
        {
            Xa.push_back(num);
        }
        getline(inFile, str);
        string_stream = std::istringstream(str);
        while(string_stream >> num)
        {
            Xb.push_back(num);
        }
        inFile >> point;
    }
    
    std::vector<double> coef = LagrangeInterpolation(Xa);
    std::cout << "\nLagrange Interpolating Polynomial: ";
    NiceVectorPrint(coef);
    printf("The value of Lagrange Integrating Polynomial in point X*: %.4f\n", LagrangePolynomial(coef, Xa, point));

    std::cout << "\nNewton Interpolating Polynomial: ";
    coef = NewtonInterpolation(Xa);
    NiceVectorPrint(coef);
    printf("The value of Newton Integrating Polynomial in point X*: %.4f\n", NewtonPolynomial(coef, Xa, point));
}

// -3 -1 1 3
// -3 0 1 3
// -0.5

//y=0.0038\left(x-3\right)\left(x+1\right)\left(x-1\right)+0.0848\left(x+3\right)\left(x-3\right)\left(x-1\right)-0.1116\left(x+3\right)\left(x-3\right)\left(x+1\right)+0.0692\left(x+3\right)\left(x+1\right)\left(x-1\right)
//y=-0.1802+0.7682\left(x+3\right)-0.1384\left(x+3\right)\left(x+1\right)+0.0461\left(x+3\right)\left(x+1\right)\left(x-1\right)