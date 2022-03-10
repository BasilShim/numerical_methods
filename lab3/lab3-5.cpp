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
        printf("%.8f ", vector[i]);
    }
    printf("\n");
}

std::vector<double> MakeSequence(double x0, double xk, double step)
{
    std::vector<double> result;

    for(double x = x0; x <= xk; x += step)
    {
        result.push_back(x);
    }

    return result;
}

double Function(double x)
{
    double result;

    result = x * sqrt(2 * x + 3);
    //result = x / ((3 * x + 4) * (3 * x + 4));

    return result;
}

double RectangleMethod(std::vector<double> x, double step)
{
    double result = 0;
    int n = x.size();

    for(int i = 0; i < n - 1; i++)
    {
        result += step * Function((x[i] + x[i + 1]) / 2);
    }

    return result;
}

double TrapezoidalRule(std::vector<double> x, double step)
{
    double result = 0;
    int n = x.size();

    for(int i = 0; i < n - 1; i++)
    {
        result += step / 2 * (Function(x[i]) + Function(x[i + 1]));
    }

    return result;
}

double SimpsonsMethod(std::vector<double> x, double step)
{
    double result = 0;
    int n = x.size();
    step /= 6;

    for(int i = 0; i < n - 1; i++)
    {
        result += step * (Function(x[i]) + 4 * Function((x[i] + x[i + 1]) / 2) + Function(x[i + 1]));
    }

    return result;
}

void RungeRombergRichardson(std::vector<double> x, double step, double acc_val)
{
    std::vector<double> temp_vec = MakeSequence(x[0], x[x.size() - 1], step / 2);

    double result = RectangleMethod(x, step);
    result = result + (result - RectangleMethod(temp_vec, step / 2)) / (0.5 * 0.5 - 1);
    printf("\n\nRectangle Method with h=%.8f using Runge-Romberg-Richardson Method for accuracy: %.8f", step, result);
    printf("\n Error: %.8f", abs(result - acc_val));

    result = TrapezoidalRule(x, step);
    result = result + (result - TrapezoidalRule(temp_vec, step / 2)) / (0.5 * 0.5 - 1);
    printf("\nTrapezoidal Rule with h=%.8f using Runge-Romberg-Richardson Method for accuracy: %.8f", step, result);
    printf("\n Error: %.8f", abs(result - acc_val));

    result = SimpsonsMethod(x, step);
    result = result + (result - SimpsonsMethod(temp_vec, step / 2)) / (powf(0.5, 4) - 1);
    printf("\nSimpson's Method with h=%.8f using Runge-Romberg-Richardson Method for accuracy: %.8f", step, result);
    printf("\n Error: %.8f", abs(result - acc_val));
}

int main()
{
    double x0, xk, step1, step2, acc_val;

    // Read system's coefficients from file
    //-------------------------------------
    std::ifstream inFile("lab3-5_params.txt");

    inFile >> x0;
    inFile >> xk;
    inFile >> step1;
    inFile >> step2;
    inFile >> acc_val;

    //printf("\n %.8f %.8f %.8f %.8f", x0, xk, step1, step2);
    std::vector<double> x_seq1 = MakeSequence(x0, xk, step1);
    std::vector<double> x_seq2 = MakeSequence(x0, xk, step2);

    //Rectangle Method
    //----------------
    double integral = RectangleMethod(x_seq1, step1);
    printf("\nRectangle Method with h=%.8f: %.8f", step1, integral); 
    integral = RectangleMethod(x_seq2, step2);
    printf("\nRectangle Method with h=%.8f: %.8f", step2, integral);

    //Trapezoidal Rule
    //----------------
    integral = TrapezoidalRule(x_seq1, step1);
    printf("\n\nTrapezoidal Rule with h=%.8f: %.8f", step1, integral);
    integral = TrapezoidalRule(x_seq2, step2);
    printf("\nTrapezoidal Rule with h=%.8f: %.8f", step2, integral);

    //Simpson's Method
    //----------------
    integral = SimpsonsMethod(x_seq1, step1);
    printf("\n\nSimpson's Method with h=%.8f: %.8f", step1, integral);
    integral = SimpsonsMethod(x_seq2, step2);
    printf("\nSimpson's Method with h=%.8f: %.8f", step2, integral);

    //Each method with increased accuracy using Runge-Romberg-Richardson Method
    //-------------------------------------------------------------------------
    RungeRombergRichardson(x_seq1, step1, acc_val);
    RungeRombergRichardson(x_seq2, step2, acc_val);
}

// -1 1 0.5 0.25