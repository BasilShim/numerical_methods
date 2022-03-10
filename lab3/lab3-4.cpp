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

double FirstDerivative(std::vector<double> x, std::vector<double> y, double point)
{
    int n = x.size();
    double result, lh_der, rh_der;
    for(int i = 0; i < n; i++)
    {
        if(x[i] == point && i > 0)
        {
            lh_der = (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
            printf("\nLeft-hand Derivative: %.4f", lh_der);
            rh_der = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
            printf("\nRight-hand Derivative: %.4f", rh_der);
            result = lh_der + (rh_der - lh_der) / (x[i + 1] - x[i - 1]) * (2 * point - x[i - 1] - x[i]);
            break;
        }
    }

    return result;
}

double SecondDericative(std::vector<double> x, std::vector<double> y, double point)
{
    int n = x.size();
    double result;

    for(int i = 0; i < n; i++)
    {
        if(x[i] == point && i > 0)
        {
            result = 2 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1])) / (x[i + 1] - x[i - 1]);
            //printf("\n2 * ((%.4f - %.4f) / (%.4f - %.4f) - (%.4f - %.4f) / (%.4f - %.4f)) / (%.4f - %.4f)", 
            //y[i + 1], y[i], x[i + 1], x[i], y[i], y[i - 1], x[i] , x[i - 1], x[i + 1], x[i - 1]);
            break;
        }
    }

    return result;
}

int main()
{
    std::vector<double> x, y;
    double point;

    // Read points from file
    //-------------------------------------
    std::ifstream inFile("lab3-4_points.txt");

    if (inFile.good())
    {
        std::string str;
        getline(inFile, str);
        std::istringstream string_stream(str);
        double num;
        while(string_stream >> num)
        {
            x.push_back(num);
        }
        getline(inFile, str);
        string_stream = std::istringstream(str);
        while(string_stream >> num)
        {
            y.push_back(num);
        }
        inFile >> point;
    }

    NiceVectorPrint(x);
    NiceVectorPrint(y);
    std::cout << "\n" << point << "\n"; 
    double first_derivative = FirstDerivative(x, y, point);
    printf("\nSecond order of accuracy First Derivative: %.4f", first_derivative);

    double second_derivative = SecondDericative(x, y, point);
    printf("\nSecond Derivative: %.4f", second_derivative);
}

// 1.0 1.2 1.4 1.6 1.8   
// 1.0 0.69444 0.5102 0.39062 0.30864
// 1.4