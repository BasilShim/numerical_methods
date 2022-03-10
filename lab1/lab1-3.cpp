// Шимко Василий
// M8О-305Б-18
// Вариант 22

#include <iostream>
#include <vector>
#include <fstream>
#include <string>

void NiceMatrixPrint(std::vector<std::vector<double>> matrix, int n)
{
    printf("\n");
    for(int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.4f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void NiceVectorPrint(std::vector<double> vector, int n)
{
    printf("\n");
    for(int i = 0; i < n; i++)
    {
        printf("%.4f ", vector[i]);
    }
    printf("\n");
}

std::vector<double> CheckRoots(std::vector<std::vector<double>> matrix, std::vector<double> vector, int n)
{
    std::vector<double> result(n, 0);
    double temp = 0;

    for(int i = 0; i < n; i++)
    {
        result[i] = 0;
        for(int j = 0; j < n; j++)
        {
            result[i] += matrix[i][j] * vector[j];
        }
    }
    return result;
}

std::vector<double> JacobiMethod(std::vector<std::vector<double>> matrix, std::vector<double> b_vector, std::vector<double> x_vector, int n, double eps)
{
    std::vector<double> prev_x = x_vector;
    double convergence;
    int iteration = 1;
    double sum = 0;

    std::cout << "\nJacobi's method\n";
    do
    {
        std::cout << "step: " << iteration << "\n";
        prev_x = x_vector;
        for(int i = 0; i < n; i++)
        {
            sum = 0;
            for(int j = 0; j < n; j++)
            {
                if(j != i)
                {
                    sum += matrix[i][j] * x_vector[j]; 
                }
            }
            x_vector[i] = (b_vector[i] - sum) / matrix[i][i];
        }

        convergence = std::abs(x_vector[0] - prev_x[0]);
        for(int i = 1; i < n; i++)
        {
            if(std::abs(x_vector[i] - prev_x[i]) > convergence)
            {
                convergence = x_vector[i] - prev_x[i];
            }
        }
        iteration++;
    } while (convergence > eps);
    
    std::cout << "\nNumber of iterations of Jacobi's Method:\n" << iteration - 1;

    return x_vector;
}

std::vector<double> GaussSeidelMethod(std::vector<std::vector<double>> matrix, std::vector<double> b_vector, std::vector<double> x_vector, int n, double eps)
{
    std::vector<double> prev_x = x_vector;
    double convergence;
    int iteration = 1;
    double sum = 0;

    std::cout << "\nGauss-Seidel's Method\n";
    do
    {
        std::cout << "step: " << iteration << "\n";
        prev_x = x_vector;
        for(int i = 0; i < n; i++)
        {
            sum = 0;
            for(int j = 0; j < i; j++)
            {
                sum += matrix[i][j] * x_vector[j];
            }
            for(int j = i + 1; j < n; j++)
            {
                sum += matrix[i][j] * prev_x[j];
            }
            x_vector[i] = (b_vector[i] - sum) / matrix[i][i];
        }

        convergence = std::abs(x_vector[0] - prev_x[0]);
        for(int i = 1; i < n; i++)
        {
            if(std::abs(x_vector[i] - prev_x[i]) > convergence)
            {
                convergence = x_vector[i] - prev_x[i];
            }
        }
        iteration++;
    } while (convergence > eps);
    
    std::cout << "\nNumber of iterations of Gauss-Seidel's Method:\n" << iteration - 1;

    return x_vector;
}

int main()
{
    int n = 4;
    double eps;
    std::vector<std::vector<double>> systemCoef;

    // Read system's coefficients from file
    //-------------------------------------
    std::ifstream inFile("lab1-3_coef.txt");

    double temp = 0;
    for(int i = 0; i < n; i++)
    {
        std::vector<double> tempVec;
        for(int j = 0; j < n; j++)
        {
            inFile >> temp;
            tempVec.push_back(temp);
        }
        systemCoef.push_back(tempVec);
    }

    NiceMatrixPrint(systemCoef, n);

    // Read epsilon
    //-------------
    std::cout << "Please, enter epsilon value\n";
    std::cin >> eps;

    // Read vector b
    //--------------

    std::vector<double> b_vector(n, 0);
    std::cout << "\nPlease enter individually each element of the vector b:\n";
    for(int i = 0; i < n; i++)
    {
        std::cin >> b_vector[i];
    }

    std::cout << "\nThe vector you've entered is:\n";
    NiceVectorPrint(b_vector, n);

    // Solve the system, using Jacobi's Method
    //----------------------------------------

    std::vector<double> x_vector(n, 1);
    x_vector = JacobiMethod(systemCoef, b_vector, x_vector, n, eps);

    std::cout << "\nThe result vector is:\n";
    NiceVectorPrint(x_vector, n);

    // Check result
    //-------------
    std::cout << "\nLet's put the roots into the system to check the results:\n";
    NiceVectorPrint(CheckRoots(systemCoef, x_vector, n), n);

    // Solve the system, using Gauss-Seidel's Method
    //----------------------------------------

    std::fill(x_vector.begin(), x_vector.end(), 1);
    x_vector = GaussSeidelMethod(systemCoef, b_vector, x_vector, n, eps);

    std::cout << "\nThe result vector is:\n";
    NiceVectorPrint(x_vector, n);

    // Check result
    //-------------
    std::cout << "\nLet's put the roots into the system to check the results:\n";
    NiceVectorPrint(CheckRoots(systemCoef, x_vector, n), n);
}

// -60
// -109
// -103
// -33