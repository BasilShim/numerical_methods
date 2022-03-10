#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <math.h> 

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

std::vector<std::vector<double>> IdentityMatrix(int n)
{
    std::vector<std::vector<double>> identityMatrix(n, std::vector<double>(n, 0));
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if(i == j)
            {
                identityMatrix[i][j] = 1;
            }
        }
    }
    return identityMatrix;
}

std::vector<std::vector<double>> MultiplyMatrices(std::vector<std::vector<double>> matrix1, std::vector<std::vector<double>> matrix2, int n)
{
    std::vector<std::vector<double>> result(n, std::vector<double>(n, 0));

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            for(int k = 0; k < n; k++)
            {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    return result;
}

std::vector<double> MaxAbsoluteValue(std::vector<std::vector<double>> matrix, int n)
{
    double max = 0;
    std::vector<double> result(2, 0);
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if(i != j && max < std::abs(matrix[i][j]))
            {
                max = std::abs(matrix[i][j]);
                result[0] = i;
                result[1] = j;
            }
        }
    }

    return result;
}

std::vector<std::vector<double>> TransposeMatrix(std::vector<std::vector<double>> matrix, int n)
{
    std::vector<std::vector<double>> result = matrix;

    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            result[j][i] = matrix[i][j];
        }
    }

    return result;
}

std::vector<std::vector<double>> JacobiRotation(std::vector<std::vector<double>> matrix, int n, double eps)
{
    std::vector<double> eigenvalues(n, 0);
    std::vector<std::vector<double>> eigenvectors = IdentityMatrix(n);
    std::vector<std::vector<double>> rotation_matrix(n, std::vector<double>(n, 0));
    std::vector<std::vector<double>> a_matrix;
    std::vector<double> max_position;
    double convergence = eps + 1, phi = 0;

    while(convergence > eps)
    {
        convergence = 0;
        a_matrix = matrix;
        max_position = MaxAbsoluteValue(matrix, n);

        if(a_matrix[max_position[0]][max_position[0]] != a_matrix[max_position[1]][max_position[1]])
        {
            phi = std::atan(2 * a_matrix[max_position[0]][max_position[1]] / (a_matrix[max_position[0]][max_position[0]] - a_matrix[max_position[1]][max_position[1]])) / 2;
        }
        else
        {
            phi = M_PI / 4;
        }

        for (auto &i : rotation_matrix)
        {
            std::fill(i.begin(), i.end(), 0);
        }

        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if(i == max_position[0] && j == max_position[1])
                {
                    rotation_matrix[i][j] = -std::sin(phi);
                    rotation_matrix[i][i] = std::cos(phi);
                    rotation_matrix[j][j] = std::cos(phi);
                    rotation_matrix[j][i] = std::sin(phi);
                }
                else if(i == j && i != max_position[0] && j != max_position[1])
                {
                    rotation_matrix[i][j] = 1;
                }
            }
        }
        matrix = MultiplyMatrices(MultiplyMatrices(TransposeMatrix(rotation_matrix, n), a_matrix, n), rotation_matrix, n);
        eigenvectors = MultiplyMatrices(eigenvectors, rotation_matrix, n);
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if(i < j)
                {
                    convergence += a_matrix[i][j] * a_matrix[i][j];
                }
            }
        }
        convergence = std::sqrt(convergence);
    }
    for(int i = 0; i < n; i++)
    {
        printf("Lambda %i = %.4f\n", i + 1, matrix[i][i]);
    }
    return eigenvectors;
}

int main()
{
    int n = 3;
    double eps;
    std::vector<std::vector<double>> matrix;

    // Read symmetric matrix from file
    //-------------------------------------
    std::ifstream inFile("lab1-4_matrix.txt");

    double temp = 0;
    for(int i = 0; i < n; i++)
    {
        std::vector<double> tempVec;
        for(int j = 0; j < n; j++)
        {
            inFile >> temp;
            tempVec.push_back(temp);
        }
        matrix.push_back(tempVec);
    }

    NiceMatrixPrint(matrix, n);

    // Read epsilon
    //-------------
    std::cout << "Please, enter epsilon value\n";
    std::cin >> eps;

    // Find eigenvalues and eigenvectors, using Jacobi's rotation algorithm
    //---------------------------------------------------------------------
    NiceMatrixPrint(JacobiRotation(matrix, n, eps), n);
}

