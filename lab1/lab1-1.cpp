// Шимко Василий
// M8О-305Б-18
// Вариант 22

#include <iostream>
#include <vector>
#include <fstream>
#include <string>

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

std::vector<double> ForwardSubstitution(std::vector<std::vector<double>> l_matrix, std::vector<double> b_vector, int n)
{
    std::vector<double> z_vector(n, 0);

    for(int i = 0; i < n; i++)
    {
        z_vector[i] = b_vector[i];
        for(int j = 0; j < i; j++)
        {
            z_vector[i] -= l_matrix[i][j] * z_vector[j];
        }
    }
    return z_vector;
}

std::vector<double> BackwardSubstitution(std::vector<std::vector<double>> u_matrix, std::vector<double> z_vector, int n)
{
    std::vector<double> x_vector(n, 0);

    for(int i = n - 1; i > -1; i--)
    {
        x_vector[i] = z_vector[i];
        for(int j = i + 1; j < n; j++)
        {
            x_vector[i] -= u_matrix[i][j] * x_vector[j];
        }
        x_vector[i] /= u_matrix[i][i];
    }

    return x_vector;
}

double Determinant(std::vector<std::vector<double>> l_matrix, std::vector<std::vector<double>> u_matrix, int n)
{
    double result = l_matrix[0][0];
    for(int i = 0; i < n; i++)
    {
        result *= l_matrix[i][i];
    }

    for(int i = 0; i < n; i++)
    {
        result *= u_matrix[i][i];
    }
    return result;
}

std::vector<std::vector<double>> InverseMatrix(std::vector<std::vector<double>> l_matrix, std::vector<std::vector<double>> u_matrix, int n)
{
    std::vector<std::vector<double>> result(n, std::vector<double>(n, 0));
    std::vector<double> temp_b(n, 0);
    std::vector<double> temp_x(n, 0);
    std::vector<double> temp_z(n, 0);

    for(int k = 0; k < n; k++)
    {
        std::fill(temp_b.begin(), temp_b.end(), 0);
        temp_b[k] = 1.0f;
        temp_z = ForwardSubstitution(l_matrix, temp_b, n);
        temp_x = BackwardSubstitution(u_matrix, temp_z, n);
        for(int i = 0; i < n; i++)
        {
            result[i][k] = temp_x[i];
        }
    }

    return result;
}

int main()
{
    int n = 4;
    std::vector<std::vector<double>> systemCoef;

    // Read system's coefficients from file
    //-------------------------------------
    std::ifstream inFile("lab1-1_coef.txt");

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

    std::vector<std::vector<double>> l_matrix(n, std::vector<double>(n, 0));
    std::vector<std::vector<double>> u_matrix(n, std::vector<double>(n, 0));

    l_matrix = IdentityMatrix(n);

    // LU-decomposition
    //-----------------
    std::vector<std::vector<double>> tempMatrix = systemCoef;
    for(int k = 0; k < n; k++)
    {
        u_matrix[k][k] = tempMatrix[k][k];
        for(int i = k + 1; i < n; i++)
        {
            l_matrix[i][k] = tempMatrix[i][k] / u_matrix[k][k];
            u_matrix[k][i] = tempMatrix[k][i];
        }
        for(int i = k + 1; i < n; i++)
        {
            for (int j = k + 1; j < n; j++)
            {
                tempMatrix[i][j] = tempMatrix[i][j] - l_matrix[i][k] * u_matrix[k][j];
            }
        }
    }

    // Print L and U matrices
    //-----------------------
    std::cout << "\nResult L and U matrices are:\n";
    NiceMatrixPrint(u_matrix, n);
    NiceMatrixPrint(l_matrix, n);

    // Check LU-decomposition by performing LU=A
    //-----------------------
    std::cout << "\nCheck L and U matrices:\n";
    NiceMatrixPrint(MultiplyMatrices(l_matrix, u_matrix, n), n);

    std::vector<double> b_vector(n, 0);

    std::cout << "\nPlease, enter individually each element of the vector b:\n";
    for(int i = 0; i < n; i++)
    {
        std::cin >> b_vector[i];
    }

    std::cout << "\nThe vector you've entered is:\n";
    NiceVectorPrint(b_vector, n);

    // Finding z from Lz=b
    //--------------------
    std::vector<double> z_vector = ForwardSubstitution(l_matrix, b_vector, n);

    // Finding x from Ux=b
    //--------------------
    std::vector<double> x_vector = BackwardSubstitution(u_matrix, z_vector, n);

    std::cout << "\nThe result vector is:\n";
    NiceVectorPrint(x_vector, n);

    // Check result
    //-------------
    std::cout << "\nLet's put the roots into the system to check the results:\n";
    NiceVectorPrint(CheckRoots(systemCoef, x_vector, n), n);

    // Determinant of matrix
    //----------------------
    std::cout << "\nDeterminant of matrix:\n" << Determinant(l_matrix, u_matrix, n) << "\n";

    // Inverse matrix
    //---------------
    std::cout << "\nInverse matrix:\n";
    NiceMatrixPrint(InverseMatrix(l_matrix, u_matrix, n), n);
}


// 120
// 31
// 6
// 25