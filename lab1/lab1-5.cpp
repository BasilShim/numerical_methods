#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h> 

struct QR
{
    std::vector<std::vector<double>> q_matrix;
    std::vector<std::vector<double>> r_matrix;
};

template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
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
            printf("%.4f ", *col);
        }
        printf("\n");
    }
    printf("\n");
}

void NiceVectorPrint(std::vector<double> vector)
{
    std::vector<double>::iterator col;
    printf("\n");
    for (col = vector.begin(); col != vector.end(); col++)
    {
            printf("%.4f ", *col);
    }
    printf("\n");
}

void PrintDecomposition(struct QR qr_decomp)
{
    std::cout << "\nQ:";
    NiceMatrixPrint(qr_decomp.q_matrix);
    std::cout << "\nR:";
    NiceMatrixPrint(qr_decomp.r_matrix);
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

std::vector<std::vector<double>> MultiplyMatrices(std::vector<std::vector<double>> matrix1, std::vector<std::vector<double>> matrix2)
{
    std::vector<std::vector<double>> result(matrix1[0].size(), std::vector<double>(matrix2[0].size(), 0));
    for(int i = 0; i < matrix1[0].size(); i++)
    {
        for (int j = 0; j < matrix2[0].size(); j++)
        {
            for(int k = 0; k < matrix2.size(); k++)
            {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    return result;
}

std::vector<std::vector<double>> ColumnByRow(std::vector<double> vector)
{
    std::vector<std::vector<double>> result(vector.size(), std::vector<double>(vector.size(), 0));
    for(int i = 0; i < vector.size(); i++)
    {
        for (int j = 0; j < vector.size(); j++)
        {
            result[i][j] = vector[i] * vector[j];
        }
    }

    return result;
}

double RowByColumn(std::vector<double> vector)
{
    double result = 0;

    for(int i = 0; i < vector.size(); i++)
    {
        result += vector[i] * vector[i];
    }

    return result;
}

std::vector<std::vector<double>> MatrixDifference(std::vector<std::vector<double>> matrix1, std::vector<std::vector<double>> matrix2)
{
    std::vector<std::vector<double>> result(matrix1[0].size(), std::vector<double>(matrix1[0].size(), 0));

    for(int i = 0; i < matrix1[0].size(); i++)
    {
        for(int j = 0; j < matrix1[0].size(); j++)
        {
            result[i][j] = matrix1[i][j] - matrix2[i][j];
        }
    }

    return result;
}

std::vector<std::vector<double>> DivideByNumber(std::vector<std::vector<double>> matrix, double num)
{
    std::vector<std::vector<double>> result = matrix;

    for(int i = 0; i < result[0].size(); i++)
    {
        for(int j = 0; j < result[0].size(); j++)
        {
            result[i][j] /= num;
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

std::vector<std::vector<double>> TransposeMatrix(std::vector<std::vector<double>> matrix)
{
    std::vector<std::vector<double>> result(matrix[0].size(), std::vector<double>(matrix[0].size(), 0));
    for(int i = 0; i < matrix[0].size(); i++)
    {
        for (int j = 0; j < matrix[0].size(); j++)
        {
            result[j][i] = matrix[i][j];
        }
    }

    return result;
}

double VectorNorm(std::vector<double> vector)
{
    double result = 0;

    for(int i = 0; i < vector.size(); i++)
    {
        result += vector[i] * vector[i];
    }

    return sqrt(result);
}

std::vector<std::vector<double>> HouseholderTransformation(std::vector<std::vector<double>> matrix, int col)
{
    std::vector<std::vector<double>> result = TransposeMatrix(matrix);
    std::vector<double> vector = result[col];
    std::vector<std::vector<double>> temp(matrix[0].size(), std::vector<double>(matrix[0].size(), 0));

    if  (col > 0)
    {
        vector[col - 1] = 0;
    }
    vector[col] += sign(vector[col]) * VectorNorm(vector);

    temp = DivideByNumber(DivideByNumber(ColumnByRow(vector), RowByColumn(vector)), 0.5);
    result = MatrixDifference(IdentityMatrix(matrix[0].size()), temp);

    return result;
}

struct QR QR_Decomposition(std::vector<std::vector<double>> matrix)
{
    QR result;
    result.q_matrix = IdentityMatrix(matrix[0].size());
    result.r_matrix = matrix;
    std::vector<std::vector<double>> householder(matrix[0].size(), std::vector<double>(matrix[0].size(), 0));

    for(int i = 0; i < matrix[0].size() - 1; i++)
    {
        householder = HouseholderTransformation(result.r_matrix, i);
        result.r_matrix = MultiplyMatrices(householder, result.r_matrix);
        result.q_matrix = MultiplyMatrices(result.q_matrix, householder);
    }

    return result;
}

bool CheckConvergence(std::vector<std::vector<double>> matrix, double eps)
{
    bool result = true;
    for(int i = 0; i < matrix[0].size(); i++)
    {
        for(int j = 0; j < matrix[0].size(); j++)
        {
            if(i != j)
            {
                if(abs(matrix[i][j]) > eps)
                {
                    result = false;
                }
            }
        }
    }

    return result;
}

std::vector<double> Eigenvalues(std::vector<std::vector<double>> matrix)
{
    int n = matrix[0].size();
    std::vector<double> result(n, 0);
    std::vector<double> vec;
    double eps;

    // Read epsilon
    //-------------
    std::cout << "Please, enter epsilon value\n";
    std::cin >> eps;
    double convergence;

    do
    {
        QR qr_decomp = QR_Decomposition(matrix);
        matrix = MultiplyMatrices(qr_decomp.r_matrix, qr_decomp.q_matrix);
        qr_decomp = QR_Decomposition(matrix);
        vec = std::vector<double>();
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if(i > j)
                {
                    vec.push_back(matrix[i][j]);
                }
            }
        }
        convergence = VectorNorm(vec);

    } while(convergence > eps);

    for(int i = 0; i < matrix[0].size(); i++)
    {
        result[i] = matrix[i][i];
    }
    std::cout << "\nEigenvalues:";

    return result;
}

int main()
{
    int n = 3;
    std::vector<std::vector<double>> matrix;

    // Read symmetric matrix from file
    //-------------------------------------
    std::ifstream inFile("lab1-5_matrix.txt");

    if (inFile.good())
    {
        std::string str;
        while(getline(inFile, str)) 
        {
            std::istringstream string_stream(str);
            double num;
            std::vector<double> temp_vector;
            while(string_stream >> num)
            {
                temp_vector.push_back(num);
            }
            matrix.push_back(temp_vector);
        }
    }

    NiceMatrixPrint(matrix);

    // QR-Decomposition
    //-----------------
    QR qr_decomop = QR_Decomposition(matrix);
    PrintDecomposition(qr_decomop);

    // Let's check the QR-Decomposition by performing QR=A
    //----------------------------------------------------
    std::cout << "QR=A";
    NiceMatrixPrint(MultiplyMatrices(qr_decomop.q_matrix, qr_decomop.r_matrix));

    // Finding eigenvalues
    //--------------------
    NiceVectorPrint(Eigenvalues(matrix));
}

// -1 8 5
// 8 -4 4
// 2 9 -2