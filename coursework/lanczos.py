# Шимко Василий
# Группа: М8О-305Б-18
# Вариант 3

# Import numpy for most stuff and scipy to generate a sparse matrix
import numpy as np
import timeit
import scipy
from scipy.sparse import random
from scipy import stats

np.set_printoptions(formatter={'float': lambda x: "{0:0.8f}".format(x)}, linewidth=200)

# Class for CSR representation
class CSR:
    values = np.array([])
    row_ptr = np.array([])
    col = np.array([])

# Function for multiplication of a CSR matrix by a vector
def CSRByVector(csr, vector):
    result = np.zeros(len(vector))

    for i in range(len(result)):
        for j in range(int(csr.row_ptr[i]), int(csr.row_ptr[i + 1])):
            result[i] = result[i] + csr.values[j] * vector[int(csr.col[j])]
    
    return result

def PrintCSR(csr):
    print("Values:\n", csr.values)
    print("Row pointers:\n", csr.row_ptr)
    print("Columns:\n", csr.col)

def Lanczos(csr, v_curr, n):
    # Starting calculations
    v_matrix = np.zeros((n, n))
    v_matrix[0,:] = v_curr
    tridiagonal = np.zeros((n, n))
    v_prev = np.zeros(n)
    beta = 0

    for i in range(n - 1):
        # Calculate all the greek letters and vectors
        w_vec = CSRByVector(csr, v_curr)
        alpha = w_vec @ v_curr
        w_vec = w_vec - alpha * v_curr - beta * v_prev
        beta = np.sqrt(w_vec @ w_vec) 
        v_prev = v_curr
        v_curr = w_vec / beta 
        # Shove greek letters into the tridiagonal matrix
        tridiagonal[i][i] = alpha 
        tridiagonal[i][i + 1] = beta
        tridiagonal[i + 1][i] = beta
        # Shove V vector into V matrix
        v_matrix[i,:] = v_curr
    
    # Final calculations
    w_vec = CSRByVector(csr, v_curr)
    alpha = w_vec @ v_curr
    w_vec = w_vec - alpha * v_curr - beta * v_prev
    tridiagonal[n - 1, n - 1] = alpha
    v_matrix[n - 1] = w_vec / np.sqrt(w_vec @ w_vec)

    return tridiagonal, v_matrix


def QRAlgorithm(matrix):
    n = len(matrix[0])
    eps = 0.0001
    convergence = eps + 1
    eig_vecs = np.identity(n)
    i = 0

    while (convergence > eps):
        q_matrix, r_matrix = np.linalg.qr(matrix)
        vec = np.zeros(n)
        matrix = r_matrix @ q_matrix
        eig_vecs = eig_vecs @ q_matrix
        # Vector, containing elements below the main diagonal
        for i in range(n):
            for j in range(n):
                if i > j:
                    vec[i] = matrix[i][j]
        
        convergence = np.sqrt(vec @ vec)
    # Vector of eigenvalues, which consists of the elements on the main diagonal
    eigenvalues = np.zeros(n)
    for i in range(n):
        eigenvalues[i] = matrix[i][i]

    return eigenvalues, eig_vecs

def main(): 
    # Generate and display random sparse matrix
    n = 10
    matrix = random(n, n, density=0.2)
    matrix = matrix.toarray() # Convert to numpy array
    matrix = matrix @ matrix.T # Make symmetrical matrix
    print("\nRandom sparse matrix:\n")
    print(matrix)

    # Calculate eigenvalues using numpy.eig() on starting sparse matrix
    numpy_ev, numpy_evecs = np.linalg.eig(matrix)
    numpy_ev = np.sort(numpy_ev)[::-1]

    # Save the matrix to the txt file and delete it
    np.savetxt("rand_matrix.txt", matrix)
    del matrix

    # Read sparse matrix from txt file into a CSR structure
    count = 0
    csr_matrix = CSR()
    csr_matrix.row_ptr = np.append(csr_matrix.row_ptr, count)

    with open("rand_matrix.txt") as file:
        lines = file.readlines()
        for line in lines:
            line = line.split()
            for i in range(len(line)):
                if np.float64(line[i]) != 0:
                    count = count + 1
                    csr_matrix.values = np.append(csr_matrix.values, np.float64(line[i]))
                    csr_matrix.col = np.append(csr_matrix.col, i)
            csr_matrix.row_ptr = np.append(csr_matrix.row_ptr, count)

    print("\nSparse matrix converted to CSR representation:\n")
    PrintCSR(csr_matrix)

    # Generate starting random V vector
    v_vec   = np.random.rand(n)
    v_vec /= np.sqrt(v_vec @ v_vec)

    # Do the Lanczos algorithm
    tridiagonal, v_matrix = Lanczos(csr_matrix, v_vec, n)

    # Calculate eigenvalues of tridiagonal matrix using numpy.eig()
    tridiag_ev, tridiag_evecs = np.linalg.eig(tridiagonal)
    tridiag_ev = np.sort(tridiag_ev)[::-1]

    # Calculate eigenvalues of tridiagonal matrix using QR decomposition
    qr_eigen, qr_evecs = QRAlgorithm(tridiagonal)
    qr_eigen = np.sort(qr_eigen)[::-1]

    # Print results
    print("Tridiagonal matrix:\n")
    print(tridiagonal)
    print("\nEigenvalues:\n")
    print("Using numpy.eig() on the sparse matrix:\n", numpy_ev)
    print("Using numpy.eig() on the tridiagonal matrix:\n", tridiag_ev)
    print("Using QR decomposition on the tridiagonal matrix:\n", qr_eigen)

main()