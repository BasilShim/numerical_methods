import numpy as np
import matplotlib.pyplot as plt
import sys
np.set_printoptions(threshold=sys.maxsize)

np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)}, linewidth=200)

def func(x):
    result = x + np.sin(np.pi * x) 
    return result

def ExactSolution(x, t, a):
    result = np.zeros(len(x))
    for i in range(len(x)):
        temp = x[i] + np.exp(-np.pi * np.pi * a * t[i]) * np.sin(np.pi * x[i])
        result[i] = temp

    return result

def TridiagSolve(a_vec, b_vec, c_vec, d_vec):

    n = len(d_vec)
    roots = np.zeros(len(d_vec))
    p = np.zeros(n)
    q = np.zeros(n)

    p[0] = -c_vec[0] / b_vec[0]
    q[0] = d_vec[0] / b_vec[0]

    for i in range(1, n - 1):
        temp = (b_vec[i] + a_vec[i - 1] * p[i - 1])
        p[i] = - c_vec[i] / temp
        q[i] = (d_vec[i] - a_vec[i - 1] * q[i - 1]) / temp

    q[-1] = (d_vec[-1] - a_vec[-1] * q[-2]) / (b_vec[-1] + a_vec[-1] * p[-2])

    roots[-1] = q[-1]
    for i in range(len(q) - 2, -1, -1):
        roots[i] = p[i] * roots[i + 1] + q[i]

    return roots

def TridiagMatrix(a_vec, b_vec, c_vec):
    matrix = np.zeros((len(b_vec), len(b_vec)))

    matrix[0][0] = b_vec[0]

    for i in range (1, len(matrix)):
        matrix[i][i] = b_vec[i]
        matrix[i][i - 1] = a_vec[i - 1]
        matrix[i - 1][i] = c_vec[i - 1]

    print(matrix)

    return matrix

# Explicit Finite Difference Method

def ExFinDiff(x, t, sigma):
    prev_layer = np.ones(len(x))
    grid = np.zeros((len(t), len(x)))

    # First Layer
    for i in range(0, len(x) - 1):
        prev_layer[i] = func(x[i])
    grid[0] = prev_layer

    # Other Layers

    temp_layer = np.zeros(len(prev_layer))
    temp_layer[-1] = 1
    for i in range(1, len(t)):
        for j in range(1, len(x) - 1):
            temp_layer[j] = (1 - 2 * sigma) * prev_layer[j] + sigma * (prev_layer[j - 1] + prev_layer[j + 1])
        grid[i] = temp_layer
        prev_layer = temp_layer
        
    print(grid.shape)
    np.savetxt("exp_grid.txt", grid, fmt='%.5f')
    result = grid[-1]

    return result

# Implicit Finite Difference Method

def ImFinDiff(x, t, sigma):
    d_vec = np.ones(len(x))
    grid = np.zeros((len(t), len(x)))
    grid[:,-1] = 1

    
    d_vec = func(x)[1:-1]
    d_vec[-1] = d_vec[-1] + sigma
    grid[0, 1:-1] = d_vec

    a_vec = np.zeros(len(d_vec) - 1)
    c_vec = np.zeros(len(d_vec) - 1)
    b_vec = np.zeros(len(d_vec))
    a_vec.fill(-sigma)
    c_vec.fill(-sigma)
    b_vec.fill(1 + 2 * sigma)

    for i in range (1, len(t)):
        d_vec = TridiagSolve(a_vec, b_vec, c_vec, d_vec)
        if (i < len(t) - 1):
            d_vec[-1] = d_vec[-1] + sigma
        grid[i, 1:-1] = d_vec

    print(grid.shape)
    np.savetxt("imp_grid.txt", grid, fmt='%.5f')
    result = grid[-1]

    return result

# Crank-Nicolson Method

def CrankNicolson(x, t, sigma):
    sigma = sigma / 2
    d_vec = np.ones(len(x))
    grid = np.zeros((len(t), len(x)))
    grid[:,-1] = 1
    
    d_vec = func(x)[1:-1]
    d_vec[-1] = d_vec[-1]
    grid[0, 1:-1] = d_vec

    a_vec = np.zeros(len(d_vec) - 1)
    c_vec = np.zeros(len(d_vec) - 1)
    b_vec = np.zeros(len(d_vec))
    a_vec.fill(-sigma)
    c_vec.fill(-sigma)
    b_vec.fill(2 + 2 * sigma)

    temp = np.zeros(len(d_vec))
    for i in range(1, len(t)):
        temp[0] = (2 - 2 * sigma) * d_vec[0] + sigma * (d_vec[1])
        for j in range(1, len(temp) - 1):
            temp[j] = (2 - 2 * sigma) * d_vec[j] + sigma * (d_vec[j - 1] + d_vec[j + 1])
        temp[-1] = (2 - 2 * sigma) * d_vec[j] + sigma * (d_vec[j - 1])
        d_vec = TridiagSolve(a_vec, b_vec, c_vec, temp)
        temp = np.zeros(len(temp))
        grid[i, 1:-1] = d_vec

    print(grid.shape)
    np.savetxt("cn_grid.txt", grid, fmt='%.5f')
    result = grid[-1]

    return result

def main():

    T = 0.01
    N = 25
    K = 50
    a = 1

    l = 1
    tau = T / K
    h = l / N

    sigma = (a * a * tau) / (h * h)
    print("Sigma:\n", sigma)

    x = np.arange(0, l + h, h)
    t = np.arange(0, T + tau, tau)

    print(x)
    print(t)

    analitical = ExactSolution(x, t, a)
    print("Analitical Solution:\n", analitical)

    explicit_fd = ExFinDiff(x, t, sigma)
    print("Explicit Finite Difference:\n", explicit_fd)

    implicit_fd = ImFinDiff(x, t, sigma)
    print("Implicit Difference:\n", implicit_fd)

    crank_nicol = CrankNicolson(x, t, sigma)
    print("Crank-Nicolson:\n", crank_nicol)

    plt.plot(analitical, label = "Analitical Solution", color="lime")
    plt.plot(explicit_fd, label = "Explicit Solution", color="red")
    plt.plot(implicit_fd, label = "Implicit Solution", color="blue")
    plt.plot(crank_nicol, label = "Crank-Nicolson", color="magenta")
    plt.autoscale(True)
    plt.legend()
    plt.show()

    return


main()