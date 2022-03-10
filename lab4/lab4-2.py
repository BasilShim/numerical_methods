import numpy as np
np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)}, linewidth=200)

def PFunction(x):
    result = -(2 * x + 4) / (x * (x + 4))
    return result

def QFunction(x):
    result = 2 / (x * (x + 4))
    return result

def FiniteDifference(a, b, h, y, y_right):
    n = int((b - a) / h)
    x = np.arange(a, b + h, h)
    coef = np.zeros((n, n))
    d_vec = np.zeros(n)

    coef[0][0] = -2 + h * h * QFunction(x[1])
    coef[0][1] = 1 + PFunction(x[1]) * h / 2
    d_vec[0] = (-1 + h * PFunction(x[1]) / 2) * y
    for i in range(1, n - 1):
        for j in range(1, n - 1):
            d_vec[i] = 0
            if i == j:
                coef[i][i] = -2 + h * h * QFunction(x[i + 1])
                coef[i][i - 1] = (1 - h * PFunction(x[i + 1]) / 2)
                coef[i][i + 1] = (1 + h * PFunction(x[i + 1]) / 2)
    coef[n - 1][n - 1] = -2 + h * h * QFunction(x[n - 2])
    coef[n - 1][n - 2] = (1 - h * PFunction(x[n - 2]) / 2)
    d_vec[n - 1] = (-1 - h * PFunction(x[n - 2]) / 2) * y_right
    y_vec = np.array([])
    y_vec = np.append(y_vec, y)
    y_vec = np.append(y_vec, np.linalg.solve(coef, d_vec))
    #y_vec = np.append(y_vec, y_right)
    
    return y_vec

def ExactSolution(x):
    result = x * x + x + 2
    return result

def ExactResult(a, b, h):
    result = np.array([])
    for i in np.arange(a, b + h, h):
        result = np.append(result, ExactSolution(i))

    return result

def Error(result, exact_result):
    error = np.absolute(result - exact_result)

    return error

def RungeRombergRichardson(a, b, h, y, y_right, exact_result):
    step = FiniteDifference(a, b, h, y, y_right)
    halfstep = FiniteDifference(a, b, h / 2, y, y_right)
    result = step

    for i in range(len(np.arange(a, b + h, h))):
        result[i] = step[i] + (step[i] - halfstep[i * 2]) / (0.5 * 0.5 - 1)

    print("\nFinite Difference method with h={:.2f} using Runge-Romberg-Richardson method for accuracy: ".format(h))
    print(result)
    print(" Error:")
    print(Error(result, exact_result))
    return None

def main():
    a = 0
    b = 2
    h = 0.1
    y = ExactSolution(a)
    y_right = ExactSolution(b)

    exact_result = ExactResult(a, b, h)
    print("\nThe results of exact solution:")
    print(np.arange(a, b + h, h))
    print(exact_result)
    print("The results of Finite Difference method:")
    result = FiniteDifference(a, b, h, y, y_right)
    print(result)
    RungeRombergRichardson(a, b, h, y, y_right, exact_result)

main()