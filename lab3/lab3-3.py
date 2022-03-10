import sys
import numpy as np
import matplotlib.pyplot as plt

def LeastSquares(points, type):
    n = len(points[0])
    coef = np.zeros(type + 1)
    system_coef = np.zeros((type + 1, type + 1))
    b_vec = np.zeros(type + 1)

    if type == 1:
        system_coef[0][0] = n
        for i in range(n):
            b_vec[0] += points[1][i]
            b_vec[1] += points[0][i] * points[1][i]

            system_coef[0][1] += points[0][i]
            system_coef[1][1] += points[0][i] * points[0][i]
        system_coef[1][0] = system_coef[0][1]

        coef = np.linalg.solve(system_coef, b_vec)
    else:
        system_coef[0][0] = n
        for i in range(n):
            b_vec[0] += points[1][i]
            b_vec[1] += points[0][i] * points[1][i]
            b_vec[2] += points[0][i] * points[0][i] * points[1][i]

            system_coef[0][1] += points[0][i]
            system_coef[0][2] += points[0][i] * points[0][i]
            system_coef[1][2] += points[0][i] * points[0][i] * points[0][i]
            system_coef[2][2] += points[0][i] * points[0][i] * points[0][i] * points[0][i]
        system_coef[1][0] = system_coef[0][1]
        system_coef[1][1] = system_coef[0][2]
        system_coef[2][0] = system_coef[1][1]
        system_coef[2][1] = system_coef[1][2]

        coef = np.linalg.solve(system_coef, b_vec)

    return coef

def Polynomial(coef, points, type):
    result = np.zeros(len(points[0]))

    if type == 1:
        for i in range(len(points[0])):
            result[i] = coef[0] + coef[1] * points[0][i]
    else:
        for i in range(len(points[0])):
            result[i] = coef[0] + coef[1] * points[0][i] + coef[2] * points[0][i] * points[0][i]

    return result

def RSS(fx, points):
    rss = 0

    for i in range(len(fx)):
        rss += (fx[i] - points[1][i]) * (fx[i] - points[1][i])

    return rss

def main():
    # Read points from file
    #----------------------
    points = np.zeros((2, 6))
    with open('lab3-3_points.txt') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            points[i] = np.fromstring(lines[i], dtype=float, sep=' ')
    print("\nThe points of the function:")
    print(points)

    plt.title("First and Second order polynomials")
    plt.xlabel("x") 
    plt.ylabel("y")
    plt.grid()
    plt.plot(points[0], points[1], label = "The given function", color="royalblue")

    # First order polynomial
    #-----------------------
    coef = LeastSquares(points, 1)
    print("\nCoefficients of first order polynomial:")
    print(coef)

    fx = Polynomial(coef, points, 1)
    print("\nValues of polynomial in given points")
    print(fx)

    rss = RSS(fx, points)
    print("\nThe residual sum of squares:")
    print(rss)

    plt.plot(points[0], fx, label = "First order polynomial", color="orangered")

    # Second order polynomial
    #------------------------
    coef = LeastSquares(points, 2)
    print("\nCoefficients of second order polynomial:")
    print(coef)

    fx = Polynomial(coef, points, 2)
    print("\nValues of polynomial in given points")
    print(fx)

    rss = RSS(fx, points)
    print("\nThe residual sum of squares:")
    print(rss)
  
    plt.plot(points[0], fx, label = "Second order polynomial", color="lime")
    plt.autoscale(True)
    plt.legend()
    plt.show()
    
main()

# -5.0 -3.0 -1.0 1.0 3.0 5.0
# -2.0558 -0.18016 1.3562 1.7854 3.3218 5.1974