import numpy as np
np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)})

def CubicSpline(points, point):
    n = len(points[0])
    coef = np.zeros((n - 2, n - 2))
    c_coef = np.zeros(n - 3)
    d_coef = np.zeros(n - 2)
    value = 0

    for i in range(2, n):
        for j in range(2, n):
            h1 = points[0][i - 1] - points[0][i - 2]
            h2 = points[0][i] - points[0][i - 1]
            d_coef[i - 2] = 3 * ((points[1][i] - points[1][i - 1]) / h2 - (points[1][i - 1] - points[1][i - 2]) / h1)
            if i == j:
                coef[i - 2][i - 2] = 2 * (h1 + h2)
                if i < n - 1:
                    coef[i - 2][i - 1] = h1
                    coef[i - 1][i - 2] = h2
        
    c_coef = np.linalg.solve(coef, d_coef)

    coef = np.zeros((n - 1, n - 1))
    coef[2][1:] = c_coef
    coef[0] = points[1][:n - 1]

    for i in range(n - 1):
        h1 = points[0][i + 1] - points[0][i]
        if i == n - 2:
            coef[1][i] = (points[1][i + 1] - points[1][i]) / h1 - 2 * h1 / 3 * coef[2][i]
            coef[3][i] = - coef[2][i] / (3 * h1)
            break
        coef[1][i] = (points[1][i + 1] - points[1][i]) / h1 - h1 / 3 * (coef[2][i + 1] + 2 * coef[2][i])
        coef[3][i] = (coef[2][i + 1] - coef[2][i]) / (3 * h1)
    
    for i in range(n - 2):
        if point < points[0][i + 1] and point > points[0][i]:
            value = coef[0][i] + coef[1][i] * (point - points[0][i]) + coef[2][i] * (point - points[0][i]) * (point - points[0][i]) + coef[3][i] * (point - points[0][i]) * (point - points[0][i]) * (point - points[0][i])
            print("{:f} {:f} {:f} {:f}".format(coef[0][i], coef[1][i], coef[2][i], coef[3][i]))


    return coef, value

def main():
    # Read points from file
    #----------------------
    points = np.zeros((2, 5))
    with open('lab3-2_points.txt') as f:
        lines = f.readlines()
        for i in range(len(lines) - 1):
            points[i] = np.fromstring(lines[i], dtype=float, sep=' ')
        point = np.float64(lines[len(lines) - 1])
    print("\nThe points of the function:")
    print(points)

    # Do the magic
    #-------------
    coef, value = CubicSpline(points, point)
    print("\nThe matrix of coefficients:")
    print(coef)
    print("\nValue of cubic spline in given point:")
    print("{:.5f}".format(value))


main()