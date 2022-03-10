import numpy as np
np.set_printoptions(formatter={'float': lambda x: "{0:0.5f}".format(x)}, linewidth=200)


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

def main():
    a_vec = np.array([-5,-5,-9,1])
    b_vec = np.array([8,22,-11,-15,7])
    c_vec = np.array([4,8,1,1])
    d_vec = np.array([48,125,-43,18,-23])

    roots = TridiagSolve(a_vec, b_vec, c_vec, d_vec)

    print(roots)

main()