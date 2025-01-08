#   cdate:      2024-12-28 00:02
#   mdate:      2024-12-28 10:15
#   author:     209818061
#   desc:       Homework assignment 8, implementation of interpolation polynomial

from math import sin, pi
# TODO: import here additional Python modules, as you desire

# TODO: define your utility functions here

# TODO: DO NOT DEFINE NEW FUNCTIONS BELOW
def eqpts(a: float, b: float, n: int) -> list:
    """
    Generate a list of n+1 equidistant points in the interval [a, b], starting with a and ending with b.
    Returns a list with n+1 equidistant points.
    """
    out = []
    
    diff = (b - a) / (n - 1)
    for i in range(n + 1):
        out.append(a + i * diff)

    return out

def ipn(x: list, y: list, x_eval: list) -> list:
    """
    Evaluate the polynomial interpolating the samples (x[k], y[k]) (k = 0, ..., len(x)-1)
    using Newton's basis.

    x: a list of floats; interpolation points.
    y: a list of floats; values on interpolation points.
    x_eval: a list of floats; points on which to evaluate the interpolating polynomial.

    x and y must have the same length.

    Returned value is a list of floats, such that 
        out[j] = p(x_eval[j]) for j = 0, ..., len(x_eval)-1
    where p(x) is the interpolation polynomial satisfying
        y[k] = p(x[k]) for k = 0, ..., len(x)-1.
    TODO interpoation polynomial Newton form
    """
    if not (len(x) == len(y)):
        raise ValueError("x and y must be lists of identical lengths.")

    n = len(x)
    divided_diff = [[0] * n for _ in range(n)]

    for i in range(n):
        divided_diff[i][0] = y[i]

    for j in range(1, n):
        for i in range(n - j):
            divided_diff[i][j] = (divided_diff[i + 1][j - 1] - divided_diff[i][j - 1]) / (x[i + j] - x[i])

    out = []
    for x0 in x_eval:
        result = divided_diff[0][0]
        term = 1
        for i in range(1, n):
            term *= (x0 - x[i - 1])
            result += divided_diff[0][i] * term
        out.append(result)

    return out


def ipl(x: list, y: list, x_eval: list) -> list:
    """
    Evaluate the polynomial interpolating the samples (x[k], y[k]) (k = 0, ..., len(x)-1)
    using Lagrange basis.

    x: a list of floats; interpolation points.
    y: a list of floats; values on interpolation points.
    x_eval: a list of floats; points on which to evaluate the interpolating polynomial.

    x and y must have the same length.

    Returned value is a list of floats, such that 
        out[j] = p(x_eval[j]) for j = 0, ..., len(x_eval)-1
    where p(x) is the interpolation polynomial satisfying
        y[k] = p(x[k]) for k = 0, ..., len(x)-1.
    """
    if not (len(x) == len(y)):
        raise ValueError("x and y must be lists of identical lengths.")

    n = len(x)

    out = []
    for x0 in x_eval:
        p_x0 = 0
        for i in range(n):
            L_i = 1
            for j in range(n):
                if i != j:
                    L_i *= (x0 - x[j]) / (x[i] - x[j])
            p_x0 += y[i] * L_i
        out.append(p_x0)

    return out

if __name__=="__main__":
    n = [6, 30]
    ip = {"Lagrange": ipl, "Newton": ipn}
    x_eval = [pi/4]
    for n_ in n:
        x = eqpts(0, pi, n_)
        y = [sin(x_) for x_ in x]
        for ip_ in ip.keys():
            approx = ip[ip_](x, y, x_eval)[0]
            err = abs(approx - sin(x_eval[0]))
            print(f"Method: {ip_:10}, n={n_:2d}:\tsin(pi/4) = {sin(x_eval[0])},\tp(pi/4) = {approx},\t|p(pi/4) - sin(pi/4)|={err:.4e}.")
        
