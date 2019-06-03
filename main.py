#!/usr/bin/env python3
import scipy.special as special
from scipy.integrate import quad
from numpy import sqrt, power as pow, log, exp
from typing import Callable

# return information given theta, a, b, and c
def information(theta: float, a: float, b: float, c: float) -> float:
    numer = (2.89 * pow(a, 2)) * (1 - c)
    denom1 = c + exp( (1.7 * a) * (theta - b))
    denom2 = pow(1 + exp( (-1.7 * a) * (theta-b)), 2)
    return (numer) / (denom1 * denom2)


# returns the theta* value
def theta_star(a: float, b: float, c: float) -> float:
    return (b + (1/1.7) * log((1 + sqrt(1 + 8*c))/2))

def ebi_weight(theta: float, a: float, b: float, c:float, std_err: float) -> float:
    lb = theta - 2 * std_err
    ub = theta + 2 * std_err
    return quad(information, lb, ub, args=(a, b, c))


def ebi(theta_hat: float, a: float, b: float, c: float, std_err: float) -> float:
    ts = theta_star(a, b, c)
    weight = ebi_weight(theta_hat, a, b, c, std_err)[0]
    return (1 + (1/information(ts, a, b, c))) * weight




def main():
    it = float( input("θ:\t"))
    a = float( input("a:\t"))
    b = float( input("b:\t"))
    c = float( input("c:\t"))
    se = float (input("standard error:\t"))
    ts = theta_star(a, b, c)
    weight = ebi_weight(it, a, b, c, se)
    ebi_ = ebi(it, a, b, c, se)


    print(f"θ*:\t{ts}\n")
    print(f"weight:\t\t\t{weight[0]}\n")
    print(f"integration error:\t{weight[1]}")
    print(f"\nEBI:\t{ebi_}")
    exit()


main()
