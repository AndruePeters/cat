#!/usr/bin/env python3

# Author: Andrue Peters
# Date: 6/2/19
# Description:  Develops test data for evaluating EBI for 3 parameter information function
#               Allows user to vary theta, a, b, c, and the standard error for the ebi function
#                   and outputs the information to a file
#
#               Follow the menu upon running.


import scipy.special as special
from scipy.integrate import quad
from numpy import sqrt, power as pow, log, exp, linspace
from typing import Callable
import random as rand
import matplotlib.pyplot as plt

#########################################################################
#       Returns the information of the parameters                       #
#########################################################################
def information(theta: float, a: float, b: float, c: float) -> float:
    numer = (2.89 * pow(a, 2)) * (1 - c)
    denom1 = c + exp( (1.7 * a) * (theta - b))
    denom2 = pow(1 + exp( (-1.7 * a) * (theta-b)), 2)
    return (numer) / (denom1 * denom2)


#########################################################################
#       Returns theta*                                                  #
#########################################################################
def theta_star(a: float, b: float, c: float) -> float:
    return (b + (1/(1.7 * a)) * log((1 + sqrt(1 + 8*c))/2))

#########################################################################
#       Calculates the EBI weight using Guassian quadrature             #
#########################################################################
def ebi_weight(theta: float, a: float, b: float, c:float, std_err: float) -> float:
    lb = theta - 2 * std_err
    ub = theta + 2 * std_err
    return quad(information, lb, ub, args=(a, b, c))

#########################################################################
#       Returns the final EBI                                           #
#########################################################################
def ebi(theta_hat: float, a: float, b: float, c: float, std_err: float) -> float:
    ts = theta_star(a, b, c)
    weight = ebi_weight(theta_hat, a, b, c, std_err)[0]
    info_ts = information(ts, a, b, c)
    if info_ts == 0:
        info_ts = 0.00000000001
    return (1 + (1/info_ts)) * weight


#########################################################################
#       Holds all values but theta constant and generates test data     #
#########################################################################
def vary_theta(file_name: str):
    print("θ varies\n")
    lb = float(input("Lower bound:\t"))
    ub = float(input("Upper bound θ:\t"))
    num = int(input ("Number of steps:\t"))
    a = float( input("a:\t"))
    b = float( input("b:\t"))
    c = float( input("c:\t"))
    se = float (input("std err:\t"))
    vals = linspace(lb, ub, num,endpoint=True)
    f = open(file_name, "w")
    write_header(f)
    for x in vals:
        info = information(x, a, b, c)
        eb = ebi(x, a, b, c, se)
        f.write(f"{x},{a},{b},{c},{se},{info},{eb}\n")
    f.close()

#########################################################################
#       Holds all values but a constant and generates test data         #
#########################################################################
def vary_a(file_name: str):
    print("a varies\n")
    theta = float(input("θ:\t"))
    lb = float(input("Lower bound:\t"))
    ub = float(input("Upper bound:\t"))
    num = int(input("Number of steps:\t"))
    b = float( input("b:\t"))
    c = float( input("c:\t"))
    se = float (input("std err:"))
    vals = linspace(lb, ub, num,endpoint=True)
    f = open(file_name, "w")
    write_header(f)
    for x in vals:
        info = information(theta, x, b, c)
        eb = ebi(theta, x, b, c, se)
        f.write(f"{theta},{x},{b},{c},{se},{info},{eb}\n")
    f.close()


#########################################################################
#       Holds all values but b constant and generates test data         #
#########################################################################
def vary_b(file_name: str):
    print("b varies\n")
    theta = float(input("θ:\t"))
    a = float( input("a:\t"))
    lb = float(input("Lower bound:\t"))
    ub = float(input("Upper bound:\t"))
    num = int(input("Number of steps:\t"))
    c = float( input("c:\t"))
    se = float (input("std err:"))
    vals = linspace(lb, ub, num,endpoint=True)
    f = open(file_name, "w")
    write_header(f)
    for x in vals:
        info = information(theta, a, x, c)
        eb = ebi(theta, a, x, c, se)
        f.write(f"{theta},{a},{x},{c},{se},{info},{eb}\n")
    f.close()

#########################################################################
#       Holds all values but c constant and generates test data         #
#########################################################################
def vary_c(file_name: str):
    print("c varies\n")
    theta = float(input("θ:\t"))
    a = float( input("a:\t"))
    b = float( input("b:\t"))
    lb = float(input("Lower bound:\t"))
    ub = float(input("Upper bound:\t"))
    num = int(input("Number of steps:\t"))
    se = float (input("std err:"))
    vals = linspace(lb, ub, num,endpoint=True)
    f = open(file_name, "w")
    write_header(f)
    for x in vals:
        info = information(theta, a, b, x)
        eb = ebi(theta, a, b, x, se)
        f.write(f"{theta},{a},{b},{x},{se},{info},{eb}\n")
    f.close()

#########################################################################
#       Holds all values but standard error constant                    #
#       and generates test data                                         #
#########################################################################
def vary_se(file_name: str):
    print("std err varies\n")
    theta = float(input("θ:\t"))
    a = float( input("a:\t"))
    b = float( input("b:\t"))
    c = float( input("c:\t"))
    lb = float(input("Lower bound:\t"))
    ub = float(input("Upper bound:\t"))
    num = int(input("Number of steps:\t"))
    vals = linspace(lb, ub, num,endpoint=True)
    f = open(file_name, "w")
    write_header(f)
    for x in vals:
        info = information(theta, a, b, c)
        eb = ebi(theta, a, b, c, x)
        f.write(f"{theta},{a},{b},{c},{x},{info},{eb}\n")
    f.close()


#########################################################################
#       Holds theta, and se, but varries a, b, and c                    #
#########################################################################
def vary_a_b_c(file_name: str):
    print("Parameters a, b, and c all vary")
    theta = float( input("θ:\t"))
    se = float (input("std_err:\t"))
    lb_a = float (input ("Lower bound of a:\t"))
    ub_a = float (input ("Upper bound of a:\t"))
    lb_b = float (input ("Lower bound of b:\t"))
    ub_b = float (input ("Upper bound of b:\t"))
    lb_c = float (input ("Lower bound of c:\t"))
    ub_c = float (input ("Upper bound of c:\t"))
    num = int( input ("Number of steps:\t"))
    vals_a = linspace(lb_a, ub_a, num, endpoint=True)
    vals_b = linspace(lb_b, ub_b, num, endpoint=True)
    vals_c = linspace(lb_c, ub_c, num, endpoint=True)
    f = open(file_name, "w")
    write_header(f)
    for a in vals_a:
        for b in vals_b:
            for c in vals_c:
                info = information(theta, a, b, c)
                eb = ebi(theta, a, b, c, se)
                f.write(f"{theta},{a},{b},{c},{se},{info},{eb}\n")
    f.close()


#########################################################################
#       Holds theta, and se, but varries a, b, and c randomly           #
#########################################################################
def vary_a_b_c_rand(file_name: str):
    print("Parameters a, b, and c are randomly generated using normal dist.")
    theta = float ( input("θ:\t"))
    se = float ( input("std_err:\t"))
    mu_a = float ( input("mean a:\t"))
    sigma_a = float ( input("std_dev a:\t"))
    mu_b = float (input ("mean b:\t"))
    sigma_b = float (input ("std_dev b:\t"))
    mu_c = float ( input ("mean c:\t"))
    sigma_c = float ( input("std_dev c:\t"))
    num = int ( input("Number of test points:\t"))
    f = open(file_name, "w")
    write_header(f)
    for i in range(num):
        a = rand.gauss(mu_a, sigma_a)
        b = rand.gauss(mu_b, sigma_b)
        c = rand.gauss(mu_c, sigma_c)
        info = information(theta, a, b, c)
        eb = ebi(theta, a, b, c, se)
        f.write(f"{theta},{a},{b},{c},{se},{info},{eb}\n")
    f.close()


#########################################################################
#       Helper function. Writes the first line of csv file              #
#       *** Assumes file_obj has already been opened                    #
#########################################################################
def write_header(file_obj):
    file_obj.write("Theta, a, b, c, std err, info, ebi\n")

#########################################################################
#       Helper function. Prints the menu.                               #
#########################################################################
def print_menu():
    print("\n\nEnter the the menu item number.\n")
    print("1. θ varies")
    print("2. a varies")
    print("3. b varies")
    print("4. c varies")
    print("5. std_err varies")
    print("6. a,b, and c vary")
    print("7. a, b, and c vary randomly on a normal distribution")
    print("0. exit")

def main():
    choice = 1
    while choice != 0:
        print_menu()
        choice = int(input("menu number:\t"))
        file_name = input("file name (include .csv):\t")
        if choice == 0:
            return
        if choice == 1:
            vary_theta(file_name)
        elif choice == 2:
            vary_a(file_name)
        elif choice == 3:
            vary_b(file_name)
        elif choice == 4:
            vary_c(file_name)
        elif choice == 5:
            vary_se(file_name)
        elif choice == 6:
            vary_a_b_c(file_name)
        elif choice == 7:
            vary_a_b_c_rand(file_name)
        else:
            print(f"{choice} is an invalid response")
        print("Completed.")

"""plt.plot([1,2,3,4],[1,2,3,4])
plt.ylabel('some numbers')
plt.show()"""
main()
