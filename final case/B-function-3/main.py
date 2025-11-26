import numpy as np
import random as rand
import cvxpy as cp

import math

from coefficients_B import *
from Conterexample import *

def gridstate(Range0, Range1, Range2, Range3, size):
    low0 = Range0[0]
    low1 = Range1[0]
    low2 = Range2[0]
    low3 = Range3[0]
    l = []
    while (low0 <= Range0[1] - size):
        # Reset low2
        low1 = Range1[0]
        while (low1 <= Range1[1] - size):
            # Reset low3
            low2 = Range2[0]
            while (low2 <= Range2[1] - size):
                # Reset low4
                low3 = Range3[0]
                while (low3 <= Range3[1] - size):
                    l.append([low0, low1, low2, low3])
                    low3 = low3 + size
                low2 = low2 + size
            low1 = low1 + size
        low0 = low0 + size
    return (l)


def gridinput(U, size):
    low0 = U[0]
    low1 = U[0]
    low2 = U[0]
    low3 = U[0]
    l = []
    while (low0 <= U[1] - size):
        # Reset low2
        low1 = U[0]
        while (low1 <= U[1] - size):
            # Reset low3
            low2 = U[0]
            while (low2 <= U[1] - size):
                # Reset low4
                low3 = U[0]
                while (low3 <= U[1] - size):
                    l.append([low0, low1, low2, low3])
                    low3 = low3 + size
                low2 = low2 + size
            low1 = low1 + size
        low0 = low0 + size
    return l

K = 3
delta = 1

epsilon_l = 0.001

# lower and upper bounds
X_min = 20
X_max = 30
X_R = [X_min, X_max]

U = [0, 1]

ul = U[0]
uu = U[1]

inter_max=50

# Fault
Xf = [24, 26, X_min, X_max]
# Xf_2=[20,26]

# Sample from state and input
i=0
X_dis = gridstate(X_R, X_R, X_R, X_R, 2)  # l.append([low0, low1, low2, low3])

U_dis = gridinput(U, 0.5) # l.append([low0, low1, low2, low3])

# Actions in DFA
Act = ['sig_1', 'sig_2', 'sig_3', 'sig_4', 'sig_5', 'sig_6']

# states of DFA
DFA = dfa_state(K)



while 1 :
    con, coe = com_coff(X_dis, U_dis,DFA,K)

    objective = cp.Minimize(1)
    prob = cp.Problem(objective, con)
    prob.solve(solver=cp.SCS)

    print("************************************************* \n")
    print("The ", i, "th iteration:\n")
    print("status:", prob.status)
    if (prob.status == cp.INFEASIBLE):
        print("The barrier with the given template doesn’t exist. \n")
        break

    else:
        outval = print_coff(prob, coe,K)
        # find counterexamples
        counter_s, counter_u, flag=compute_counter(outval)
        if flag==0:
            break
        else:
            X_dis.append(counter_s)
            U_dis.append(counter_u)
            i=i+1
    if i>inter_max:
        print("Last interation!")
        break


