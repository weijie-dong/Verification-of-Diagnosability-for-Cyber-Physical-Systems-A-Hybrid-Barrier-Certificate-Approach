import numpy as np
import random as rand
import cvxpy as cp
import time
import math

# import eventlet
from func_timeout import func_set_timeout
import func_timeout



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

print("Starting the problem: \n")
K = 3
delta = 0.5
# K = 3
# delta = 1

epsilon_l = 0.001

# lower and upper bounds
X_min = 15
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
print("Start the sampling: \n")
i=0
sam=1
X_dis = gridstate(X_R, X_R, X_R, X_R, sam)  # l.append([low0, low1, low2, low3])

U_dis = gridinput(U, 0.5) # l.append([low0, low1, low2, low3])

# Actions in DFA
Act = ['sig_1', 'sig_2', 'sig_3', 'sig_4', 'sig_5', 'sig_6']

# states of DFA
DFA = dfa_state(K)

print("Start the compute coff: \n")
con, coe = com_coff(X_dis, U_dis,DFA,K)


print("Start the optimazation: \n")
objective = cp.Minimize(1)
prob = cp.Problem(objective, con)
prob.solve(solver=cp.SCS)

outval = []
for i in range(K + 3):
    outval.append(coe[i][0].value)
    # print(outval[i])
    
print("The final barrier certificate is: ")

print(outval)

txtfile='data-V.txt'
with open(txtfile, 'w') as file:
     for item in outval:
        print(item)
        for j in item:
            print(j)
            file.write(str(j)+'\t')
        file.write('\n')

# print_outval(outval,K)

### Counter example
###################################################################
# while 1 :
#     con, coe = com_coff(X_dis, U_dis,DFA,K)

#     objective = cp.Minimize(1)
#     prob = cp.Problem(objective, con)
#     prob.solve(solver=cp.SCS)
    
#     # print("************************************************* \n")
#     # print("The ", i, "th iteration:\n")
#     # print("status:", prob.status)
#     if (prob.status == cp.INFEASIBLE):
#         print("The barrier with the given template doesn’t exist. \n")
#         break

#     else:
#         outval = []
#         for i in range(K + 3):
#             outval.append(coe[i][0].value)
#             # print(outval[i])

#         try:
#             counter_s, counter_u, flag = compute_process(outval)
#             if flag == 0:
#                 print("The final barrier certificate is: \n")
#                 print_outval(outval,K)
#                 break
#             else:
#                 X_dis.append(counter_s)
#                 U_dis.append(counter_u)
#                 i = i + 1
#         except func_timeout.exceptions.FunctionTimedOut:
#             print("Counterexample cannot be found!")
#             print("The final barrier certificate is: ")
#             print_outval(outval,K)
#             break

#     if i>inter_max:
#         print("Last interation!")
#         print("The final barrier certificate is: \n")
#         print_outval(outval,K)
#         break
