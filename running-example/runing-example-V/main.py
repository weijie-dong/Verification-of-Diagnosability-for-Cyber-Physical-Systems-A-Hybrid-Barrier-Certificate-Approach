import numpy as np
import random as rand
from z3 import *
import cvxpy as cp
import math


def tran(q, u):  # q=x , u=u
    ind = int(q)
    if ind == 1:
        if u == 'a':
            temp = '2'
        elif u == 'b':
            temp = '6'
        else:
            print('input error')
    elif ind < 5 and ind >= 2:
        if u == 'c':
            temp = str(ind + 1)
        else:
            print('input error')
    elif ind < 9 and ind >= 6:
        if u == 'c':
            temp = str(ind + 1)
        else:
            print('input error')
    elif ind == 5 or ind == 9:
        temp = q
    return temp




X_min = 0
X_max = 10
X_R = [X_min, X_max]
X_R_b = [X_min - 30, X_max + 30]

# # degree: 0:3
# def monomial(par, par_h):
#     # 1
#     #      x1
#     #      x2
#     #    x1^2
#     #   x1*x2
#     #    x2^2
#     #    x1^3
#     # x1^2*x2
#     # x1*x2^2
#     #    x2^3
#     x1 = par
#     x2 = par_h
#     mon = [1, x1, x2, x1 ** 2, x1 * x2, x2 ** 2, x1 ** 3, (x1 ** 2) * x2, x1 * (x2 ** 2), x2 ** 3]  # 1,10
#     # for i in range(m+1):
#     #     mon.append(pow(x,i))
#     C = np.array(mon)
#     return C


# degree: 0:4
def monomial(par,par_h,mode):
    if mode==4:
    #     1
    #     x1
    #     x2
    # x1 ^ 2
    # x1 * x2
    # x2 ^ 2
    # x1 ^ 3
    # x1 ^ 2 * x2
    # x1 * x2 ^ 2
    # x2 ^ 3
    # x1 ^ 4
    # x1 ^ 3 * x2
    # x1 ^ 2 * x2 ^ 2
    # x1 * x2 ^ 3
    # x2 ^ 4
        x1=par
        x2 = par_h
        mon=[1, x1, x2, x1**2, x1*x2, x2**2, x1**3, (x1**2)*x2, x1*(x2**2), x2**3,
             x1 ** 4, x1 ** 3 * x2, x1 ** 2 * x2 ** 2, x1 * x2 ** 3, x2 ** 4] # 1,15
        # for i in range(m+1):
        #     mon.append(pow(x,i))
        C = np.array(mon)
        return C
    elif mode==3:


#     # 1
#     #      x1
#     #      x2
#     #    x1^2
#     #   x1*x2
#     #    x2^2
#     #    x1^3
#     # x1^2*x2
#     # x1*x2^2
#     #    x2^3
        x1=par
        x2 = par_h
        mon=[1, x1, x2, x1**2, x1*x2, x2**2, x1**3, (x1**2)*x2, x1*(x2**2), x2**3] # 1,10
        # for i in range(m+1):
        #     mon.append(pow(x,i))
        C = np.array(mon)
        return C
    elif mode==2:


    #     # 1
    #     #      x1
    #     #      x2
    #     #    x1^2
    #     #   x1*x2
    #     #    x2^2
        x1=par
        x2 = par_h
        mon=[1, x1, x2, x1**2, x1*x2, x2**2] # 1,6
        # for i in range(m+1):
        #     mon.append(pow(x,i))
        C = np.array(mon)
        return C

def dfa_state(k):
    DFA = []
    for i in range(k + 1):
        DFA.append(str(i))
    DFA.append('trap')
    DFA.append('F')
    return DFA


# DFA = ['0', '1', '2', '3', 'trp', 'F']

def statepair(q):  # q is the state set of system rather than DFA
    l = []
    for i in q:
        for j in q:
            l.append([i, j])
    return l


def inputpair(U):  # q is the state set
    l = []
    for i in U:
        for j in U:
            l.append([i, j])
    return l


def P1(x, delta):  # x=[x1,x2], x1:str
    ob = output(x[0]) - output(x[1])
    b = math.pow(ob, 2) <= math.pow(delta, 2)
    return b


def P2(x, xf):  # x=[x1,x2], x1:str, Xf = [0, 0.8]
    b = (output(x[0]) <= xf[1]) and (output(x[0]) >= xf[0])
    return b


def P3(x, xf):  # x=[x1,x2], x1:str, Xf = [0, 0.8]
    b = (output(x[1]) <= xf[1]) and (output(x[1]) >= xf[0])
    return b


def R(x, X_R):  # x=[x1,x2], x1:str, X_R = [X_min, X_max]
    b = (output(x[0]) <= X_R[1]) and (output(x[0]) >= X_R[0]) \
        and (output(x[1]) <= X_R[1]) and (output(x[1]) >= X_R[0])
    return b


def sigma(x, X_R, xf, delta):  # x=[x1,x2], x1:str
    act = []
    if (P1(x, delta) and R(x, X_R) and (not P2(x, xf)) and (not P3(x, xf))):
        act.append('sig_1')
    if ((R(x, X_R) and (not P1(x, delta))) or (P3(x, xf) and P1(x, delta))):
        act.append('sig_2')
    if (P1(x, delta) and P2(x, xf) and (not P3(x, xf)) and R(x, X_R)):
        act.append('sig_3')
    return act


def q_0(sig):
    if (sig == 'sig_1'):
        q_n = '0'
    elif (sig == 'sig_2'):
        q_n = 'trap'
    elif (sig == 'sig_3'):
        q_n = '1'
    else:
        print('sigma in 0 is wrong')
    return (q_n)


def q_t(sig):
    if (sig == 'sig_1' or sig == 'sig_2' or sig == 'sig_3'):
        q_n = 'trap'
    else:
        print('sigma in trap is wrong')
    return (q_n)


def q_F(sig):
    if (sig == 'sig_1' or sig == 'sig_2' or sig == 'sig_3'):
        q_n = 'F'
    else:
        print('sigma in F is wrong')
    return (q_n)


def q_i(q, sign, K):
    ind = q2int(q, K)
    if (sign == 'sig_2'):
        q_n = 'trap'
    elif ind < K and (sign == 'sig_1' or sign == 'sig_3'):
        ind = ind + 1
        q_n = str(ind)
    elif (ind == K and (sign == 'sig_1' or sign == 'sig_3')):
        q_n = 'F'
    else:
        print('sig_i error')
    return q_n


def q2int(q, k):
    if (q == 'trap'):
        temp = k + 1
    elif (q == 'F'):
        temp = k + 2
    else:
        temp = int(q)
    return temp


def next(q, sig, K):
    if (q == '0'):
        q_n = q_0(sig)
    elif (q == 'trap'):
        q_n = q_t(sig)
    elif (q == 'F'):
        q_n = q_F(sig)
    elif (q2int(q, K) <= K and q2int(q, K) >= 1):
        q_n = q_i(q, sig, K)
    return (q_n)


# Obtain the values out of X_R: X_R_b[0] < X_R[0] < X_R[1] < X_R_b[1]
def outregion(X_R,X_R_b,siz):
    l=[]

    i=X_R[0]
    while(i <= X_R[0] and i >=X_R_b[0]):
        i=i-siz
        j = X_R[0]
        while(j <= X_R[0] and j >=X_R_b[0]):
            j = j - siz
            l.append([i,j])
        j = X_R[1]
        while (j >= X_R[1] and j <= X_R_b[1]):
            j = j + siz
            l.append([i, j])

    i = X_R[1]
    while (i >= X_R[1] and i <= X_R_b[1]):
        i = i + siz
        j = X_R[1]
        while (j >= X_R[1] and j <= X_R_b[1]):
            j = j + siz
            l.append([i, j])

        j = X_R[0]
        while (j <= X_R[0] and j >= X_R_b[0]):
            j = j - siz
            l.append([i, j])

    return l





K = 2
delta = 1

# Parameters

# epsilon_1=0;





U = ['a', 'b', 'c'];
# U2_all = (u2-U_min)*(U_max-u2);

# Fault: less than 0.8
Xf = [1, 1.5]

def output(q):
    # out = [0.5,
    #        1.3, 2.8, 6, 7,
    #        1.8, 3, 6, 10]
    # out = [0,
    #        2, 2.8, 6, 6.5,
    #        2.5, 3, 7, 9]
    out = [0,
           1.2, 3.2, 5.2, 7.2,
           2.2, 4.2, 6.2, 9]
    ind = int(q) - 1
    return out[ind]

System = ['1', '2', '3', '4', '5', '6', '7', '8', '9']

# Actions in DFA, new version
Act = ['sig_1', 'sig_2', 'sig_3']

DFA = dfa_state(K)

# Generate coefficients
coe = []

# degree number
# dn=15
# mode=4
# dn = 10
# mode=3

dn = 6
mode=2

siz=1



for i in range(K + 3):
    coe.append(cp.Variable((1, dn)))

con = []  # initialize the constraints

state_pair = statepair(System)
input_pair = inputpair(U)

# Condition 1
for x_c in state_pair:
    mono = monomial(output(x_c[0]), output(x_c[1]),mode)
    act = sigma(x_c, X_R, Xf, delta)
    # print(act)
    q_ini = '0'  # initial state of DFA
    Q_ini = next(q_ini, act[0], K)

    ind = DFA.index(Q_ini)
    con.append(coe[ind] @ mono <= 0)

# Condition 2
# X_min = 0
# X_max = 10
# X_R = [X_min, X_max]
# X_R_b = [X_min - 30, X_max + 30]

outre=outregion([-0.1,10.1],[-0.3,10.3],0.1)

Q2=DFA[0:-1]

epsilon_l = 1e-3

for x_c in outre:
    for q in Q2:
        ind=Q2.index(q)
        con.append(coe[ind] @ mono >= epsilon_l)


# Condition 3
def is_val(x_c,u): #x_c=[x1,x2], u=[u1,u2]
    if (x_c[0]=='1' and u[0]=='c') or (x_c[1]=='1' and u[1]=='c'):
        return False
    elif (int(x_c[0])>=2 and int(x_c[0]) <=9 and u[0] != 'c') \
        or (int(x_c[1])>=2 and int(x_c[1]) <=9 and u[1] != 'c'):
        return False
    else:
        return True

Q3=DFA[0:-1]
state_pair = statepair(System)
input_pair = inputpair(U)

for x_c in state_pair:
    for u in input_pair:
        if is_val(x_c,u):
            x1_n=tran(x_c[0],u[0])
            x2_n=tran(x_c[1],u[1])
            x_n = [x1_n, x2_n]
            mono = monomial(output(x_c[0]), output(x_c[1]),mode)
            mono_n = monomial(output(x1_n), output(x2_n),mode)
            act = sigma(x_n, X_R, Xf, delta)

            # Qc = DFA[0:-2]  # states without 'trap' and 'F'
            for q_c in Q3: # Go through the possible current states
                q_n=next(q_c, act[0], K)
                ind=DFA.index(q_c)
                ind_n = DFA.index(q_n)
                con.append(coe[ind] @ mono >= coe[ind_n] @ mono)

objective = cp.Minimize(1)
prob = cp.Problem(objective, con)
prob.solve(solver=cp.SCS)
# prob.solve()
print("status:", prob.status)
print("optimal value", prob.value)
print("Optimal var")
# print(coe[0].value)
#
# x_t=[-5,-5]
# mono_t=monomial(output(x_t[0]), output(x_t[1]),mode)
# test=np.dot(coe[1].value,mono_t)
# print(test)

outval=[]
for i in range(K+3):
    outval.append(coe[i][0].value)
    # print(outval[i])

for i in range(K+3):
    # print('V{0}='.format(i,coe[i].value))
    print('V{0}='.format(i))
    # print('{0:.4f}, {1:.4f}x1, {2:.4f}x2, {3:.4f}x1**2, {4:.4f}x1*x2,{5:.4f}x2**2, {6:.4f}x1**3, {7:.4f}(x1**2)*x2, {8:.4f}x1*(x2**2), {9:.4f}x2**3'.format(outval[i][0], outval[i][1],outval[i][2], outval[i][3],outval[i][4],outval[i][5],outval[i][6],outval[i][7],outval[i][8],outval[i][9]))
    # print(
    #     '{0}, {1}x1, {2}x2, {3}x1**2, {4}x1*x2,{5}x2**2, {6}x1**3, {7}(x1**2)*x2, {8}x1*(x2**2), {9}x2**3'.format(
    #         outval[i][0], outval[i][1], outval[i][2], outval[i][3], outval[i][4], outval[i][5], outval[i][6],
    #         outval[i][7], outval[i][8], outval[i][9]))
    # print(
    #     '{0:.5f}, {1:.5f}x1, {2:.5f}x2, {3:.5f}x1**2, {4:.5f}x1*x2,{5:.5f}x2**2, {6:.5f}x1**3, {7:.5f}(x1**2)*x2, {8:.5f}x1*(x2**2), {9:.5f}x2**3'.format(
    #         outval[i][0], outval[i][1], outval[i][2], outval[i][3], outval[i][4], outval[i][5], outval[i][6],
    #         outval[i][7], outval[i][8], outval[i][9]))
    label=[]
    for val in outval[i]:
        if val > 0:
            label.append("+")
        elif val <0:
            label.append("-")

    if mode==3:
        print(
            '{0:.2e} {1} {2:.2e}x1 {3} {4:.2e}x2 {5} {6:.2e}x1^2 {7} {8:.2e}x1*x2 {9} {10:.2e}x2^2 {11} {12:.2e}x1^3 {13} {14:.2e}(x1^2)*x2 {15} {16:.2e}x1*(x2^2) {17} {18:.2e}x2^3'.format(
                outval[i][0], label[1], abs(outval[i][1]), label[2], abs(outval[i][2]), label[3],  abs(outval[i][3]), label[4], abs(outval[i][4]), label[5], abs(outval[i][5]), label[6], abs(outval[i][6]), label[7],
                abs(outval[i][7]), label[8], abs(outval[i][8]), label[9], abs(outval[i][9])))
    elif mode==2:
        print(
            '{0:.2e} {1} {2:.2e}x1 {3} {4:.2e}x2 {5} {6:.2e}x1^2 {7} {8:.2e}x1*x2 {9} {10:.2e}x2^2'.format(outval[i][0], label[1], abs(outval[i][1]), label[2], abs(outval[i][2]), label[3],  abs(outval[i][3]),
                label[4], abs(outval[i][4]), label[5], abs(outval[i][5])) )