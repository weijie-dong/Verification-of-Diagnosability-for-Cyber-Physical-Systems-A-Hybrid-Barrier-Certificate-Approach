import numpy as np
import random as rand
from z3 import *
import cvxpy as cp
import math


def tran(q, u): # q=x , u=u
    ind=int(q)
    if ind==1:
        if u=='a':
            temp='2'
        elif u=='b':
            temp = '6'
        else:
            print('input error')
    elif ind< 5 and ind>=2:
        if u=='c':
            temp=str(ind+1)
        else:
            print('input error')
    elif ind<9 and ind>=6:
        if u=='c':
            temp=str(ind+1)
        else:
            print('input error')
    elif ind==5 or ind==9:
        temp=q
    return temp



# degree: 0:3
# def monomial(par,par_h):
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
#     x1=par
#     x2 = par_h
#     mon=[1, x1, x2, x1**2, x1*x2, x2**2, x1**3, (x1**2)*x2, x1*(x2**2), x2**3] # 1,10
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


def dfa_state(k):
    DFA=[]
    for i in range(k+1):
        DFA.append(str(i))
    DFA.append('trap')
    DFA.append('F')
    return DFA
# DFA = ['0', '1', '2', '3', 'trp', 'F']

def statepair(q): # q is the state set
    l=[]
    for i in q:
        for j in q:
            l.append([i,j])
    return l

def inputpair(U): # q is the state set
    l=[]
    for i in U:
        for j in U:
            l.append([i,j])
    return l

def P1(x, delta): # x=[x1,x2], x1:str
    ob = output(x[0]) - output(x[1])
    b = math.pow(ob, 2) <= math.pow(delta, 2)
    return b


def P2(x, xf): # x=[x1,x2], x1:str, Xf = [0, 0.8]
    b = (output(x[0]) <= xf[1]) and (output(x[0]) >= xf[0])
    return b


def P3(x, xf): # x=[x1,x2], x1:str, Xf = [0, 0.8]
    b = (output(x[1]) <= xf[1]) and (output(x[1]) >= xf[0])
    return b


def R(x, X_R): # x=[x1,x2], x1:str, X_R = [X_min, X_max]
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

def next(q, sig,K):
    if (q == '0'):
        q_n = q_0(sig)
    elif (q == 'trap'):
        q_n = q_t(sig)
    elif (q == 'F'):
        q_n = q_F(sig)
    elif (q2int(q, K) <= K and q2int(q, K) >= 1):
        q_n = q_i(q, sig, K)
    return (q_n)

K = 3
delta = 1

# Parameters

# epsilon_1=0;
epsilon_l = 0.001

X_min = 0
X_max = 10
X_R = [X_min, X_max]

U = ['a', 'b','c'];
# U2_all = (u2-U_min)*(U_max-u2);

# Fault: less than 0.8
Xf = [1, 1.5]

def output(q):
    out=[0,
           1.2, 3.2, 5.2, 7.2,
           2.2, 4.2, 6.2, 9]
    ind=int(q)-1
    return out[ind]

System=['1', '2', '3', '4', '5','6', '7', '8', '9']

# Actions in DFA, new version
Act = ['sig_1', 'sig_2', 'sig_3']

DFA=dfa_state(K)


#Generate coefficients
coe=[]

#degree number
# dn=15
dn=10

mode=3

for i in range(K + 3):
   coe.append(cp.Variable((1,dn)))
   # print(coe)
   # print(coe[i])
   # print(coe[i][0])

con = [] #initialize the constraints

state_pair=statepair(System)
input_pair=inputpair(U)

# Condition 1

# print(state_pair)
# act = sigma([1,9], X_R, Xf, delta)

for x_c in state_pair:
    mono=monomial(output(x_c[0]),output(x_c[1]),mode)
    act = sigma(x_c, X_R, Xf, delta)
    print(act)
    q_ini = '0'  # initial state of DFA
    Q_ini = next(q_ini, act[0],K)

    ind=DFA.index(Q_ini)
    con.append(coe[ind] @ mono <= 0)

# Condition 2
for x_c in state_pair:
    mono = monomial(output(x_c[0]), output(x_c[1]),mode)
    con.append(coe[K+2] @ mono >= epsilon_l)

# Condition 3
def is_val(x_c,u): #x_c=[x1,x2], u=[u1,u2]
    if (x_c[0]=='1' and u[0]=='c') or (x_c[1]=='1' and u[1]=='c'):
        return False
    elif (int(x_c[0])>=2 and int(x_c[0]) <=9 and u[0] != 'c') \
        or (int(x_c[1])>=2 and int(x_c[1]) <=9 and u[1] != 'c'):
        return False
    else:
        return True

for x_c in state_pair:
    for u in input_pair:
        if is_val(x_c,u):
            x1_n=tran(x_c[0],u[0])
            x2_n=tran(x_c[1],u[1])
            x_n=[x1_n,x2_n]
            mono = monomial(output(x_c[0]), output(x_c[1]),mode)
            mono_n = monomial(output(x1_n), output(x2_n),mode)
            act = sigma(x_n, X_R, Xf, delta)

            Qc = DFA[0:-2]  # states without 'trap' and 'F'
            for q_c in Qc: # Go through the possible current states
                q_n=next(q_c, act[0], K)
                ind=DFA.index(q_c)
                ind_n = DFA.index(q_n)
                con.append(coe[ind] @ mono >= coe[ind_n] @ mono)

objective = cp.Minimize(1)
prob = cp.Problem(objective, con)
prob.solve(solver=cp.ECOS)
# prob.solve()
print("status:", prob.status)
print("optimal value", prob.value)
print("Optimal var")

outval=[]
for i in range(K+3):
    outval.append(coe[i][0].value)
    # print(outval[i])

for i in range(K+3):
    # print('B{0}='.format(i,coe[i].value))
    print('B{0}='.format(i))
    label = []
    for val in outval[i]:
        if val > 0:
            label.append("+")
        elif val < 0:
            label.append("-")

    print(
          '{0:.2f} {1} {2:.2f}x1 {3} {4:.2f}x2 {5} {6:.2f}x1^2 {7} {8:.2f}x1*x2 {9} {10:.2f}x2^2 {11} {12:.2f}x1^3 {13} {14:.2f}(x1^2)*x2 {15} {16:.2f}x1*(x2^2) {17} {18:.2f}x2^3'.format(
          outval[i][0], label[1], abs(outval[i][1]), label[2], abs(outval[i][2]), label[3], abs(outval[i][3]), label[4], abs(outval[i][4]), label[5], abs(outval[i][5]), label[6],
          abs(outval[i][6]), label[7], abs(outval[i][7]), label[8], abs(outval[i][8]), label[9], abs(outval[i][9])))


    # print('{0:.4f}, {1:.4f}x1, {2:.4f}x2, {3:.4f}x1**2, {4:.4f}x1*x2,{5:.4f}x2**2, {6:.4f}x1**3, {7:.4f}(x1**2)*x2, {8:.4f}x1*(x2**2), {9:.4f}x2**3'.format(outval[i][0], outval[i][1],outval[i][2], outval[i][3],outval[i][4],outval[i][5],outval[i][6],outval[i][7],outval[i][8],outval[i][9]))
    #
    # print(
    #     '{0:.2e}, {1:.2e}x1, {2:.2e}x2, {3:.2e}x1**2, {4:.2e}x1*x2,{5:.2e}x2**2, {6:.2e}x1**3, {7:.2e}(x1**2)*x2, {8:.2e}x1*(x2**2), {9:.2e}x2**3'.format(
    #         outval[i][0], outval[i][1], outval[i][2], outval[i][3], outval[i][4], outval[i][5], outval[i][6],
    #         outval[i][7], outval[i][8], outval[i][9]))

# for x_c in X_dis:
#     mono = monomial(x_c)
#     act = sigma(x_c, X_R, Xf, delta)
#
#     q_ini = '0'  # initial state of DFA
#     Q_ini = []  # set of states enabled by x_c
#     if (Act[0] in act):
#         Q_ini.append(next(q_ini, Act[0]))
#     elif (Act[1] in act):
#         Q_ini.append(next(q_ini, Act[1]))
#     elif (Act[2] in act):
#         Q_ini.append(next(q_ini, Act[2]))
#
#     if (DFA[0] in Q_ini):
#         exec('con.append(C0 @ mono <= 0)')
#     if (DFA[1] in Q_ini):
#         exec('con.append(C1 @ mono <= 0)')
#     if (DFA[K + 1] in Q_ini):
#         exec('con.append(C{0} @ mono <= 0)'.format(K + 1))
