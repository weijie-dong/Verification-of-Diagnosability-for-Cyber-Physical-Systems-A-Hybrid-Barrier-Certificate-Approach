import numpy as np
import random as rand
import cvxpy as cp
from z3 import *
import math


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


def stateset(K):
    q = ['0']
    for i in range(1, K + 1):
        q.append(str(i))
    q.append('F')
    q.append('trp')
    return q


'''
K-step diagnosability has K+3 states:
Delta Funcation: q_0 ... q_i ... q_K q_trap F
                  0       i        K  k+1  k+2    

Here assume K=3
'''


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


def q_0(sig):
    if (sig == 'sig_1'):
        q_n = '0'
    elif (sig == 'sig_2'):
        q_n = 'trap'
    elif (sig == 'sig_3'):
        q_n = '1'
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


'''
x1^2 , x1 * x2 , x1 * x12 , x1 * x22 , x1 , x2 ^ 2 
, x2 * x12 ,  x2 * x22 , x2 , x12 ^ 2 , x12 * x22 
    x12 , x22 ^ 2 , x22 , 1
'''


def monomial(X):
    # x1^2 , x1 * x2 , x1 * x12 , x1 * x22 , x1 , x2 ^ 2
    # , x2 * x12 ,  x2 * x22 , x2 , x12 ^ 2 , x12 * x22
    #     x12 , x22 ^ 2 , x22 , 1
    x1 = X[0]
    x2 = X[1]
    x12 = X[2]
    x22 = X[3]
    C = np.array([x1 ** 2, x1 * x2, x1 * x12, x1 * x22, x1, x2 ** 2
                     , x2 * x12, x2 * x22, x2, x12 ** 2, x12 * x22,
                  x12, x22 ** 2, x22, 1])
    return C


def output(x):
    C = np.array([0, 1])
    X = np.array(x)
    out = C.dot(X)
    return out


def P1(x, delta):
    ob = output(x[0:2]) - output(x[2:4])
    b = math.pow(ob, 2) <= math.pow(delta, 2)
    return b


def P2(x, xf):
    b = (x[0] <= xf[1]) and (x[0] >= xf[0]) \
        and (x[1] <= xf[3]) and (x[1] >= xf[2])
    return b


def P3(x, xf):
    b = (x[2] <= xf[1]) and (x[2] >= xf[0]) \
        and (x[3] <= xf[3]) and (x[3] >= xf[2])
    return b


def R(x, X_R):
    b = (x[0] <= X_R[1]) and (x[0] >= X_R[0]) \
        and (x[1] <= X_R[1]) and (x[1] >= X_R[0]) \
        and (x[2] <= X_R[1]) and (x[2] >= X_R[0]) \
        and (x[3] <= X_R[1]) and (x[3] >= X_R[0])
    return b


def sigma(x, X_R, xf, delta):  # x=[x1,x2], x1:str, Xf = [x1,x2,x3,x4]
    act = []
    if (P1(x, delta) and R(x, X_R) and (not P2(x, xf)) and (not P3(x, xf))):
        act.append('sig_1')
    if ((R(x, X_R) and (not P1(x, delta))) or (P3(x, xf) and P1(x, delta))):
        act.append('sig_2')
    if (P1(x, delta) and P2(x, xf) and (not P3(x, xf)) and R(x, X_R)):
        act.append('sig_3')
    return act


def dynamics(x0, x1, u0, u1):
    # Parameters
    tau_s = 5
    alph = 0.01
    theta = 0.04  # alpha_e
    mu = 0.145  # alpha_h
    v0 = u0
    v1 = u1
    a_00 = 1 - 2 * alph - theta - mu * v0
    a_11 = 1 - 2 * alph - theta - mu * v1
    a_01 = alph
    a_10 = alph
    A = np.array([[a_00, a_01], [a_10, a_11]])
    T_e = 10
    T_h = 50
    next_temp = np.dot(A, [x0, x1]) + [mu * T_h * v0, mu * T_h * v1] + [theta * T_e, theta * T_e]
    return (next_temp)


def dfa_state(k):
    DFA = []
    for i in range(k + 1):
        DFA.append(str(i))
    DFA.append('trap')
    DFA.append('F')
    return DFA


def com_coff(X_dis, U_dis):
    con = []  # initialize the constraints

    # generate parameter metric K+2: 1*15
    # for i in range(K + 3):
    #     exec('C{0}=cp.Variable((1,15))'.format(i))
    coe = []
    for i in range(K + 3):
        coe.append(cp.Variable((1, dn)))

    # Condition 1
    for x_c in X_dis:
        mono = monomial(x_c)
        act = sigma(x_c, X_R, Xf, delta)

        q_ini = '0'  # initial state of DFA
        Q_ini = []  # set of states enabled by x_c
        Q_ini = next(q_ini, act[0], K)

        ind = DFA.index(Q_ini)
        con.append(coe[ind] @ mono <= 0)

    # Condition 2
    for x_c in X_dis:
        mono = monomial(x_c)
        # exec('con.append(C{0} @ mono >= epsilon_l)'.format(K + 2))
        con.append(coe[K + 2] @ mono >= epsilon_l)

    # Condition 3
    for x_c in X_dis:
        for u in U_dis:
            # x = [25, 23, 23, 22.5]
            x1_n = dynamics(x_c[0], x_c[1], u[0], u[1])
            x2_n = dynamics(x_c[2], x_c[3], u[2], u[3])

            x_n = x1_n.tolist() + x2_n.tolist()  ############################################### 串联
            mono = monomial(x_c)
            mono_n = monomial(x_n)
            act = sigma(x_n, X_R, Xf, delta)
            # print(act)
            if act == []:
                # print(x_n)
                continue

            Qc = DFA[0:-2]  # states without 'trap' and 'F'
            #
            for q_c in Qc:  # Go through the possible current states
                q_n = next(q_c, act[0], K)
                ind = DFA.index(q_c)
                ind_n = DFA.index(q_n)
                con.append(coe[ind] @ mono >= coe[ind_n] @ mono)

    return con, coe


K = 5
delta = 0.5

epsilon_l = 0.001

# lower and upper bounds
X_min = 20
X_max = 30
X_R = [X_min, X_max]

U = [0, 1]

ul = U[0]
uu = U[1]

# Fault
Xf = [24, 26, X_min, X_max]
# Xf_2=[20,26]

# Sample from state and input
X_dis = gridstate(X_R, X_R, X_R, X_R, 2)  # l.append([low0, low1, low2, low3])

U_dis = gridinput(U, 0.5)

# Actions in DFA
Act = ['sig_1', 'sig_2', 'sig_3', 'sig_4', 'sig_5', 'sig_6']

# states of DFA
DFA = dfa_state(K)

dn = 15

con, coe = com_coff(X_dis, U_dis)  #

objective = cp.Minimize(1)
prob = cp.Problem(objective, con)
prob.solve(solver=cp.SCS)
print("status:", prob.status)
print("optimal value", prob.value)
print("Optimal var")

outval = []
for i in range(K + 3):
    outval.append(coe[i][0].value)
    # print(outval[i])

for i in range(K + 3):
    # print('B{0}='.format(i,coe[i].value))
    print('B{0}='.format(i))
    # print('{0:.4f}, {1:.4f}x1, {2:.4f}x2, {3:.4f}x1**2, {4:.4f}x1*x2,{5:.4f}x2**2, {6:.4f}x1**3, {7:.4f}(x1**2)*x2, {8:.4f}x1*(x2**2), {9:.4f}x2**3'.format(outval[i][0], outval[i][1],outval[i][2], outval[i][3],outval[i][4],outval[i][5],outval[i][6],outval[i][7],outval[i][8],outval[i][9]))
    print('{0:.4f} x1 ** 2, {0:.4f} x1 * x2, {0:.4f} x1 * x12, {0:.4f} x1 * x22, {0:.4f} x1,'.format(outval[i][0],
                                                                                                     outval[i][1],
                                                                                                     outval[i][2],
                                                                                                     outval[i][3],
                                                                                                     outval[i][4]),
          '{0:.4f} x2 ** 2, {0:.4f} x2 * x12, {0:.4f} x2 * x22, {0:.4f} x2, {0:.4f} x12 ** 2,'.format(outval[i][5],
                                                                                                      outval[i][6],
                                                                                                      outval[i][7],
                                                                                                      outval[i][8],
                                                                                                      outval[i][9]),
          '{0:.4f} x12 * x22, {0:.4f} x12, {0:.4f} x22 ** 2, {0:.4f} x22, {0:.4f}'.format(outval[i][10], outval[i][11],
                                                                                          outval[i][12], outval[i][13],
                                                                                          outval[i][14]), )

    print('B{0}_e='.format(i))
    print('{0:.2e} x1 ** 2, {0:.2e} x1 * x2, {0:.2e} x1 * x12, {0:.2e} x1 * x22, {0:.2e} x1,'.format(outval[i][0],
                                                                                                     outval[i][1],
                                                                                                     outval[i][2],
                                                                                                     outval[i][3],
                                                                                                     outval[i][4]),
          '{0:.2e} x2 ** 2, {0:.2e} x2 * x12, {0:.2e} x2 * x22, {0:.2e} x2, {0:.2e} x12 ** 2,'.format(outval[i][5],
                                                                                                      outval[i][6],
                                                                                                      outval[i][7],
                                                                                                      outval[i][8],
                                                                                                      outval[i][9]),
          '{0:.2e} x12 * x22, {0:.2e} x12, {0:.2e} x22 ** 2, {0:.2e} x22, {0:.2e}'.format(outval[i][10], outval[i][11],
                                                                                          outval[i][12], outval[i][13],
                                                                                          outval[i][14]), )


###################################################################################################################
# Conterexample

def is_R(x, y, x_h, y_h, X_R):
    return If(And(And(And((x <= X_R[1]), (x >= X_R[0])), \
                      And((y <= X_R[1]), (y >= X_R[0]))), \
                  And(And((x_h <= X_R[1]), (x_h >= X_R[0])), \
                      And((y_h <= X_R[1]), (y_h >= X_R[0])))), True, False)


def is_P1(x, y, x_h, y_h, delta):
    return If(y ** 2 - y_h ** 2 <= delta, True, False)


def is_P2(x, y, x_h, y_h, xf):
    return If(And(And((x <= xf[1]), (x >= xf[0])), \
                  And((y <= xf[3]), (y >= xf[2]))), True, False)


def is_P3(x, y, x_h, y_h, xf):
    return If(And(And((x_h <= xf[1]), (x_h >= xf[0])), \
                  And((y_h <= xf[3]), (y_h >= xf[2]))), True, False)


def is_sigma1(x, y, x_h, y_h, X_R, xf, delta):
    return If(And(And(is_P1(x, y, x_h, y_h, delta), is_R(x, y, x_h, y_h, X_R)),
                  And(Not(is_P2(x, y, x_h, y_h, xf)), Not(is_P3(x, y, x_h, y_h, xf)))), True, False)


def is_sigma2(x, y, x_h, y_h, X_R, xf, delta):
    return If(Or(Not(is_P1(x, y, x_h, y_h, delta)), And(is_P1(x, y, x_h, y_h, delta), is_P3(x, y, x_h, y_h, xf))), True,
              False)


def is_sigma3(x, y, x_h, y_h, X_R, xf, delta):
    return If(And(And(is_P1(x, y, x_h, y_h, delta), is_R(x, y, x_h, y_h, X_R)),
                  And(is_P2(x, y, x_h, y_h, xf), Not(is_P3(x, y, x_h, y_h, xf)))), True, False)


x, y, x_h, y_h, u1, u2, u1_h, u2_h = Reals("x y x_h y_h u1 u2 u1_h u2_h")


def B(x, y, x_h, y_h, c):
    para = np.array([x ** 2, x * y, x * x_h, x * y_h, x, y ** 2
                        , y * x_h, y * y_h, y, x_h ** 2, x_h * y_h,
                     x_h, y_h ** 2, y_h, 1])
    # exec('B{0}=np.dot(C{1}.value, para)'.format(n, n))
    b=np.dot(c, para)
    return b


# for i in range(K + 3):
#     exec('B{0}=np.dot(C{1}.value, para)'.format(i, i))

quer_i = And(True)
quer_p = And(True)
# Region of parameters
quer_i = And(quer_i, x >= X_min)
quer_i = And(quer_i, x <= X_max)
quer_i = And(quer_i, y >= X_min)
quer_i = And(quer_i, y <= X_max)
quer_i = And(quer_i, x_h >= X_min)
quer_i = And(quer_i, x_h <= X_max)
quer_i = And(quer_i, y_h >= X_min)
quer_i = And(quer_i, y_h <= X_max)

quer_i = And(quer_i, u1 >= ul)
quer_i = And(quer_i, u1 <= uu)

quer_i = And(quer_i, u2 >= ul)
quer_i = And(quer_i, u2 <= uu)

quer_i = And(quer_i, u1_h >= ul)
quer_i = And(quer_i, u1_h <= uu)

quer_i = And(quer_i, u2_h >= ul)
quer_i = And(quer_i, u2_h <= uu)

# Condition 1
sig1 = is_sigma1(x, y, x_h, y_h, X_R, Xf, delta)
sig2 = is_sigma2(x, y, x_h, y_h, X_R, Xf, delta)
sig3 = is_sigma3(x, y, x_h, y_h, X_R, Xf, delta)
# b0=B0[0]>epsilon_l
# quer_p = And(quer_p, Implies(sig1, B0[0]>=epsilon_l))
# quer_p = And(quer_p, Implies(sig2, B4[0]>=epsilon_l))
# quer_p = And(quer_p, Implies(sig3, B1[0]>=epsilon_l))

quer_p = And(quer_p, Implies(sig1, B(x, y, x_h, y_h, outval[0]) >= epsilon_l))
quer_p = And(quer_p, Implies(sig2, B(x, y, x_h, y_h, outval[K + 1]) >= epsilon_l))
quer_p = And(quer_p, Implies(sig3, B(x, y, x_h, y_h, outval[1]) >= epsilon_l))

# Condition 2
quer_p = And(quer_p, B(x, y, x_h, y_h, outval[K + 2]) <= 0)

# Condition 3
alph = 0.01
theta = 0.04
mu = 0.145
# v0 = u0
# v1 = u1

a_00 = 1 - 2 * alph - theta - mu * u1
a_11 = 1 - 2 * alph - theta - mu * u2
a_00_h = 1 - 2 * alph - theta - mu * u1_h
a_11_h = 1 - 2 * alph - theta - mu * u2_h

a_01 = alph
a_10 = alph
# A = np.array([[a_00, a_01], [a_10, a_11]])
T_e = 15
T_h = 40
# next dynamic (x,y,u1,u2) -> (x_n,y_n)   (x_h,y_h,u1_h,u2_h) -> (x_h_n,y_h_n)
x_n = a_00 * x + a_01 * y + mu * T_h * u1 + theta * T_e
y_n = a_10 * x + a_11 * y + mu * T_h * u2 + theta * T_e

x_h_n = a_00_h * x_h + a_01 * y_h + mu * T_h * u1_h + theta * T_e
y_h_n = a_10 * x_h + a_11_h * y_h + mu * T_h * u2_h + theta * T_e
# next_temp = np.dot(A, [x0, x1]) + [mu * T_h * v0, mu * T_h * v1] + [theta * T_e, theta * T_e]

sig1 = is_sigma1(x_n, y_n, x_h_n, y_h_n, X_R, Xf, delta)
sig2 = is_sigma2(x_n, y_n, x_h_n, y_h_n, X_R, Xf, delta)
sig3 = is_sigma3(x_n, y_n, x_h_n, y_h_n, X_R, Xf, delta)

# current state: 0
quer_p = And(quer_p, Implies(sig1, B(x_n, y_n, x_h_n, y_h_n, outval[0]) <= B(x, y, x_h, y_h, outval[0])))
quer_p = And(quer_p, Implies(sig2, B(x_n, y_n, x_h_n, y_h_n, outval[K + 1]) <= B(x, y, x_h, y_h, outval[0])))
quer_p = And(quer_p, Implies(sig2, B(x_n, y_n, x_h_n, y_h_n, outval[1]) <= B(x, y, x_h, y_h, outval[0])))

# current state: 1-K-1
for i in range(1, K):
    quer_p = And(quer_p, Implies(Or(sig1, sig3), B(x_n, y_n, x_h_n, y_h_n, outval[i + 1]) <= B(x, y, x_h, y_h, outval[i])))
    quer_p = And(quer_p, Implies(sig2, B(x_n, y_n, x_h_n, y_h_n, outval[K + 1]) <= B(x, y, x_h, y_h, outval[i])))

# current state: K
quer_p = And(quer_p, Implies(Or(sig1, sig3), B(x_n, y_n, x_h_n, y_h_n, outval[K + 2]) <= B(x, y, x_h, y_h, outval[K])))
quer_p = And(quer_p, Implies(sig2, B(x_n, y_n, x_h_n, y_h_n, outval[K + 1]) <= B(x, y, x_h, y_h, outval[K])))

s = Solver()

s.add(And(quer_i, Not(quer_p)))
print("Checking")
print(s.check())
print(s.model())
