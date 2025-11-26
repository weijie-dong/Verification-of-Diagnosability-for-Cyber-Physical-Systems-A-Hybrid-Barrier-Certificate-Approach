import cvxpy as cp
import numpy as np
import math


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


def com_coff(X_dis, U_dis,DFA,K):
    con = []  # initialize the constraints
    #parameters
    K = 5
    delta = 0.5
    # K = 3
    # delta = 1

    epsilon_l = 0.001

    # lower and upper bounds
    X_min = 20
    X_max = 30
    X_R = [X_min, X_max]

    U = [0, 1]

    ul = U[0]
    uu = U[1]

    inter_max = 50

    # Fault
    Xf = [24, 26, X_min, X_max]
    dn = 15
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


# def coefficients_B(X_dis,U_dis,DFA):
#     con, coe = com_coff(X_dis, U_dis,DFA)  #
#
#     objective = cp.Minimize(1)
#     prob = cp.Problem(objective, con)
#     prob.solve(solver=cp.SCS)
#     print("status:", prob.status)

def print_coff(prob,coe,K):
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
              '{0:.4f} x12 * x22, {0:.4f} x12, {0:.4f} x22 ** 2, {0:.4f} x22, {0:.4f}'.format(outval[i][10],
                                                                                              outval[i][11],
                                                                                              outval[i][12],
                                                                                              outval[i][13],
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
              '{0:.2e} x12 * x22, {0:.2e} x12, {0:.2e} x22 ** 2, {0:.2e} x22, {0:.2e}'.format(outval[i][10],
                                                                                              outval[i][11],
                                                                                              outval[i][12],
                                                                                              outval[i][13],
                                                                                              outval[i][14]), )

    return outval

def print_outval(outval,K):
    for i in range(K + 3):
        print('B{0}='.format(i))
        label = []
        for val in outval[i]:
            if val > 0:
                label.append("+")
            elif val < 0:
                label.append("-")

        print('{0:.4f} x1 ^ 2 {1} {2:.4f} x1 * x2 {3} {4:.4f} x1 * x12 {5} {6:.4f} x1 * x22 {7} {8:.4f} x1 {9}'.format(outval[i][0], label[1],
                                                          abs(outval[i][1]), label[2],
                                                          abs(outval[i][2]), label[3],
                                                          abs(outval[i][3]), label[4],
                                                          abs(outval[i][4]),label[5]),
              '{0:.4f} x2 ^ 2 {1} {2:.4f} x2 * x12 {3} {4:.4f} x2 * x22 {5} {6:.4f} x2 {7} {8:.4f} x12 ^ 2 {9}'.format(outval[i][5], label[6],
                                                          abs(outval[i][6]),label[7],
                                                          abs(outval[i][7]),label[8],
                                                          abs(outval[i][8]),label[9],
                                                          abs(outval[i][9]),label[10]) ,
              '{0:.4f} x12 * x22 {1} {2:.4f} x12 {3} {4:.4f} x22 ^ 2 {5} {6:.4f} x22 {7} {8:.4f}'.format(outval[i][10], label[11],
                                                          abs(outval[i][11]),label[12],
                                                          abs(outval[i][12]),label[13],
                                                          abs(outval[i][13]),label[14],
                                                          abs(outval[i][14]) ), )

