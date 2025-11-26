from z3 import *
import numpy as np
import threading
import time
from func_timeout import func_set_timeout
import func_timeout

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
@func_set_timeout(60)  # 设定函数超执行时间_
def compute_process(outval):
    K = 3
    delta = 0.5
    # K = 3
    # delta = 1
    epsilon_l = 0.001
    X_min = 15
    X_max = 30
    X_R = [X_min, X_max]

    U = [0, 1]

    ul = U[0]
    uu = U[1]

    inter_max = 50

    # Fault
    Xf = [24, 26, X_min, X_max]

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
    for i in range(0, K+1):
        quer_p = And(quer_p, B(x, y, x_h, y_h, outval[i]) <= 0)

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
    quer_p = Or(quer_p, Implies(sig1, B(x_n, y_n, x_h_n, y_h_n, outval[0]) <= B(x, y, x_h, y_h, outval[0])))
    quer_p = Or(quer_p, Implies(sig2, B(x_n, y_n, x_h_n, y_h_n, outval[K + 1]) <= B(x, y, x_h, y_h, outval[0])))
    quer_p = Or(quer_p, Implies(sig2, B(x_n, y_n, x_h_n, y_h_n, outval[1]) <= B(x, y, x_h, y_h, outval[0])))

    # current state: 1-K-1
    for i in range(1, K):
        quer_p = Or(quer_p, Implies(Or(sig1, sig3), B(x_n, y_n, x_h_n, y_h_n, outval[i + 1]) <= B(x, y, x_h, y_h, outval[i])))
        quer_p = Or(quer_p, Implies(sig2, B(x_n, y_n, x_h_n, y_h_n, outval[K + 1]) <= B(x, y, x_h, y_h, outval[i])))

    # current state: K
    quer_p = Or(quer_p, Implies(Or(sig1, sig3), B(x_n, y_n, x_h_n, y_h_n, outval[K + 2]) <= B(x, y, x_h, y_h, outval[K])))
    quer_p = Or(quer_p, Implies(sig2, B(x_n, y_n, x_h_n, y_h_n, outval[K + 1]) <= B(x, y, x_h, y_h, outval[K])))

    s = Solver()

    s.add(And(quer_i, Not(quer_p)))
    print("Checking: \n")
    o = s.check()
    model = s.model()
    if o == sat:
        print("反例结果：", model) #x y x_h y_h u1 u2 u1_h u2_h
        result_x = float(model.eval(x).as_fraction())
        result_y = float(model.eval(y).as_fraction())
        result_x_h = float(model.eval(x_h).as_fraction())
        result_y_h = float(model.eval(y_h).as_fraction())
        result_u1 = float(model.eval(u1).as_fraction())
        result_u2 = float(model.eval(u2).as_fraction())
        result_u1_h = float(model.eval(u1_h).as_fraction())
        result_u2_h = float(model.eval(u2_h).as_fraction())

        return [result_x,result_y,result_x_h,result_y_h], [result_u1,result_u2,result_u1_h,result_u2_h], 1
    else:
        # print("Counterexample cannot be found!")
        return [], [], 0

##############################################################
# from func_timeout import func_set_timeout
# import func_timeout
# @func_set_timeout(10)#设定函数超执行时间_
# def compute_counter(outval):
#     try:
#         return compute_process(outval)
#     # 若调用函数超时自动走异常(可在异常中写超时逻辑处理)
#     except func_timeout.exceptions.FunctionTimedOut:
#         print("Counterexample cannot be found!")
#         return [], [], 0


##############################################################
import time
import eventlet  #导入eventlet这个模块
def compute_counter(outval):
    eventlet.monkey_patch()   #必须加这条代码
    with eventlet.Timeout(10,False):   #设置超时时间为2秒
       counter_s, counter_u, flag=compute_process(outval)
       return counter_s, counter_u, flag
    print("Counterexample cannot be found!")
    return [], [], 0