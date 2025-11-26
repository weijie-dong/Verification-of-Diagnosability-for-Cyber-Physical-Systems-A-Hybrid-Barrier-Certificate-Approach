import numpy as np
import random as rand
import cvxpy as cp
from z3 import *
import math

# x=cp.Variable()
# y=cp.Variable()
# x_h=cp.Variable()
# y_h=cp.Variable()
# u1=cp.Variable()
# u2=cp.Variable()
# u1_h=cp.Variable()
# u2_h=cp.Variable()
#
# mono=np.array([x ** 2, x * y, x * x_h, x * y_h, x, y ** 2
#                      , y * x_h, y * y_h, y, x_h ** 2, x_h * y_h,
#                   x_h, y_h ** 2, y_h, 1])
# print(mono.shape)
# A=np.random.randn(1,15)
# print(A)
# B=A@ np.transpose(mono)
# print(B.shape)
# print(B)
# C=A@ mono
# print(C)
# # numpy_array = np.array(matrix)
# # A = np.random.randn(3, 15)
# con=[C[0]]
# # print(con)
#
# obj = cp.Maximize(con[0])
# prob = cp.Problem(obj)
# prob.solve()  # Returns the optimal value.
# print("status:", prob.status)
# print("optimal value", prob.value)
# print("optimal var", x.value, y.value)

import time
import eventlet#导入eventlet这个模块
eventlet.monkey_patch()#必须加这条代码
with eventlet.Timeout(2,False):#设置超时时间为2秒
  time.sleep(4)
  print('没有跳过这条输出')
print('跳过了输出')