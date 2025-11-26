from z3 import *
import numpy as np
import cvxpy as cp


# x, y = Reals('x y')
# s = Solver()
# s.add(x - y == 3, 3 * x - 8 * y == 4)
#
# o=s.check()
# model=s.model()
# if o == sat:
#   print("结果：", s.model())
#
# result_x=float(model.eval(x).as_fraction())
# print("result_x: ", result_x)
# print("type of result_x: ", type(result_x))

print(cp.INFEASIBLE)