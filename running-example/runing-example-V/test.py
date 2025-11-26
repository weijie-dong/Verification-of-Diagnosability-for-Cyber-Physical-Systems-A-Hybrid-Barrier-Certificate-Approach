out = [1, 2, -3, 4]
i = 0
print(out[0], end="")
i=i+1
while i < 4:
    if out[i] >0:
        print(" +", out[i], end="")
    elif out[i] <0:
        print(" -", abs(out[i]), end="")
    i=i+1

print("\n")
val=[]
val.append("+")
print(val[0])