'''张子扬202104061202'''
import numpy as np
import matplotlib.pyplot as plt
import math
def lagrange(arr_x,arr_y,_x):
    len_x = len(arr_x)
    len_y = len(arr_y)
    l = []
    result = 0
    for i in range(len_x):
        denominator = 1
        molecular = 1
        for j in range(len_y):
            if i != j:
                denominator *= arr_x[i] - arr_x[j]
                molecular *= _x - arr_x[j]
        l.append(molecular / denominator) 
        result += l[i] * arr_y[i]
    return result

def original(x):
    result = math.exp(x)
    return result


original_x = np.arange(1,10)
original_y = []
for i in range(len(original_x)):
    original_y.append(original(original_x[i]))

x = np.arange(1,10)
y = []

for i in range(len(x)):
    y.append(lagrange(original_x,original_y,x[i]))
e = (original(11)-lagrange(original_x,original_y,10))/original(11)
print(e)
plt.plot(original_x, original_y, label='f(x) = exp(x)', color='black')
plt.plot(x,y,label='p(x)', color='blue')
plt.scatter(original_x,original_y,label='The interpolation points', color='red')
plt.title('Lagrange interpolation')
plt.legend(loc='upper left')
plt.xlabel("x")
plt.ylabel("y")
plt.show()