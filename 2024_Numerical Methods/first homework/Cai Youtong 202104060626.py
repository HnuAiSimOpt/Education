#蔡佑桐 车辆2101 202104060626
def lagrange_interpolation(x, y, xi):
    result = 0
    for i in range(len(y)):
        term = y[i]
        for j in range(len(x)):
            if j != i:
                term = term * (xi - x[j]) / (x[i] - x[j])
        result += term
    return result# -*- coding: utf-8 -*-
x=[1,2,3,5]
y=[1,4,9,25]
xi=2.5
result=lagrange_interpolation(x,y,xi)
print(result)
"""
Spyder Editor

This temporary script file is located here:
C:\Users\Administrator\.spyder2\.temp.py
"""

