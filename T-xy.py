# We have a system of Acetone(1)/water(2)
# we have to make a plot between T-xy
import math
import matplotlib.pyplot as plt

P=101.33
R= 8.314
V1,V2= 74.05, 18.07
a12,a21= 291.27, 1448.01

A1,B1,C1= 14.3145, 2756.22, -45.090
A2,B2,C2= 16.3872, 3885.7, -42.98
T1sat = (B1/(A1 - math.log(P))) - C1
T2sat = (B2/(A2 - math.log(P))) - C2

T = (T1sat + T2sat)/2

P1sat = math.exp(A1 - (B1/(T + C1)))
P2sat = math.exp(A2 - (B2/(T + C2)))

lambda12 = (V2/V1)*math.exp(-a12/(R*T))
lambda21 = (V1/V2)*math.exp(-a21/(R*T))

x1_val = []
y1_val = []
T_val  = []
x1 = 0.0
while x1 <= 1.0:
    if x1 == 0 or x1 == 1:
        x1 += 0.01
        continue
    x2 = 1 - x1
    gamma1 = math.exp(-math.log(x1 + lambda12*x2) 
                      + x2 * ((lambda12/(x1+lambda12*x2)) - (lambda21/(x2+lambda21*x1))))

    gamma2 = math.exp(-math.log(x2 + lambda21*x1) 
                      + x1 * ((lambda21/(x2+lambda21*x1)) - (lambda12/(x1+lambda12*x2))))

    y1 = (x1 * gamma1 * P1sat) / P
    y2 = (x2 * gamma2 * P2sat) / P

    x1_val.append(x1)
    y1_val.append(y1)
    T_val.append(T)
    if P < x1*gamma1*P1sat+ x2*gamma2*P2sat:
        T = (T + T1sat) / 2
    elif P > x1*gamma1*P1sat+ x2*gamma2*P2sat:
        T = (T + T2sat) / 2
    lambda12 = (V2/V1)*math.exp(-a12/(R*T))
    lambda21 = (V1/V2)*math.exp(-a21/(R*T))
    P1sat = math.exp(A1 - (B1/(T + C1)))
    P2sat = math.exp(A2 - (B2/(T + C2)))

    x1 += 0.01
plt.plot(x1_val, T_val, label="T vs x1")
plt.plot(y1_val, T_val, label="T vs y1")

plt.xlabel("x1 / y1")
plt.ylabel("Temperature (K)")
plt.title("T - x,y Diagram (Your Wilson Model)")
plt.grid(True)
plt.legend()
plt.show()
