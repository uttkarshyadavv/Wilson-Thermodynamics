# Doing It for Acetone(1)/water(2) Mixture
# Wilson Parameter== delta
import math
import matplotlib.pyplot as plt
V1 = 74.05 
V2 = 18.07 
A12 = 291.27
A21 = 1448.01
R = 8.314
T = 333.15

x1_val=[]
P_val=[]
gamma1_val=[]
gamma2_val=[]
x2_val=[]
y1_val=[]
y2_val=[]

delta12 = (V2 / V1) * math.exp(-A12 / (R * T))
delta21 = (V1 / V2) * math.exp(-A21 / (R * T))
print("delta12, delta21 =", delta12, delta21)

# Antoine coefficients (kept in the form you originally used)
A1, B1, C1 = 14.3145, 2756.22, -45.090
A2, B2, C2 = 16.3872, 3885.70, -42.980

# Using the same Antoine expression you wrote (math.exp). 
# If your source uses log10 form, swap to 10 ** (...) later.
P1sat = math.exp(A1 - (B1 / (T + C1)))
P2sat = math.exp(A2 - (B2 / (T + C2)))

# loop in 0.01 steps (keeps your original stepping and structure)
x = 0.0
step = 0.01
while x <= 1.0:
    x1 = round(x, 6)
    # avoid exact 0 or 1 to prevent division-by-zero in Wilson denominators
    if x1 == 0.0 or x1 == 1.0:
        x += step
        continue

    x2 = 1.0 - x1

    # Wilson denominators
    denom1 = x1 + delta12 * x2
    denom2 = x2 + delta21 * x1

    # ln(gamma) form first, then exponentiate (keeps your math style)
    ln_gamma1 = -math.log(denom1) + x2 * ((delta12 / denom1) - (delta21 / denom2))
    ln_gamma2 = -math.log(denom2) + x1 * ((delta21 / denom2) - (delta12 / denom1))

    gamma1 = math.exp(ln_gamma1)
    gamma2 = math.exp(ln_gamma2)

    P = x1 * gamma1 * P1sat + x2 * gamma2 * P2sat

    x1_val.append(x1)
    P_val.append(P)
    gamma1_val.append(gamma1)
    gamma2_val.append(gamma2)
    x2_val.append(x2)

    x += step

# compute vapor mole fractions (y1, y2) using stored lists
for i in range(len(x1_val)):
    y1 = (x1_val[i] * P1sat * gamma1_val[i]) / P_val[i]
    y2 = (x2_val[i] * P2sat * gamma2_val[i]) / P_val[i]
    y1_val.append(y1)
    y2_val.append(y2)

# quick printout of first few results
for i in range(min(10, len(x1_val))):
    print(f"x1={x1_val[i]:.3f}, y1={y1_val[i]:.3f}, P={P_val[i]:.6f}, gamma1={gamma1_val[i]:.6f}")

# --- 1) x-y diagram (liquid vs vapor) ---
plt.figure(figsize=(6,5))
plt.plot(x1_val, y1_val, label='y1 vs x1')
plt.plot([0,1],[0,1],'--', label='x=y', linewidth=0.8)
plt.xlabel('x1 (liquid acetone mole fraction)')
plt.ylabel('y1 (vapor acetone mole fraction)')
plt.title('x-y Diagram (Wilson model)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# --- 2) Pressure vs x1 (bubble curve) ---
plt.figure(figsize=(6,4))
plt.plot(x1_val, P_val)
plt.xlabel('x1 (liquid acetone mole fraction)')
plt.ylabel('Bubble pressure (same units as P_sat)')
plt.title('P_bubble vs x1')
plt.grid(True)
plt.tight_layout()
plt.show()

# --- 3) Activity coefficients vs x1 ---
plt.figure(figsize=(6,4))
plt.plot(x1_val, gamma1_val, label='gamma1 (acetone)')
plt.plot(x1_val, gamma2_val, label='gamma2 (water)')
plt.xlabel('x1 (liquid acetone mole fraction)')
plt.ylabel('Activity coefficient')
plt.title('gamma vs x1')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
