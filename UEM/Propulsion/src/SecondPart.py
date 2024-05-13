import sympy as sp
import math

M = 0.5
b_1 = math.radians(51)  # math.radians(51)
# b_2 = math.radians(51)
h_c = 1.5
w2_w1 = 0.7
v = 450
nu = 1.73e-5
diff_factor = 0.45
a_3 = 0
pi_c = 18
m = 66.67
P1 = 101325
Utip = 250

# Atmosfera estandar
Tt1 = 288
gamma = 1.4
R = 287
Cp = 287 * 1.4 / (1.4 - 1)
# Ángulos relativos y absolutos de la corriente de aire, deflexión de la corriente en rotor y estátor (en línea media)
a_1 = 0
T1 = Tt1 / (1 + (gamma-1)/2 * M**2)
Ca_1 = M * math.sqrt(gamma*R*T1)


W1 = Ca_1 / math.cos(b_1)
U = Ca_1 / math.tan(b_1)

W2 = w2_w1 * W1
Ca_2 = Ca_1
Ca_3 = Ca_1
b_2 = math.acos(Ca_2 / W2)
W_u2 = W2 * math.sin(b_2)
V_u2 = U - W_u2
C2 = math.sqrt(V_u2**2 + Ca_2**2)
a2 = math.atan(V_u2 / Ca_2)
d12 = b_1 - b_2

T02rel = T1 + W1**2 / 2 / Cp
T2 = T02rel - W2**2 / 2 / Cp
Tt3 = T2 + C2**2 / 2 / Cp

pt3pt1 = (Tt3 / Tt1)**(gamma / (gamma - 1))
n_esc = pi_c / pt3pt1

rho1 = P1 / R / T1

A = m / rho1 / Ca_1

R, Re, Ri= sp.symbols('R Re Ri')

eq1 = sp.Eq(A, math.pi * (Re**2 - Ri**2))
eq2 = sp.Eq(R, (Re + Ri) / 2)
eq3 = sp.Eq(R / Re, U / Utip)
(R, Re, Ri) = sp.solve([eq1, eq2, eq3], (R, Re, Ri))[1]

w = U / R


c = (Re - Ri) / h_c
sigma = abs(W_u2 - U) / 2 / W1 / (diff_factor - 1 + W2 / W1)
S = c / sigma
N = math.ceil(2 * math.pi * R / S)

Reynolds = Ca_1 * c / nu
T3 = Tt3 - Ca_3**2 / 2 / Cp
lam = (T2 - T1) / (T3 - T1)



print("alfa_1", math.degrees(a_1))
print("alfa_2", math.degrees(a2))
print("alfa_3", math.degrees(a_3))
print("beta_1", math.degrees(b_1))
print("beta_2", math.degrees(b_2))
print("beta_3", math.degrees(b_1))
print("delta_12", math.degrees(d12))
print("n_esc", n_esc)
print("salto de presion", pt3pt1)
print("Area de la tobera", A)
print("Ri", Ri)
print("Re", Re) 
print("R", R)
print("w", w)
print("c", c)
print("N", N)
print("reynolds", Reynolds)
print("lambda", lam)