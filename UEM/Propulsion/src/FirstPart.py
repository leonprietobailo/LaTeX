import sympy as sp
import math

ma = 66.7        # kg/s
pi_i = 0.96     # 1
pi_c = 18       # 1
n_c = 0.99      # 1
Tt4 = 1456      # K
n_cc = 0.96     # 1
pi_cc = 0.96    # 1
n_t = 0.91      # 1
L = 42.8e6      # J/kg

#? Suposiciones
M0 = 0
gamma = 1.4
R = 287

# Etapa 0
# Presion y atmosfera estandar a nivel de mar.
T0 = 288
P0 = 101325

Tt0 = T0 * (1 + (gamma-1)/2 * M0 ** 2)
Pt0 = P0 * (1 + (gamma-1)/2 * M0 ** 2) ** (gamma / (gamma - 1))

# Etapa 2
Pt2 = Pt0 * pi_i
Tt2 = Tt0 * pi_i ** ((gamma-1)/gamma)

# Etapa 3
Pt3 = Pt2 * pi_c
Tt3 = sp.symbols("Tt3")
equality = sp.Eq(n_c, (pi_c ** ((gamma-1)/gamma)-1) / (Tt3/Tt2-1))
Tt3 = sp.solve(equality, Tt3)[0]

# Etapa 4
Pt4 = pi_cc * Pt3
Cp = R * gamma / (gamma - 1)
Tt4p = sp.symbols("Tt4p")
equality = sp.Eq(n_cc, (Tt4-Tt3) / (Tt4p - Tt3))
Tt4p = sp.solve(equality, Tt4p)[0]
f = (Tt4p/Tt3 - 1)/(L/Cp/Tt3 - Tt4p/Tt3)
mp = ma * (1 + f)


# Etapa 5

Tt5 = sp.symbols("Tt5")
equality = sp.Eq(ma * Cp * (Tt3 - Tt2), mp * Cp * (Tt4 - Tt5))
Tt5 = sp.solve(equality, Tt5)[0]

Pt5 = sp.symbols("Pt5")
equality = sp.Eq(n_t, (1-Tt5/Tt4)/(1-(Pt5/Pt4)**((gamma-1)/gamma)))
Pt5 = sp.solve(equality, Pt5)[0]


# Etapa 7
Pt7 = Pt5
Tt7 = Tt5

P_asterisco = (1 + (gamma - 1) / 2) ** (gamma/(gamma - 1))

T7 = None
P7 = None
M7 = None
if(Pt7 / P0 > P_asterisco): # Condiciones criticas
    T7 = Tt7 / (1 + (gamma - 1)/2)
    P7 = Pt7 / (1 + (gamma - 1)/2) ** (gamma/(gamma-1))
    M7 = 1

a7 = math.sqrt(gamma*R*T7)
v7 = M7 * a7
rho_7 = P7 / R / T7
A7 = mp / rho_7 / v7

T = mp * v7 + (P7 - P0) * A7
P = ma*Cp*(Tt3-Tt2)

mf = ma * f
consumo_sp = mf / T * 1e6
rendimiento_termico = P / mf / L


print("Pt0", Pt0)
print("Tt0", Tt0)

print("Pt2", Pt2)
print("Tt2", Tt2)

print("Pt3", Pt3)
print("Tt3", Tt3)

print("Pt4", Pt4)
print("Tt4", Tt4)

print("Pt5", Pt5)
print("Tt5", Tt5)

print("Pt7", Pt7)
print("Tt7", Tt7)

print("P7", P7)
print("T7", T7)

print("Pt7/P0", Pt7 / P0)
print("V7",v7)
print("Rho_7", rho_7)
print("A7",A7)
print("T",T)
print("TSFC", consumo_sp)

print("P", P)
print("f", f)
print("Tt4'", Tt4p)
print("Rendimiento Termico", rendimiento_termico)