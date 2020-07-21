import sympy as syp
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

V_a = 100
V_b_space = np.linspace(0,130,500)
C_a = 0.01
C_b = 0.01
K_a = 1.76*10**-5
K_w = 1.0*10**-14
y = syp.symbols("y")
ph_space = np.array([])
i = 0
for V_b in V_b_space:
    i += 1
    if V_b <= V_a:
        print(i,"bef")
        equ = syp.Eq(y * (C_b*V_b + y) / (C_a*V_a - C_b*V_b - y),K_a*(V_a + V_b))
        h3o = syp.solve(equ,y)
        ph = -np.log10(float(h3o[1])/(V_a + V_b))
        ph_space = np.append(ph_space,[ph])

    if V_b > V_a:
        print(i,"aft")
        equ = syp.Eq(y * (C_b*V_b - C_b*V_a + y) / (C_b*V_a - y),K_w / K_a * (V_a + V_b))
        h3o = syp.solve(equ, y)
        ph = -np.log10(K_w / (C_b*V_b -C_b*V_a +float(h3o[1])) * (V_a + V_b) )
        ph_space = np.append(ph_space, [ph])


plt.plot(V_b_space,ph_space)
plt.show()