import math
from matplotlib import pyplot as plt
import numpy

"""Dados pré determinados"""
L1 = 7
L2 = 2
L3 = 6
L4 = 4
angulo2 = 90
angulo2_radiano = math.radians(angulo2)

"""Valores iniciais para Theta3 e Theta4 """
angulo3 = 20
angulo3_radiano = math.radians(angulo3)
angulo4 = 100
angulo4_radiano = math.radians(angulo4)

erro = 10**(-13)
X = 1

while (numpy.linalg.norm(X) > erro):
    f1 = L2*math.cos(angulo2_radiano) + L3*math.cos(angulo3_radiano) - L1 - L4*math.cos(angulo4_radiano)
    f2 = L2*(math.sin(angulo2_radiano)) + L3*(math.sin(angulo3_radiano))-L4*math.sin(angulo4_radiano)
    J=[[-L3*math.sin(angulo3_radiano), L4*math.sin(angulo4_radiano)], [L3*math.cos(angulo3_radiano), -L4*math.cos(angulo4_radiano)]]
    b=[[-f1], [-f2]]
    X=numpy.linalg.inv(J)*b
    angulo3_radiano = angulo3_radiano + X[0][0]
    angulo4_radiano = angulo4_radiano + X[1][0]

angulo3_final = math.degrees(angulo3_radiano)
angulo4_final = math.degrees(angulo4_radiano)
print(f"Para um theta2 de entrada igual a {angulo2}, obtêm-se theta3 igual a {angulo3_final} e theta4 igual a {angulo4_final}")



