#@Codigo produzido por Thiago Jose da Silva
#@Data: 15/12/2021

import numpy
import math
from matplotlib import pyplot as plt

def equacao_biela_manivela(R2, R3, angulo):
    """Retorna o valor de R1"""
    return (R2*math.cos(angulo)) + (R3*(numpy.sqrt(1-((-R2/R3)*math.sin(angulo))**2)))

R2 = 3
R3 = 4
angulo_theta_1 = list(range(361))
angulo_radiano = []

for angulo in angulo_theta_1:
    angulo2 = math.radians(angulo)
    angulo_radiano.append(angulo2)

movimento = []

for theta in angulo_radiano:
    movimento_insta = equacao_biela_manivela(R2, R3, theta)
    movimento.append(movimento_insta)

plt.figure(1)
plt.title('Movimento linear do pistão')
plt.plot(angulo_theta_1,movimento,'b',label='$R_1$')
plt.grid()
plt.ylabel('$R_1$')
plt.xlabel('Ângulo')
plt.legend()
plt.ylim(0,8)
plt.xlim(0,max(angulo_theta_1))
plt.show()

