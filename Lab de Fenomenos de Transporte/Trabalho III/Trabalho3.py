#@Laboratório de Fenômenos de Transporte
#@Autores: Thiago José da Silva, Luiza Gomes de Castro e Sá e Victor Alves Morais
#@Data: 11/02/2022

import numpy as np                          #Importando a biblioteca numpy
from matplotlib import pyplot as plt        #Importando a biblioteca matplotlib responsável por plotar gráficos
import statistics as sta                    #Importando a biblioteca statistic
import math

valores_vazao = list(range(0, 40, 1))
valores_Hman = []
for Q in valores_vazao:
    eq_altura_man = 48.6 - 0.001740762*Q**2
    valores_Hman.append(eq_altura_man)

"""
Plotando a curva característica do sistema
"""
plt.figure(1)
plt.title("Curva Característica do Sistema")
plt.plot(valores_vazao,valores_Hman,'b', label='Eq. Altura Manométrica')
plt.xlabel("$Q$")
plt.ylabel("$H_{man}$")
plt.legend()
plt.grid()
plt.show()