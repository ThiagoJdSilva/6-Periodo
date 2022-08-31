#@Análise de Sistemas Lineares
#@Autores: Thiago José da Silva e Luiza Gomes de Castro e Sá
#@Data: 31/01/2022

import numpy as np                          #Importando a biblioteca numpy
from matplotlib import pyplot as plt        #Importando a biblioteca matplotlib responsável por plotar gráficos
import control as ct                        #Importando a biblioteca control para realizar o controle do sistema
import statistics as sta                    #Importando a biblioteca statistic
import math

# """
# ----------------------------------------------------------------------------------------------------------------------------------------------
#                                                                 Questão 1
# ----------------------------------------------------------------------------------------------------------------------------------------------
# """

# valor = list(np.arange(-100, 100,0.1)) #Criando uma lista de valores para plotar a parabola

# Ki = 0 * np.ones(len(valor)) # Valor mínimo de Ki

# Kp = [] # array para alocar os valores de Kp

# intersecao_x = [] # array área de estabilidade em x
# intersecao_y = [] # array área de estabilidade em y

# for kp in valor:
#     ki = -0.16*kp**2 + 4.31*kp + 4.47
#     if kp > -1 and kp < 0:
#         intersecao_x.append(kp + 0.1)   #valor para Kp estavavel
#         intersecao_y.append(ki + 0.1)    # valor para Ki estavel
#     elif kp > 0 and kp < 13:
#         intersecao_x.append(kp + 0.1)   #valor para Kp estavavel
#         intersecao_y.append(ki - 0.1)    # valor para Ki estavel
#     elif kp > 13 and kp < 27.52:
#         intersecao_x.append(kp - 0.1)   #valor para Kp estavavel
#         intersecao_y.append(ki - 0.1)    # valor para Ki estavel
#     Kp.append(ki)


# plt.figure(1)
# plt.title("Região estável")
# plt.plot(valor,Ki,'--k', label='Minimo de $K_i$')
# plt.plot(valor,Kp,'b', label='Valores $K_p$')
# plt.plot(valor)
# plt.xlim(-50,60)
# plt.ylim(-50,50)
# for n in list(range(-1,0)):
#     plt.vlines(x=intersecao_x[n], ymin=0, ymax=intersecao_y[n], color='m', label='Região estavel')
# for n in list(range(0,len(intersecao_y))):
#     plt.vlines(x=intersecao_x[n], ymin=0, ymax=intersecao_y[n], color='m')
# plt.xlabel("$K_p$")
# plt.ylabel("$K_i$")
# plt.legend()
# plt.grid()

# x0 = 27
# xf = 28
# y0 = 0
# yf = 2

# # ==============================
# # PLOT NORMAL
# # ==============================
# plt.plot([x0,xf], [y0,y0], 'c--')
# plt.plot([x0,xf], [yf,yf], 'c--')
# plt.plot([x0,x0], [y0,yf], 'c--')
# plt.plot([xf,xf], [y0,yf], 'c--')

# # ==============================
# # PLOT COM ZOOM
# # ==============================
# a = plt.axes([0.55, 0.18, 0.15, 0.15]) # left, bottom, width, height
# plt.xlim(x0,xf)
# plt.ylim(y0,yf)
# plt.plot(valor,Ki,'--k', label='Minimo de $K_i$')
# plt.plot(valor,Kp,'b', label='Valores $K_p$')
# for n in list(range(-1,0)):
#     plt.vlines(x=intersecao_x[n], ymin=0, ymax=intersecao_y[n], color='m', label='Região estavel')
# for n in list(range(0,len(intersecao_y))):
#     plt.vlines(x=intersecao_x[n], ymin=0, ymax=intersecao_y[n], color='m')
# plt.grid()
# plt.plot()
# plt.show()

s = ct.tf('s')          #Criando s como transferência

# tempo_final = 10   # Tempo final da análise
# periodo = 0.1
# time = np.arange(0, tempo_final, periodo) # Criando um lista que vai de 0 até 10, com período 0.1

# G_s = 83.33 / ((s+1)*(s+5)*(s+16.66))  #Definindo a função G(s)

# valores_y = [] 
# par = [] #Array para pares de Kp x Ki

# ref = 1 * np.ones(len(time))

# for k in range(len(time)):
#     if time[k] > 6 and time[k] < 8.2:
#         ref[k] = 1.2 # Definindo um degrau de 20%

# for n in intersecao_x:
#     for x in intersecao_y:
#         C_s = n + (x/s)
#         malha_fechada = (G_s*C_s)/(1+(G_s*C_s))
#         t, y = ct.forced_response(malha_fechada, time, ref, 0)
#         valores_y.append(y)
#         par.append([n, x])

# print(len(valores_y[61:81]))
# valores_ISE = []
# for n in valores_y:
#     e0 = np.abs(n[61:81] - ref[61:81]) # Calculo do erro do Controlador 
#     ISE = np.sum(e0**2)
#     valores_ISE.append(ISE)

# print(min(valores_ISE))
# print(f"O par de (Kp, Ki) que apresentaram o melhor desempenho ISE é: {par[valores_ISE.index(min(valores_ISE))]}")
# """
# ----------------------------------------------------------------------------------------------------------------------------------------------
#                                                                 Questão 2
# ----------------------------------------------------------------------------------------------------------------------------------------------
# """
# tempo_final = 100   # Tempo final da análise
# periodo = 0.1
# time = np.arange(0, tempo_final, periodo) # Criando um lista que vai de 0 até 50000, com período 0.001

# ref = 1 * np.ones(time.shape)

G_S2 = 2.11 / (s*(s+100)*(s + 1.71))
G_S3 = 2.11 / (s*(100)*(s + 1.71))
plt.figure(2)
rlist, klist = ct.root_locus(G_S2, print_gain=True, grid=True)
plt.ylim([-400,400])
plt.xlim([-150,50])
out_st_magn = ct.stability_margins(G_S2)
Gain_margin = out_st_magn
plt.show()

# K = 202.1

# Malha_Fechada_2 = (G_S2*K)/(1+(G_S2*K))
# Sistema_Segunda_Ordem_3 =(G_S3*K)/(1+(G_S3*K))
# t, y1 = ct.step_response(Malha_Fechada_2, time)
# t, y2 = ct.step_response(Sistema_Segunda_Ordem_3, time)
# plt.figure(2)
# plt.plot(t, y1, 'k', label='Controlador 1 - Questão 2.1')
# plt.plot(t, y2, 'm', label='Sistema Segunda Ordem')
# plt.ylim(0.6, 1.3)
# plt.xlim(0,10)
# plt.grid()
# plt.legend()
# plt.show()

# e0 = np.abs(y1[:] - ref[:]) # Calculo do erro do Controlador Proporcional Integral
# e1 = np.abs(y2[:] - ref[:]) # Calculo do erro do Controlador Proporcional Integral

# erro = np.abs(sum(e0)- sum(e1))
# print(erro)

# """
# ----------------------------------------------------------------------------------------------------------------------------------------------
#                                                                 Questão 3
# ----------------------------------------------------------------------------------------------------------------------------------------------
# """

G_S3 = 1/(s*(s+7))
M_S3 = (s + 10)/(s + 25.52)

malha_aberta = G_S3
malha_fechada = (G_S3*M_S3)/(1 + (G_S3*M_S3))
K1 = 1
K2 = 476.3

plt.figure(3)
ct.root_locus(malha_aberta, print_gain=True, grid=True)
mag1, fase1, omega1 = ct.bode_plot(K1*malha_aberta, plot=True, dB=True, label='malha aberta')
mag2 , fase2, omega2 = ct.bode_plot(K1*malha_fechada, plot=True, dB=True, label='malha fechada')
plt.legend()
plt.show()

Mf1 = ct.feedback(K1*G_S3)
Mf2 = ct.feedback(K2*G_S3)

time = np.arange(0, 40, 0.1)
t, y = ct.step_response(Mf1, T=time)
t, y2 = ct.step_response(Mf2, T=time)

plt.figure(4)
plt.plot(t,y, 'r', label='K = 1')
plt.plot(t,y2,'b', label='K = 476.3')
plt.plot(t,0.98*np.ones(len(time)),'--k')
plt.plot(t,1.02*np.ones(len(time)),'--k')
plt.grid()
plt.legend()
plt.show()









