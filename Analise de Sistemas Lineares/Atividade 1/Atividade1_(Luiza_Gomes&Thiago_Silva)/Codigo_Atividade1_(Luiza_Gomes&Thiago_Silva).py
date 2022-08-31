# Atividade 01 - Análise de Sistema Lineares
# Thiago José da Silva
# Luiza Gomes de Castro e Sá

import numpy as np
from matplotlib import pyplot as plt
import control as ct
from numpy.core.numeric import NaN

plt.close('all')

"""
                                            Sistema 1.1 - Processo Petroquímico
"""

# # ==============================
# # Variáveis do sistema
# # ==============================

s=ct.tf("s") # Frequência

# # ==============================
# # Variéveis e vetores de tempo/entrada/referência
# # ==============================
periodo = 0.001
tempo_final = 100 # Tempo de execução do sistema
time = np.arange(0,tempo_final,periodo)  # Array do tempo de execução do sistema
U = 40 # Entrada do Sistema
u = U * np.ones(time.shape) # Array preenchido com a entrada para o ponto de operação que não será alterada
ref = 2 * np.ones(time.shape) # Array de referência para o ponto de operação

def valor_K_estacionario():
    """
    Determinando o valor de K para s inicial igual a 0 (zero).
    """
    Vm = 2
    Q = 40
    s = 0
    return (Vm *  (s**3 + 36*s**2 + 280*s + 3000)/(Q * (100*s + 150)))

K = valor_K_estacionario()

funcao_de_transferencia = (K*(100*s + 150))/((s**3) + (36*s**2) + (280*s) + 3000) # Equação de transferência
tempo_eq_transf, saida_eq_transf = ct.forced_response(funcao_de_transferencia,T=time, U=u) # saida_eq_transf é a saida do medidor e tempo_eq_transf é o tempo de integraçào

# # ==============================
# # Plot equação de tranferência com controlador proporcional
# # ==============================
plt.figure(1)
plt.subplot(211) 
plt.title('Equação de Transfêrencia')
plt.plot(tempo_eq_transf,saida_eq_transf,'b',label='$g(s)$')
plt.plot(time,ref,'--k',label='ref')
plt.grid()
plt.ylabel('$V_m(s)$')
plt.legend()
plt.ylim(-2,11)
plt.xlim(0,3)

plt.subplot(212) 
plt.plot(tempo_eq_transf,u,'b', label = '$Q(s)$')
plt.grid()
plt.ylabel('Q(s)')
plt.xlabel('Tempo(s)')
plt.legend()
plt.xlim(0,30)
plt.ylim(0,80)
plt.show()

"""
Criando os degraus positivos para 45 Kg/s
"""
u_degrau_positivo = U * np.ones(len(time)) # Foi criado um novo array para os degraus 
u_degrau_positivo = np.array_split(u_degrau_positivo, 2) # O array foi separado em 2 novos arrays
# Foi adicionado os degraus para cada intervalo
u_degrau_positivo[0].fill(U)
u_degrau_positivo[1].fill(45)
u_degrau_positivo = np.concatenate([u_degrau_positivo[0], u_degrau_positivo[1]]) # Unindo os arrays dos degraus

tempo_degrau_positivo, saida_degrau_positivo = ct.forced_response(funcao_de_transferencia,T=time, U=u_degrau_positivo) #saida_degrau_positivo é a saida do medidor e tempo_degrau_positivo é o tempo de integraçào
ref_degrau = 2.25 * np.ones(time.shape) # Array de referência para o ponto de operação

# # ==============================
# # Plot equação de tranferência com controlador proporcional degrau positivo
# # ==============================
plt.figure(2)
plt.subplot(211) 
plt.title('Equação de Transfêrencia com degrau positivo')
plt.plot(tempo_degrau_positivo,saida_degrau_positivo,'r',label='$g(s)$')
plt.plot(time,ref_degrau,'--g',label='ref degrau')
plt.plot(time,ref,'--k',label='ref')
plt.grid()
plt.ylabel('$-$')
plt.legend()
plt.ylim(1.5,3.5)
plt.xlim(49.5,53)

plt.subplot(212) 
plt.plot(tempo_degrau_positivo,u_degrau_positivo,'r', label = '$Q(s)$')
plt.grid()
plt.ylabel('-')
plt.xlabel('-')
plt.legend()
plt.ylim(35,50)
plt.xlim(49.5,53)
plt.show()

compensação_serie = 0.2/(s + 4) #Equação da compensação
tempo_comp_degrau_positivo, saida_comp_degrau_positivo = ct.forced_response(compensação_serie,T=time, U=u_degrau_positivo)

# # ==============================
# # Plot equação de tranferência com controlador proporcional degrau positivo com compensação 
# # ==============================
plt.figure(3)
plt.title('Equação de Transfêrencia com compensação degrau positivo')
plt.plot(tempo_degrau_positivo,saida_degrau_positivo,'r',label='$g(s)$')
plt.plot(tempo_comp_degrau_positivo, saida_comp_degrau_positivo,'b',label='$compensação$')
plt.plot(time,ref_degrau,'--g',label='ref degrau')
plt.plot(time,ref,'--k',label='ref')
plt.grid()
plt.ylabel('$-$')
plt.legend()
plt.ylim(1.5,3.5)
plt.xlim(49.5,53)
plt.show()

"""
Criando dos degraus negativos de 45Kg/s para 35 Kg/s
"""
u_degrau_negativo = U * np.ones(len(time)) # Foi criado um novo array para os degraus 
u_degrau_negativo = np.array_split(u_degrau_negativo, 2) # O array foi separado em 2 novos arrays
# Foi adicionado os degraus para cada intervalo
u_degrau_negativo[0].fill(45)
u_degrau_negativo[1].fill(35)
u_degrau_negativo = np.concatenate([u_degrau_negativo[0], u_degrau_negativo[1]]) # Faço o concat dos arrays dos degraus

tempo_degrau_negativo, saida_degrau_negativo = ct.forced_response(funcao_de_transferencia,T=time, U=u_degrau_negativo) #tempo_degrau_positivo é a saida do medidor e tempo_degrau_positivo é otempo de integraçào
ref_degrau = 1.75 * np.ones(time.shape) # Array de referência para o ponto de operação

tempo_comp_degrau_negativo, saida_comp_degrau_negativo = ct.forced_response(compensação_serie,T=time, U=u_degrau_negativo)
ref2 = 2.25 * np.ones(time.shape) # Array de referência para o ponto de operação
ref3 = 1.71 * np.ones(time.shape) # Array de referência para o ponto de operação
ref4 = 1.785 * np.ones(time.shape) # Array de referência para o ponto de operação

# # ==============================
# # Plot equação de tranferência com controlador proporcional degrau negativo de 45Kg/s para 35Kg/s sem compensação
# # ==============================
plt.figure(4)
plt.subplot(211)
plt.title('Equação de Transfêrencia com degrau negativo')
plt.plot(tempo_degrau_negativo,saida_degrau_negativo,'r',label='$g(s)$')
plt.plot(time,ref_degrau,'--g',label='ref degrau')
plt.plot(time,ref3,'--y',label='$-2%$')
plt.plot(time,ref4,'--c',label='$+2%$')
plt.plot(time,ref2,'--k',label='ref')
plt.grid()
plt.ylabel('$-$')
plt.legend()
plt.ylim(-0.5, 3)
plt.xlim(49.5,52)

x0 = 51.3
xf = 51.7
y0 = 1.70
yf = 1.80

# ==============================
# PLOT NORMAL
# ==============================

plt.plot([x0,xf], [y0,y0], 'k--')
plt.plot([x0,xf], [yf,yf], 'k--')
plt.plot([x0,x0], [y0,yf], 'k--')
plt.plot([xf,xf], [y0,yf], 'k--')

# ==============================
# PLOT COM ZOOM
# ==============================

a = plt.axes([0.35, 0.5, 0.15, 0.15]) # left, right, width, height
plt.xlim(x0,xf)
plt.ylim(y0,yf)
plt.grid()
plt.plot(tempo_degrau_negativo,saida_degrau_negativo,'r',label='$g(s)$')
plt.plot(time,ref_degrau,'--g',label='ref degrau')
plt.plot(time,ref3,'--y',label='$-2\\%$')
plt.plot(time,ref4,'--c',label='$+2\\%$')
plt.plot(time,ref2,'--k',label='ref')

plt.subplot(212) 
plt.plot(tempo_comp_degrau_negativo,u_degrau_negativo,'r', label = '$Q(t)$')
plt.grid()
plt.ylabel('Q(s)')
plt.xlabel('Tempo')
plt.legend()
plt.ylim(30,50)
plt.xlim(49.5,52)
plt.show()

# # ==============================
# # Plot equação de tranferência com controlador proporcional degrau negativo de 45Kg/s para 35Kg/s com compensação
# # ==============================
plt.figure(5)
plt.subplot(211)
plt.title('Equação de Transfêrencia com compensação degrau negativo')
# plt.plot(tempo_degrau_negativo,saida_degrau_negativo,'r',label='$-$')
plt.plot(tempo_comp_degrau_negativo, saida_comp_degrau_negativo,'b',label='$g(s)$')
plt.plot(time,ref_degrau,'--g',label='ref degrau')
plt.plot(time,ref3,'--y',label='$-2\\%$')
plt.plot(time,ref4,'--c',label='$+2\\%$')
plt.plot(time,ref2,'--k',label='ref')
plt.grid()
plt.ylabel('$-$')
plt.legend()
plt.ylim(-0.5, 3)
plt.xlim(49.5,52)

x0 = 50.9
xf = 51.1
y0 = 1.70
yf = 1.80

# ==============================
# PLOT NORMAL
# ==============================

plt.plot([x0,xf], [y0,y0], 'k--')
plt.plot([x0,xf], [yf,yf], 'k--')
plt.plot([x0,x0], [y0,yf], 'k--')
plt.plot([xf,xf], [y0,yf], 'k--')

# ==============================
# PLOT COM ZOOM
# ==============================

a = plt.axes([0.35, 0.5, 0.15, 0.15]) # left, right, width, height
plt.xlim(x0,xf)
plt.ylim(y0,yf)
plt.grid()
plt.plot(tempo_comp_degrau_negativo, saida_comp_degrau_negativo,'b',label='$-$')
plt.plot(time,ref_degrau,'--g',label='ref degrau')
plt.plot(time,ref2,'--k',label='ref')

plt.subplot(212) 
plt.plot(tempo_comp_degrau_negativo,u_degrau_negativo,'r', label = '$Q(t)$')
plt.grid()
plt.ylabel('Q(s)')
plt.xlabel('Tempo')
plt.legend()
plt.ylim(30,50)
plt.xlim(49.5,52)
plt.show()

###                                                Fim do sistema 1.1                                                                        ###

"""
                                                Sistema 1.2 - Interconexão de sistemas
"""
u = 20
u_degraus_2 = u * np.ones(len(time)) # Foi criado um novo array para os degraus 
u_degraus_2 = np.array_split(u_degraus_2, 2) # O array foi separado em 2 novos arrays
# Foi adicionado os degraus para cada intervalo
u_degraus_2[0].fill(20)
u_degraus_2[1].fill(25)
u_degraus_2 = np.concatenate([u_degraus_2[0], u_degraus_2[1]]) # unindo os arrays dos degraus

def funcao_transfer_questao2(tau, u, time):
    funcao_transferencia = 20 /((tau*s + 10)*(s*2 + s + 1))  #Equação de transferência
    return ct.forced_response(funcao_transferencia,T=time, U=u)

valores_tau = [0.2, 0.5, 1.0, 5.0, 10.0, 20.0, 30.0]

tempo_tau0, saida_tau0 = funcao_transfer_questao2(valores_tau[0], u_degraus_2, time)
tempo_tau1, saida_tau1 = funcao_transfer_questao2(valores_tau[1], u_degraus_2, time)
tempo_tau2, saida_tau2 = funcao_transfer_questao2(valores_tau[2], u_degraus_2, time)
tempo_tau3, saida_tau3 = funcao_transfer_questao2(valores_tau[3], u_degraus_2, time)
tempo_tau4, saida_tau4 = funcao_transfer_questao2(valores_tau[4], u_degraus_2, time)
tempo_tau5, saida_tau5 = funcao_transfer_questao2(valores_tau[5], u_degraus_2, time)
tempo_tau6, saida_tau6 = funcao_transfer_questao2(valores_tau[6], u_degraus_2, time)

# # ==============================
# # Plot equação de tranferência com diferentes valores de tau
# # ==============================
plt.figure(6)
plt.subplot(211)
plt.title('Equação de transferência para diferentes valores de $\\tau$')
plt.plot(tempo_tau0, saida_tau0,'r',label='$\\tau$ = 0.2')
plt.plot(tempo_tau1, saida_tau1,'b',label='$\\tau$ = 0.5')
plt.plot(tempo_tau2, saida_tau2,'y',label='$\\tau$ = 1')
plt.plot(tempo_tau3, saida_tau3,'g',label='$\\tau$ = 5')
plt.plot(tempo_tau4, saida_tau4,'k',label='$\\tau$ = 10')
plt.plot(tempo_tau5, saida_tau5,'c',label='$\\tau$ = 20')
plt.plot(tempo_tau6, saida_tau6,'m',label='$\\tau$ = 30')
plt.grid()
plt.ylabel('$y(t)$')
plt.legend()
plt.ylim(39,51)
plt.xlim(49,65)

plt.subplot(212)
plt.plot(time,u_degraus_2,'b',label='$u(t)$')
plt.grid()
plt.ylabel('u(t) [ ? ] ')
plt.xlabel('Tempo [s]')
plt.legend()
plt.ylim(15,30)
plt.xlim(49,65)
plt.show()
###                                                 Fim do sistema 1.2                                                                       ###



"""
                                                 Sistema 1.3 - Sistema Realimentado
"""

u = 5
"""
Realizando a criação dos degraus
"""
u_degraus_3 = u * np.ones(len(time)) # Foi criado um novo array para os degraus 
u_degraus_3 = np.array_split(u_degraus_3, 2) # O array foi separado em 2 novos arrays
# Foi adicionado os degraus para cada intervalo
u_degraus_3[0].fill(5)
u_degraus_3[1].fill(6)
u_degraus_3 = np.concatenate([u_degraus_3[0], u_degraus_3[1]]) # Unindo os arrays dos degraus

k = 20 # Valor de K para um coeficiente de amortecimento = 0.5
ref3 = 6 * np.ones(time.shape) # Array de referência para o ponto de operação

def funcao_transfer_questao3(k, u, time):
    funcao_transferência = k*5/((s**2) + (10*s) + 5*k)  #Equação de transferência
    return ct.forced_response(funcao_transferência,T=time, U=u) 

tempo_k20, saida_k20 = funcao_transfer_questao3(k, u_degraus_3, time)
tempo_k40, saida_k40 = funcao_transfer_questao3(40, u_degraus_3, time)
tempo_k60, saida_k60 = funcao_transfer_questao3(60, u_degraus_3, time)
tempo_k80, saida_k80 = funcao_transfer_questao3(80, u_degraus_3, time)
tempo_k100, saida_k100 = funcao_transfer_questao3(100, u_degraus_3, time)

# # ==============================
# # Plot equação de tranferência para diferentes valores de K
# # ==============================
plt.figure(7)
plt.subplot(211) 
plt.title('Equação de Transfêrencia para diferentes valores de K')
plt.plot(tempo_k20,saida_k20,'b',label='$K = 20$')
plt.plot(tempo_k40,saida_k40,'y',label='$K = 40$')
plt.plot(tempo_k60,saida_k60,'g',label='$K = 60$')
plt.plot(tempo_k80,saida_k80,'r',label='$K = 80$')
plt.plot(tempo_k100,saida_k100,'c',label='$K = 100$')
plt.plot(time,ref3,'--k',label='ref')
plt.grid()
plt.ylabel('$-$')
plt.xlabel('Tempo [s]')
plt.legend(loc='right')
plt.ylim(4.9,6.5)
plt.xlim(49.9,51)

plt.subplot(212) 
plt.plot(tempo_k20,u_degraus_3,'b', label = '$-$')
plt.grid()
plt.ylabel('-')
plt.xlabel('-')
plt.legend()
plt.xlim(49.9,51)
plt.ylim(2,8)
plt.show()
###                                                 Fim do sistema 1.3                                                                       ###

