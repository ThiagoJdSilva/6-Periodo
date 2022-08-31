#Laboratório de Análise de Sistemas Lineares
#@Autores: Thiago José da Silva e Luiza Gomes de Castro e Sá
#@Data: 29/01/2022

import numpy as np                          #Importando a biblioteca numpy
from matplotlib import pyplot as plt        #Importando a biblioteca matplotlib responsável por plotar gráficos
import control as ct                        #Importando a biblioteca control para realizar o controle do sistema
import statistics as sta                    #Importando a biblioteca statistic
import math

plt.close('all') #Fechando todas as abas de gráficos abertas

def model_update(t,x,u,params):
    Temp = x[0] 
    Qr = x[1] 
    Td = u[1]
    Ft = F*(1+Td)
    u = u[0] 
    dT = (1/(C*V))*(Qr+(C*Ft*(T_inicial-Temp))-Qq-(h*A_ext*(Temp-T_inicial)))
    dQr = (1/(12.5*h*np.pi))*(-Qr+(R*u)) 
    return dT,dQr # Retorna Temperatura e o Qr

def model_output(t,x,u,params):   
    # Para o caso em estudo, a saída do sistema, y, é o estado X[0], ou seja, y = theta.
    y1 = x[0]
    y2 = x[1]

    return y1,y2

def valor_Qr_Entrada():
    """
    Retorna o valor de Qr em regime permanente para o ponto de operação igual a 60ºC
    """
    return -((C*F*(T_inicial-ponto_de_operacao)) - Qq - (h*A_ext*(ponto_de_operacao - T_inicial))) #Valor de Qr para o ponto de operação

def valor_u_entrada():
    """
    Retorna o valor de u em regime permanente para o ponto de operação igual a 60ºC
    """
    return Qr/R 

def equacao(vetor1, vetor2, c_angular, k):
    """
    Retorna a equação da reta gerada pelo método de Ziegler Nichols
    """
    return c_angular*vetor1 - c_angular*vetor1[k] + vetor2[k]

def valor_K(vetor1, vetor2, ponto_de_operacao, variavel_de_controle, k):
    """
    Retorna o valor da variação de Y, variação u, e o valor de K
    """
    return vetor1[k]-ponto_de_operacao, vetor2[k]-variavel_de_controle, abs(vetor1[k]-ponto_de_operacao)/abs(vetor2[k]-variavel_de_controle)

def coeficiente_angular(valor_y, valor_x, j, k):
    """
    Retorna o coeficiente angular da reta tangente a ao degrau da curva
    """
    return abs(valor_y[j]-valor_y[k]) / abs(valor_x[j]-valor_x[k])

def saturation_output(t,x,u,params):  
    if u < 0:
        return 0  # Quando u assume valores negativos a função retorna 0
    else:
        return u  # Quando u assume valores positivos, a função retorna o próprio u

# Definindo a saturação da entrada, segundo biblioteca----------------------------------------------------------------
Saturacao = ct.NonlinearIOSystem(0, saturation_output, states=0, name='saturacao', inputs = ('u'), outputs = ('y'))

Reator_Quimico = ct.NonlinearIOSystem(model_update, model_output, states=2, name='Reator_Quimico', inputs=('u','Td'), outputs=('y1','y2'))

#----------------------------------------------------------------------------------------------------------------------------------------
" Definindo as varíaveis do sistema"
#----------------------------------------------------------------------------------------------------------------------------------------
V = 10         # Volume do reator
C = 4500       # Capacidade térmica da solução
F = 0.5        # Vazao volumetrica
h = 15         # Coeficiente de convecção
A_ext = 31.4    # Superficie do tanque   
Qq = 7900      # Energia necessária para catalização
R = 10000      # Valor resistência térmica
T_inicial = 20 # Temperatura inicial
ponto_de_operacao = 60     # Ponto de operação
Qr = valor_Qr_Entrada()
U = valor_u_entrada()
print(f"O valor do Qr de entrada é: {Qr}")
print(f"O valor do U de entrada é: {U}")

"""
Definindo os vetores de tempo
"""
periodo = 0.1
tempo_final = 30000   # Tempo final da análise
time = np.arange(0, tempo_final, periodo) # Criando um lista que vai de 0 até 50, com período 0.1
Td = np.zeros(len(time))   #criando vetor Td com zeros, de tamanho t

ref_Qr = Qr * np.ones(time.shape)
ref_U = U * np.ones(time.shape)
ref_T = ponto_de_operacao * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao ponto de operação

"""
Criando os degraus para validação:
"""
u_degrau = U * np.ones(time.shape)

for k in range(len(time)):
    if time[k] < 15000:
        u_degrau[k] =  U
    else:
        u_degrau[k] = U + 0.2*U

t, y = ct.input_output_response(Reator_Quimico, time, [u_degrau, Td], [T_inicial, 0])

#----------------------------------------------------------------------------------------------------------------------------------------
"Plot do sistema com aplicação de degraus na entrada"
#----------------------------------------------------------------------------------------------------------------------------------------
plt.figure(1)
plt.subplot(3,1,1)
plt.plot(t,ref_T,'--k', label='ref')
plt.plot(t,y[0],'b', label='${\\theta(t)}$')
plt.ylabel('${\\theta(t)}$${[\\degree C]}$')
plt.xlim(0,tempo_final)
plt.ylim(20,75)
plt.legend()
plt.grid()

plt.subplot(3,1,2)
plt.plot(t,ref_Qr,'--k', label='ref')
plt.plot(t,y[1],'b', label='${\\theta(t)}$')
plt.ylabel('$Q_r(t)$[J]')
plt.xlim(0,tempo_final)
plt.ylim(20,150000)
plt.legend()
plt.grid()

plt.subplot(3,1,3)
plt.plot(time,u_degrau,'b', label='${u(t)}$')
plt.ylabel('${u(t)[A]}$')
plt.plot(t,ref_U,'--k', label='ref')
plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,tempo_final)
plt.ylim(10,15)
plt.grid()
plt.show()
#----------------------------------------------------------------------------------------------------------------------------------------
"""
Realizando o método de Ziegler Nichols
"""
coeficiente_degrau = coeficiente_angular(y[0], time, 150500, 151000)
equacao_degrau = equacao(time, y[0], coeficiente_degrau, 150500)

plt.figure(2)
plt.plot(time,y[0],'b',label='$\\theta(t)$')
plt.plot(time,equacao_degrau,'r',label='')
plt.plot(time,y[0][149000] * np.ones(time.shape),'--k',label='ref resposta')
plt.plot(time,y[0][180000] * np.ones(time.shape),'--g',label='ref resposta degrau')
plt.grid()
plt.ylabel('$\\theta(t)[\\degree C]$')
plt.xlabel('Tempo[s]')
plt.legend()
plt.ylim(50,75)
plt.xlim(14000,18000)

x0 = 14980
xf = 15100
y0 = 59.9
yf = 60.9

x0_ = 14925
xf_ = 15150
y0_ = 59.5
yf_ = 60.8

# ==============================
# PLOT NORMAL
# ==============================
plt.plot([x0_,xf_], [y0_,y0_], 'c--')
plt.plot([x0_,xf_], [yf_,yf_], 'c--')
plt.plot([x0_,x0_], [y0_,yf_], 'c--')
plt.plot([xf_,xf_], [y0_,yf_], 'c--')

# ==============================
# PLOT COM ZOOM
# ==============================
a = plt.axes([0.55, 0.18, 0.15, 0.15]) # left, bottom, width, height
plt.xlim(x0,xf)
plt.ylim(y0,yf)
plt.grid()
plt.plot(time,y[0],'b',label='$g(s)$')
plt.plot(time,equacao_degrau,'r',label='ref degrau')
plt.plot(time,y[0][149000] * np.ones(time.shape),'--k',label='ref resposta degrau')
plt.show()

tempo_final_validacao = 60000 # Tempo de execução do sistema
tempo_validacao = np.arange(0,tempo_final_validacao,periodo)  # Array do tempo de execução do sistema

u_degrau_validacao = U * np.ones(tempo_validacao.shape)
#Criando diferentes degraus para realizar a validação
for k in range(len(tempo_validacao)):
    if tempo_validacao[k] < 12000:
        u_degrau_validacao[k] = U
    elif tempo_validacao[k] < 24000:
        u_degrau_validacao[k] = 1.1*U
    elif tempo_validacao[k] < 36000:
        u_degrau_validacao[k] = U
    elif tempo_validacao[k] < 48000:
        u_degrau_validacao[k] = 0.9*U
    else:
        u_degrau_validacao[k] = U

Td = np.zeros(len(tempo_validacao))

tau = 15654 - 15014
theta, u, K1 = valor_K(y[0], u_degrau, ponto_de_operacao, U, 295000)
num_pos = [K1]                # Numerador da função de transferência 
den_pos = [tau, 1]               # Denominador da função de transferência 
G_pos = ct.tf(num_pos,den_pos) 
print(f"A equação de 1ª ordem é: {G_pos}")
t, mz = ct.forced_response(G_pos,T= tempo_validacao, U = u_degrau_validacao - U) # Simulando a saída de um sistema linear
t_valid, ms = ct.input_output_response(Reator_Quimico, tempo_validacao, [u_degrau_validacao, Td], [ponto_de_operacao,Qr])     #resposta do sistema original para os degraus aplicados

print(f"Variação em y: {theta}")
print(f"Variação em U: {u}")
print(f"Valor de K: {K1}")
print(f"Valor de Tau: {tau}")

plt.figure(3)
plt.subplot(211)
plt.plot(t_valid, ms[0], 'm', label = 'Sistema não linear')
plt.plot(t_valid, mz+ponto_de_operacao, 'g', label = 'Modelo 1ª Ordem')
plt.plot(t_valid, ponto_de_operacao*np.ones(t_valid.shape), '--k', label = 'Ref.')
plt.ylabel('$\\theta(t) [\\degree C]$')
plt.xlim(0,60000)
plt.grid()
plt.legend()

plt.subplot(212)
plt.plot(t_valid, u_degrau_validacao, 'c',label = 'Sinal de controle')
plt.plot(t_valid, U*np.ones(t_valid.shape), '--k', label = 'Ref.')
plt.xlabel('Tempo [s]')
plt.ylabel('u(t) [A]')
plt.xlim(0,60000)
plt.grid()
plt.legend()
plt.show()

"""
Realizando a aproximação de Padé de 5ª ordem
"""
atraso = 8
n = 5

N = ct.pade(atraso,n) #aproximação de padé
Gd = ct.tf(N[0], N[1])
Atraso = ct.tf2io(Gd, name = 'atraso', inputs = 'u', outputs = 'y')    #criando bloco para o atraso
GdI = Gd * G_pos
print(GdI)

#Controlador P
Kp_chr = (0.3*tau)/(K1*atraso)       #calculando o valor do Kp para um controlador P
print(f"Valor de Kp - controlador P: {Kp_chr}")
Gp_chr = ct.tf(Kp_chr,1)                    #definindo a função de transferência do controlador P pelo método chr 
print(f"Controlador P - CHR: {Gp_chr}") 

#controlador PI
Kpi_chr = ((0.6*tau)/(K1*atraso))/10 #calculando o valor de Kp para um controlador PI
print(f"Valor de Kp - controlador PI: {Kpi_chr}")
Ti_chr = 4*atraso   #calculando o tempo integral para um controlador PI
print(f"Valor de Ti - controlador PI: {Ti_chr}")
num = [Kpi_chr*Ti_chr,Kpi_chr]    #definindo o numerador da função de transferência
den = [Ti_chr,0]                      #definindo o denominador da função de transferência
Gpi_chr = ct.tf(num,den)        #calculando a função de transferência do controlador PI pelo método chr     
print(f"Controlador PI - CHR: {Gpi_chr}") 

controlador_P = ct.tf2io(Gp_chr, name='controlador_P', inputs='u', outputs='y')   #controlar P
controlador_PI = ct.tf2io(Gpi_chr, name='controlador_PI', inputs='u', outputs='y')   #controlar PI

# Controlador 1 (P) 
SysP = ct.InterconnectedSystem(
    (controlador_P,Atraso,Saturacao,Reator_Quimico), name='sysdelayed1',
    connections = (('controlador_P.u','-Reator_Quimico.y1'),('saturacao.u','controlador_P.y'),
                    ('atraso.u','saturacao.y'),('Reator_Quimico.u','atraso.y')),
    inplist = ('controlador_P.u','Reator_Quimico.Td','saturacao.u'),
    inputs = ('ref','Td','u0'),
    outlist = ('Reator_Quimico.y1','Reator_Quimico.y2', 'atraso.u'),
    outputs = ('y1','y2','u'))

# Controlador 2 (PI) 
SysPI = ct.InterconnectedSystem(
    (controlador_PI,Atraso,Saturacao,Reator_Quimico), name='sysdelayed2',
    connections = (('controlador_PI.u','-Reator_Quimico.y1'),('saturacao.u','controlador_PI.y'),
                    ('atraso.u','saturacao.y'),('Reator_Quimico.u','atraso.y')),
    inplist = ('controlador_PI.u','Reator_Quimico.Td','saturacao.u'),
    inputs = ('ref','Td','u0'),
    outlist = ('Reator_Quimico.y1','Reator_Quimico.y2', 'atraso.u'),
    outputs = ('y1','y2','u'))

"""
    Determinando funções para os controladores P e PI
"""
#Ganho de malha
Lp =  Gp_chr * G_pos     #ganho de malha controlador P
Lpi = Gpi_chr * G_pos    #ganho de malha controlador PI
print(f"Ganho de malha controlador P: {Lp} \n Ganho de malha controlador PI: {Lpi}")

#função sensitividade
Sp = 1/(1+Lp)       #função sensitividade controlador P
Spi = 1/(1+Lpi)     #função sensitividade controlador PI
print(f"Função sensitividade controlador P: {Sp} \n Função sensitividade controlador PI: {Spi}")

#função sensitividade complementar
Cp = Lp/(1+Lp)      #função sensitividade complementar controlador P
Cpi = Lpi/(1+Lpi)   #função sensitividade complementar controlador PI
print(f"Função sensitividade complementar controlador P: {Cp} \n Função sensitividade complementar controlador PI: {Cpi}")

"""
Resposta a malha fechada
"""
# Resposta do Sistema a um impulso

novo_tempo_final = 15000  #tempo final da simulação de validação
tempo_impulso = np.arange(0, novo_tempo_final, periodo)   #vetor tempo da simulação de validação
ref = ponto_de_operacao*np.ones(tempo_impulso.shape)    #vetor com o ponto de operação
Td = np.zeros(len(tempo_impulso))   #vetor de perturbação volumetrica

u_impulso_controlador = U * np.ones(tempo_impulso.shape)

#Colocando um impulso de 2x para 1 segundo do código
for k in range(len(tempo_impulso)):
    if tempo_impulso[k] >= 5000 and tempo_impulso[k] <= 5010:
        u_impulso_controlador[k] =  2*U
    else:
        u_impulso_controlador[k] = U
#Função da biblioteca que retorna a resposta da malha fechada para impulso e degrau (Controlador 1) 
t2, y1 = ct.input_output_response(SysP, tempo_impulso, [ref, Td, u_impulso_controlador], [0,0,0,0,0,ponto_de_operacao,Qr]) 
#Função da biblioteca que retorna a resposta da malha fechada para impulso e degrau (Controlador 2) 
t2, y2 = ct.input_output_response(SysPI, tempo_impulso, [ref, Td, u_impulso_controlador], [0,0,0,0,0,0,ponto_de_operacao,Qr]) 

plt.figure(4)
plt.subplot(311)
plt.plot(t2, y1[0], 'r', label = 'Controlador P')
plt.plot(t2, y2[0], 'b', label = 'Controlador PI')
plt.plot(t2, ponto_de_operacao*np.ones(t2.shape), '--k', label = 'Ref.' )
plt.ylabel('$\\theta(t)[°C]$')
plt.xlim(0,15000)
plt.grid()
plt.legend()

plt.subplot(312)
plt.plot(t2, y1[1], 'r', label = 'Controlador P')
plt.plot(t2, y2[1], 'b', label = 'Controlador PI')
plt.plot(t2, Qr*np.ones(t2.shape), '--k', label = 'Ref.' )
plt.ylabel('Qr(t) [J]')
plt.xlim(0,15000)
plt.grid()
plt.legend()

plt.subplot(313)
plt.plot(t2, u_impulso_controlador, label="Entrada impulso")
plt.plot(t2, U*np.ones(t2.shape), '--k', label = 'Ref.')
plt.ylabel('u(t) [A]')
plt.xlabel('Tempo [s]')
plt.xlim(0,15000)
plt.grid()
plt.legend()
plt.show()

# Resposta do Sistema a um degrau com amplitude 0.25
u_degrau_controlador = U * np.ones(tempo_impulso.shape)

#Acrescentando um degrau de 25% para o tempo de 5000 até 10000.
for k in range(len(tempo_impulso)):
    if tempo_impulso[k] < 5000:
        u_degrau_controlador[k] =  U
    elif tempo_impulso[k] < 10000:
        u_degrau_controlador[k] =  1.25*U
    else:
        u_degrau_controlador[k] = U

#Função da biblioteca que retorna a resposta da malha fechada para impulso e degrau (Controlador 1) 
t2, y1 = ct.input_output_response(SysP, tempo_impulso, [ref, Td, u_degrau_controlador], [0,0,0,0,0,ponto_de_operacao,Qr]) 
#Função da biblioteca que retorna a resposta da malha fechada para impulso e degrau (Controlador 2) 
t2, y2 = ct.input_output_response(SysPI, tempo_impulso, [ref, Td, u_degrau_controlador], [0,0,0,0,0,0,ponto_de_operacao,Qr]) 

plt.figure(5)
plt.subplot(311)
plt.plot(t2, y1[0], 'r', label = 'Controlador P')
plt.plot(t2, y2[0], 'b', label = 'Controlador PI')
plt.plot(t2, ponto_de_operacao*np.ones(t2.shape), '--k', label = 'Ref.' )
plt.ylabel('$\\theta(t) [°C]$')
plt.xlim(0,15000)
plt.grid()
plt.legend()

plt.subplot(312)
plt.plot(t2, y1[1], 'r', label = 'Controlador P')
plt.plot(t2, y2[1], 'b', label = 'Controlador PI')
plt.plot(t2, Qr*np.ones(t2.shape), '--k', label = 'Ref.' )
plt.ylabel('Qr(t)[J]')
plt.xlim(0,15000)
plt.grid()
plt.legend()

plt.subplot(313)
plt.plot(t2, u_degrau_controlador, label="Entrada degrau")
plt.plot(t2, U*np.ones(t2.shape), '--k', label = 'Ref.')
plt.ylabel('u(t)[A]')
plt.xlabel('Tempo [s]')
plt.xlim(0,15000)
plt.grid()
plt.legend()
plt.show()

# Resposta do Sistema a perturbação na vazão volumétrica
F_volumetrico = F * np.ones(tempo_impulso.shape)

for k in range(len(tempo_impulso)):
    if tempo_impulso[k] < 5000:
        F_volumetrico[k] = F
    elif tempo_impulso[k] < 10000:
        F_volumetrico[k] =  1.3*F
    else:
        F_volumetrico[k] = F

Qr_vol = np.ones(len(tempo_impulso))
for j in range(len(tempo_impulso)):
    Qr_vol[j] = -((C*F_volumetrico[j]*(T_inicial-ponto_de_operacao))-Qq-(h*A_ext*(ponto_de_operacao-T_inicial)))

u_vol = np.ones(len(tempo_impulso))
for i in range(len(tempo_impulso)):
    u_vol[i] = Qr_vol[i]/R

#Função da biblioteca que retorna a resposta da malha fechada para variação da vazão volumétrica (Controlador 1) 
t3, y3 = ct.input_output_response(SysP, tempo_impulso, [ref, Td, u_vol], [0,0,0,0,0,ponto_de_operacao,Qr])   
#Função da biblioteca que retorna a resposta da malha fechada para variação da vazão volumétrica (Controlador 2) 
t3, y4 = ct.input_output_response(SysPI, tempo_impulso, [ref, Td, u_vol], [0,0,0,0,0,0,ponto_de_operacao,Qr]) 

plt.figure(6)
plt.subplot(311)
plt.plot(t3, y3[0], 'r', label = 'Controlador P')
plt.plot(t3, y4[0], 'b', label = 'Controlador PI')
plt.plot(t3, ponto_de_operacao*np.ones(t3.shape), '--k', label = 'Ref.' )
plt.ylabel('$\\theta(t)[°C]$')
plt.xlim(0,15000)
plt.grid()
plt.legend()

plt.subplot(312)
plt.plot(t3, y3[1], 'r', label = 'Controlador P')
plt.plot(t3, y4[1], 'b', label = 'Controlador PI')
plt.plot(t3, Qr*np.ones(t3.shape), '--k', label = 'Ref.' )
plt.ylabel('Qr(t) [J]')
plt.xlim(0,15000)
plt.grid()
plt.legend()

plt.subplot(313)
plt.plot(t3, F_volumetrico, 'c', label = 'Perturbação volumétrica')
plt.plot(t3, U*np.ones(t3.shape), '--k', label = 'Ref.')
plt.ylabel('u [A]')
plt.xlabel('Tempo [s]')
plt.xlim(0,15000)
plt.ylim(0, 1)
plt.grid()
plt.legend()
plt.show()
