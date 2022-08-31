# Relatório 02 - Laboratório de Análise de Sistema Lineares
# Modelagem Caixa Preta
# João Victor Tavares Santos
# Maria Clara Wasconcelos

import numpy as np
from matplotlib import pyplot as plt
import control as ct
from numpy.core.numeric import NaN
import math

plt.close('all')

# # ==============================
# # Variéveis do sistema
# # ==============================
operationPoint = 50 # Ponto de operação
operationPointRad = math.radians(operationPoint) # Ponto de operação
L_A = 0.154 # Largura da placa movel de aluminio 
L_1 = 0.155 # Comprimento da placa abaixo do eixo de rotação
L_T = 0.270 # Comprimento da placa total
D_cm = 0.02 # Distância do centro de massa da placa
M_T = 0.1 # Massa total da placa
rho = 1.23  # Massa especifica da solução
C_A = 2.05 # Coeficiente de arrasto
u = 5.0 # Coeficiente de atrito viscoso
g = 9.81 # Gravidade
intialPosition = 0 # Posição inicial
K_1 = (D_cm * rho * C_A * L_A * L_1) / ( 2 * M_T * ( ((L_T **2)/12) + (D_cm **2)) ) 
K_2 = (g * D_cm) / (((L_T * L_T)/12) + (D_cm * D_cm))
K_3 = (u * D_cm**2) / ( M_T * ( ((L_T **2)/12) + (D_cm **2)) )
U = (K_2 * math.sin(operationPointRad)) / (K_1 * math.cos(operationPointRad)**2)
print(U)

# # ==============================
# # Variéveis e vetores de tempo/entrada/referência
# # ==============================
periodo = 0.001
timeEnd = 10 # Tempo de execução do sistema
time = np.arange(0,timeEnd,periodo)  # Array do tempo de execução do sistema
u = U * np.ones(time.shape) # Array preenchido com a entrada para o ponto de operação que não será alterada
ref = operationPointRad * np.ones(time.shape) # Array de referência para o ponto de operação

# # ==============================
# # Aplicando degraus no ponto de operação 
# # ==============================
u_degrau = U * np.ones(len(time)) # Foi criado um novo array para os degraus 
# u_degrau = np.array_split(u_degrau, 3) # O array foi separado em 3 novos arrays
# Foi adicionado os degraus para cada intervalo
# u_degrau[0].fill(U)
u_degrau[6000:].fill(1.2 * U)
# u_degrau[1].fill(1.2 * U)
# u_degrau[2].fill(U) 
# u_degrau = np.concatenate([u_degrau[0], u_degrau[1], u_degrau[2]]) # Faço o concat dos arrays dos degraus

# # ==============================
# # Criando uma função para a equação diferencial do sistema
# # ==============================
def model_update(t,x,u,params):
    # Variáveis de estado
    x1 = x[0] # Posição angular
    x2 = x[1] # Velocidade angular
    return np.array([x2, K_1*(np.cos(x1)**2)*u - (K_2*np.sin(x1) + K_3*x2)]) # Retornando a EDO em forma de vetor
def model_output(t,x,u,params):
    # Para o caso em estudo, a saí­da do sistema, y, é o estado X[0], ou seja, y = theta.   
    return x[0]
#Definindo o sistema não linear da FanPlate, segundo a biblioteca-------------------------------------------------------------
FanPlate = ct.NonlinearIOSystem(model_update, model_output , states=2, name='FanPlate', inputs = ('u'), outputs = ('y'))

# # ==============================
# # Malha aberta
# # ==============================
# Simulação do sistema sem aplicação de degraus na entrada
t, y = ct.input_output_response(FanPlate, time, u, intialPosition)
# Plot do sistema sem aplicação de degraus na entrada
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,ref,'--k', label='ref')
plt.plot(t,y,'b', label='${x_1(t)}$')
plt.ylabel('${x_1(t)}$[rad]')
plt.xlim(0,10)
plt.ylim(0,1.25)
plt.legend()
#plt.title('Resposta temporal do sistema em malha aberta sem degrau')
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,u,'b', label='${u(t)}$')
plt.ylabel('${u(t)}$')
plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,100)
plt.ylim(2,4)
plt.grid()
plt.show()

# Simulação do sistema aplicando degraus na entrada 
t_degrau, y_degrau = ct.input_output_response(FanPlate, time, u_degrau, intialPosition)
# Plot do sistema com aplicação de degraus na entrada
plt.figure(2)
plt.subplot(2,1,1)
plt.plot(t_degrau,ref,'--k', label='ref(0.873rad)')
plt.plot(t_degrau,y_degrau,'b', label='$\\theta$(t)')
plt.ylabel('$\\theta$(t)[rad]')
plt.xlim(0,10)
plt.ylim(0.3,1.5)
plt.legend()
#plt.title('Resposta temporal do sistema em malha aberta com degrau')
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,u_degrau,'b', label='${u(t)}$')
plt.ylabel('${u(t)}[m^2/s^2]$')
plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,100)
plt.ylim(2,4)
plt.grid()
plt.show()

# Normalizando a curva do degrau positivo analisado
yn = y_degrau[6000:10000]
tn = t_degrau[0:len(yn)] 

thetan = (yn - min(yn))/(y_degrau[9999] - min(yn))

# Vetores de referência para os gráficos normalizados
refn = np.ones(tn.shape)               # Referência do ponto de equilíbrio
reftsp = (1+ 0.02)*np.ones(tn.shape)   # Referência da margem de +2%
reftsn = (1- 0.02)*np.ones(tn.shape)   # Referência da margem de -2%
reftsp5 = (1+ 0.05)*np.ones(tn.shape)  # Referência da margem de +5%
reftsn5 = (1- 0.05)*np.ones(tn.shape)  # Referência da margem de -5%


tp = 0.3440  # Instante de pico
b = 0.5883   # y(tp)
Mp = b/1     # Sobressinal máximo (overshoot)
tr = 0.1872  # Tempo de subida
ts = 2.502  # Tempo de acomodação (2%)
zeta = 0.172 # Constante de amortecimento
wn = 9.2714  # Frequência natural do sistema

Ep = 1 + (np.e*((-zeta)*wn*tn))/(np.sqrt(1-(zeta*2)))  # Curva envoltória superior
En = 1 - (np.e*((-zeta)*wn*tn))/(np.sqrt(1-(zeta*2)))  # Curva envoltória inferior

plt.figure(3)
plt.plot(tn,refn,'--k',label='ref')
plt.plot(tn,thetan,'b',label='$\\theta$(t)')
plt.plot(tn,reftsp,'--y',label='ref $\pm$ 2%')
plt.plot(tn,reftsn,'--y')
plt.plot(tn,reftsp5,'--g',label='ref $\pm$ 5%')
plt.plot(tn,reftsn5,'--g')
plt.plot(tn,Ep,'r', label= 'Curvas envoltórias')
plt.plot(tn,En,'r')
plt.ylabel('$\\theta(t)$ [rad]')
plt.xlabel('Tempo [s]')
plt.xlim(0,4)
plt.ylim(0,2)
plt.grid()
plt.show()

s=ct.tf("s") # Frequência

K=0.004302

eqTransf = K*(wn**2)/((s**2)+(2*s(wn*zeta))+(wn**2)) #Equação de transferência

timeTransfer, tempOut = ct.forced_response(eqTransf,T=time, U=u_degrau-U) 

plt.figure(4)
plt.subplot(211) 
#plt.title('Equação de Transfêrencia')
plt.plot(timeTransfer,tempOut+(np.radians(50)),'b',label='Modelo')
plt.plot(t_degrau,y_degrau,'y', label='Real')
plt.plot(time,ref,'--k',label='ref(0.873rad)')
plt.grid()
plt.ylabel('$\\theta$(t)[rad]')
plt.legend()
plt.xlim(0,10)
plt.ylim(0.3,1.5)

plt.subplot(212) 
plt.plot(timeTransfer,u_degrau,'b', label = '$Q(t)$')
plt.grid()
plt.ylabel('${u(t)}[m^2/s^2]$')
plt.xlabel('Tempo [s]')
plt.legend()
plt.xlim(0,100)
plt.ylim(0,10)
plt.show()