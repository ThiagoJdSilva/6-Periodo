#@Autor: Thiago José da Silva e Luiza Gomes de Castro e Sá
#@Data: 21/11/2021

import numpy as np #Importando a biblioteca numpy
from matplotlib import pyplot as plt #Importando a biblioteca matplotlib responsável por plotar gráficos
import control as ct #Importando a biblioteca control para realizar o controle do sistema

plt.close('all') #Fechando todas as abas de gráficos abertas


def model_uptade(t,T,Q,params):
    """
    Função resposável por 
    """
    V = params.get('V', 10)
    C = params.get('C', 4500)
    F = params.get('F', 3)
    T_inicial = params.get('T_inicial', 20)
    h = params.get('h', 15)
    A_ext = params.get('A_ext', 31.4)
    rho = params.get('rho', 1000)

    dT = (Q/(rho*V*C))+(h*A_ext*(T_inicial-T)/(rho*V*C))+((F*(T_inicial-T)/V))
    return dT


def model_output(t,T,Q,params):
    return T

V =  10
C = 4500
F = 3
T_inicial = 20
h = 15
A_ext = 31.4
rho = 1000

SYSTEM = ct.NonlinearIOSystem(model_uptade, model_output,
states=1, name='SYSTEM', inputs=('u'), outputs=('y'))

tempo_final = 30
time = np.linspace(0, tempo_final, tempo_final+1)
Q0 = 540018840 * np.ones(time.shape)
ref = 60 * np.ones(time.shape)
T0 = 20

T,y = ct.input_output_response(SYSTEM, time, Q0, T0, method='Radau')

#Gráfico EDO
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(T,ref,'--k', label='ref')
plt.plot(T,y,'b', label='T(t)')
plt.ylabel('T(t)[C]')
plt.xlim(0,tempo_final)
plt.legend()
plt.ylim(0,100)
plt.title('Resposta temporal de tanques acoplados em malha aberta')
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,Q0,'b', label='Q0')
plt.ylabel('Q(t)')
plt.legend()
plt.xlabel('Tempo[s]')
plt.xlim(0,tempo_final)
plt.ylim(0,1000000000)
plt.grid()
plt.show()

"""
Plotando gráfico da função de transferência
"""
s=ct.tf("s")

transformacao_de_laplace = 1/((s+(F/V+h*A_ext/(V*rho*C)))*(rho*C*V))

t, ysp = ct.forced_response(transformacao_de_laplace,T=time, U=Q-Q0) # Validação da parametrização pelo degrau de subida pelo método dos três parâmetros

ysp += T_inicial # Somando o referencial

#Gráfica função de transferência
plt.figure(2)
plt.subplot(211) 
plt.title('Equação de Transfêrencia')
plt.plot(t,ysp,'b',label='$\\theta_0$')
plt.plot(time,ref,'--k',label='ref')
plt.grid()
plt.ylabel('$\\theta (°)$')
plt.legend()
plt.ylim(15,55)
plt.xlim(0,tempo_final)
plt.subplot(212) 
plt.plot(t,Q0,'b', label = '$Q1$')
plt.grid()
plt.ylabel('Q (J)')
plt.xlabel('Tempo (s)')
plt.legend()
plt.xlim(0,tempo_final)
plt.ylim(0,1000000000)
plt.show()

#Definindo os degrus
CalorDegraus = np.array_split(CalorDegraus, 4)

CalorDegraus[0][:] = Q0[0]
CalorDegraus[1][:] = 1.15*Q0[0]
CalorDegraus[2][:] = Q0[0]
CalorDegraus[3][:] = 0.85*Q0[0]

CalorDegraus = np.concatenate([CalorDegraus[0], CalorDegraus[1], CalorDegraus[2], CalorDegraus[3]])
valor_alocado[0] = T_inicial   
    
CalorDegraus = np.concatenate([CalorDegraus[0], CalorDegraus[1], CalorDegraus[2], CalorDegraus[3]]) # Faço o concat dos arrays dos degraus
valor_alocado[0] = T_inicial # Iniciar com um valor previo, pos um Not a Number não pode ser somado na integração abaixo

# Integrando pelo metodo de Taylor
for k in range(len(time)-1): 
    valor_alocado[k+1] = valor_alocado[k]+((F*(T_inicial - valor_alocado[k])/V)+ (CalorDegraus[k]/(rho*V*C)) + ((h*A_ext)*(T_inicial - valor_alocado[k])/(rho*V*C)))

plt.figure(3)
plt.subplot(211)
plt.title('Equação Edo')
plt.plot(time,valor_alocado,'b',label='$\\t$')
plt.plot(time,ref,'--k',label='ref')
plt.grid()
plt.ylabel('$\\t (°)$')
plt.legend()
plt.ylim(15,55)
plt.xlim(0,timeEnd)

plt.subplot(212)
plt.plot(time,CalorDegraus,'b',label='Q1')
plt.grid()
plt.ylabel('Q (J) ')
plt.xlabel('Tempo (s)')
plt.legend()
plt.ylim(0,1000000000)
plt.xlim(0,timeEnd)
plt.show()