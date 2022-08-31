#importando as bibliotecas necessárias para fazer a programação pelo método 
# de Euler utilizando a biblioteca control systems 
from matplotlib import pyplot as plt
import numpy as np
import control as ct 
import statistics as st #Importando a biblioteca Statistics

plt.close('all') #para fechar todas as janelas de gráficos

#criando uma função para cálculo da equação diferencial do sistema
#declaracao dos parametros
def model_update(t, T, E_j, params):
  rho = params.get ('rho',1000)      #massa específica da solução (Kg/m^3)
  V = params.get('V',10)             #volume (m^3)
  c   = params.get ('c',4500)        #capacidade térmica da solução (J/KgºC)
  h   = params.get ('h',15)          #coeficiente de convecção natural (J/sm^2ºC)
  A   = params.get ('A',31.4)        #superfície do tanque para troca de calor por conveccao (m^2)
  T_amb = params.get ('T_amb',20)    #temperatura ambiente (ºC)
  Tr = params.get('Tr',50)           #referencial arbitrátrio de temperatura (ºC)
  F = params.get('F', 2.97)          #vazão volumétrica (m^3/s)
# A corrente de entrada i é dada pela variável interna u enquanto o estado do sistema 
# T é descrito pela variável interna x
  dT = (E_j/(rho*V*c))+(h*A*(T_amb-T)/(rho*V*c))+((F*(T_amb-T)/V))

  return dT

def model_output(t, T, E_j,params): # os dados retornados nessa aquisição estarão 
                                    # disponíveis para aquisição e tratamento externo
                                
  return T 

#Definindo o sistema não linear do tanque

SYSTEM = ct.NonlinearIOSystem(model_update, model_output,states = 1, name='System', inputs =('u'), outputs = ('y'))#construção do sistema não linear

#.....simulação.....
E_j0 = 868755615 #sinal de controle
tref = 85 # ponto de operação com relação a temperatura
tempo_final = 100 #duração da simulação
time = np.linspace(0,tempo_final,tempo_final+1) #simulação temporal com o espaçamento de 1s
E_j = E_j0*np.ones(time.shape) # sinal de controle necessário para levar o sistema para o ponto de operação desejado em malha aberta (teste)
ref = tref * np.ones(time.shape) # ponto de operação de referência

#aplicando os degraus
E_j[20:50]=1.05*E_j0 #degrau positivo
E_j[50:80]=0.95*E_j0 #degrau negativo
E_j[80:100]=E_j0 #volta a estabilidade

T0 = 20 #condição inicial do sistema em relação a temperatura

T,y = ct.input_output_response(SYSTEM, time, E_j, T0)# respostas do sistema não linear baseado nas entradas 

#-----------------PLOTANDO A VARIÁVEL TEMPERATURA ------------------------------------
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(T,y, ref, '--k', label='ref')
plt.plot(T,y,'b', label='T(t)')
plt.grid()
plt.ylabel('T(ºC)')
plt.legend()
plt.yticks(range(20, 100, 10))
plt.ylim(0,100)
plt.xlim(0,tempo_final)
plt.title('Resposta do Sistema de Acordo com a Modelagem Caixa Branca')

#---------------PLOTANDO A VARIÁVEL E_J --------------------------------

plt.subplot(2,1,2)
plt.plot(time, E_j, 'b',label='E_j0')
plt.grid()
plt.ylabel('E_j (J/s)')
plt.xlabel('tempo(s)')
plt.legend()
plt.ylim(0,1500000000)
plt.xlim(0,tempo_final)
plt.show()

#---MODELO DE FREQUÊNCIA ---
s=ct.tf("s")
rho =1000
c = 4500
F = 2.97
V = 10
h = 15
A = 31.4
T_amb = 20
G=1/((s+(F/V+h*A/(V*rho*c)))*(rho*c*V))
t, ysp = ct.forced_response(G,T=time, U = E_j) # Validação da parametrização pelo degrau de subida pelo método dos três parâmetros
ysp = ysp + T_amb # Somando o referencial

#-------------------PLOTANDO A FUNÇÃO DE TRANSFERÊNCIA ---------------------------
plt.figure(2)
plt.subplot(2,1,1)
plt.plot(T,y, ref, '--k', label='ref')
plt.plot(T,y,'b', label='T(t)')
plt.plot(t,ysp,'r',label='T(s)')
plt.grid()
plt.ylabel('T(ºC)')
plt.legend()
plt.yticks(range(20, 100, 10))
plt.ylim(20,100)
plt.xlim(0,tempo_final)
plt.title('Resposta em função da frequência do sistema de acordo com a modelagem em Caixa Branca')

#--------------PLOTANDO A VARIÁVEL E_J NOVAMENTE------------------------------

plt.subplot(2,1,2)
plt.plot(time, E_j, 'b',label='E_j0')
plt.grid()
plt.ylabel('E_j (J)')
plt.xlabel('tempo(s)')
plt.legend()
plt.ylim(0,1500000000)
plt.xlim(0,tempo_final)
plt.show()
#-----------------IREI USAR A MINHA PARTE DE CIMA ------------------------------
#----------------- MÉTODO DE 3 PARÂMETROS ----DEGRAU POSITIVO-------------------
K_positivo = (y[48]-85)/(E_j[48]-E_j0) #cálculo do ganho estático para o degrau positivo
coeficiente_angular = (y[21]-y[20])/(t[21]-t[20]) # Coeficiente angular da curva tangente ao degrau positivo
coeficiente_linear = y[20]-coeficiente_angular*t[20] # Coeficiente linear da reta tangente ao degrau negativo
reta = coeficiente_angular*t + coeficiente_linear # Reta tangente ao degrau de subida
tmetodo_3p= 23-17.5 # método dos três parâmetros
#------------------MÉTODO DE MILLER-------------------------------
tmmill= 23.7-17.5 # método de Miller

G_3_parametros=K_positivo/(tmetodo_3p*s+1) #dinâmica apresentada pelo método de 3 parâmetros com o degrau positivo
G_Miller=K_positivo/(tmmill*s+1) # dinâmica apresentada pelo método de Miller para o degrau positivo
t, yps = ct.forced_response(G_3_parametros,T=time, U=E_j)
yps=yps+tref
t, yms = ct.forced_response(G_Miller,T=time, U=E_j)
yms=yms+tref

#---------------PLOTANDO O DEGRAU POSITIVO------------------------
plt.figure(3)
plt.plot(t,y,'b',label='$T_0$')
plt.plot(t,reta,'r',label='$T_1$')
plt.plot(t,yps,'g',label='$T_2$')
plt.plot(t,yms,'c',label='$T_3$')
plt.plot(t,y[46]*np.ones(len(t)),'--r',label='$T_4$')
plt.plot(t,(0.63*(y[46]-y[20])+85)*np.ones(len(t)),'--m',label='$T_5$')
plt.grid()
plt.title('Degrau positivo')
plt.ylabel('$T(°C)$')
plt.legend()
plt.ylim(80,90)
plt.xlim(15,30)
plt.show()

#-----------MÉTODO DE 3 PARÂMETROS----DEGRAU NEGATIVO--------------------------

Kd = (y[79]-y[49])/(E_j[79]-E_j[49])

coeficiente_angular_1 = (y[52]-y[51])/(t[52]-t[51])# Coeficiente angular da curva tangente ao degrau de descida
coeficiente_linear_1 = y[51]-coeficiente_angular_1*t[51]# Coeficiente linear da reta tangente ao degrau de descida
reta1 = coeficiente_angular_1*t + coeficiente_linear_1 #função da reta tangente ao degrau de descida
tmetodo_3p_n=55.629-49
#---------MÉTODO DE MILLER-----------
tmmill_n=55.103222-49

Gpd=Kd/(tmetodo_3p_n*s+1)#dinâmica apresentada pelo método de 3 parâmetros com o degrau negativo
Gmd=Kd/(tmmill_n*s+1)#dinâmica apresentada pelo método de Miller para o degrau positivo
t, ypd = ct.forced_response(Gpd,T=time, U=E_j)
ypd=ypd+tref
t, ymd = ct.forced_response(Gmd,T=time, U=E_j)
ymd=ymd+tref

#---------------PLOTANDO O DEGRAU POSITIVO------------------------

plt.figure(4)
plt.plot(t,y,'b',label='$T_0$')
plt.plot(t,reta1,'r',label='$T_1$')
plt.plot(t,ypd,'c',label='$T_2$')
plt.plot(t,ymd,'g',label='$T_3$')
plt.plot(t,y[79]*np.ones(len(t)),'--r',label='$T_4$')
plt.plot(t,(0.63*(y[79]-y[49])+y[49])*np.ones(len(t)),'--m',label='$T_5$')
plt.grid()
plt.title('Degrau negativo')
plt.ylabel('$T (°C)$')
plt.legend()
plt.ylim(80,90)
plt.xlim(45,81)
plt.show()

# ------------------ COMPARANDO AS RESPOSTAS DO SISTEMA PELO ÍNDICE DE PERFORMANCE RMSE ------------------
RMSE = np.zeros(4)

instante_dos_degraus = y[20:80] # T no instante dos degrais

erro_degrau_subida_3_p = yps[20:80] # Erro da curva parametrizada pelo degrau de subida pelo método dos três parâmetros
RMSE[0] = np.sqrt(st.mean(((erro_degrau_subida_3_p-instante_dos_degraus)**2)/100000))

Erro_degrau_descida_3_p = ypd[20:80] # Erro da curva parametrizada pelo degrau de descida pelo método dos três parâmetros
RMSE[1] = np.sqrt(st.mean(((Erro_degrau_descida_3_p-instante_dos_degraus)**2)/100000))

erro_degrau_subida_miller = yms[20:80] # Erro da curva parametrizada pelo degrau de subida pelo método de Miller
RMSE[2] = np.sqrt(st.mean(((erro_degrau_subida_miller-instante_dos_degraus)**2)/100000))

erro_degrau_descida_miller = ymd[20:80] # Erro da curva parametrizada pelo degrau de descida pelo método de Miller
RMSE[3] = np.sqrt(st.mean(((erro_degrau_descida_miller-instante_dos_degraus)**2)/100000))

print('O RMSE DO DEGRAU DE SUBIDA PELO MÉTODO DE 3 PARÂMETROS É', RMSE[0])
print('O RMSE DO DEGRAU DE DESCIDA PELO MÉTODO DE 3 PARÂMETROS É', RMSE[1])
print('O RMSE DO DEGRAU DE SUBIDA PELO MÉTODO DE MILLER É', RMSE[2])
print('O RMSE DO DEGRAU DE DESCIDA PELO MÉTODO DE MILLER É', RMSE[3])

# -------------------------------------NORMALIZAÇÃO ---------------------------
yn=(y[0:20]-y[0])/(y[20]-y[0])
yl=(ysp[0:20]-ysp[0])/(ysp[20]-ysp[0])
yn=yn+1
# ------------------------- PLOTANDO A NORMALIZAÇÃO----------------------------
plt.figure(7)
plt.title('Normalização') 
plt.plot(t[0:20],yn,'b',label='$T_0$')
plt.plot(t[0:20],yl,'r',label='$T_1$')
plt.plot(t[0:20],0.98*yn[-1]*np.ones(len(t[0:20])),'--r',label='$T_2$') #tempo de acomodação (tempo aonde se inicia o degrau até alcançar 98% da acomodação.
plt.plot(t[0:20],0.9*yn[-1]*np.ones(len(t[0:20])),'--k',label='$T_3$') # tempo de subida  \ (90% da subida - 10% da subida)
plt.plot(t[0:20],0.1*yn[-1]*np.ones(len(t[0:20])),'--k',label='$T_4$') # tempo de subida  /
plt.grid()
plt.ylabel('$T(°C)$')
plt.legend()
plt.ylim(0,1.1)
plt.xlim(0,18)
plt.show()

#--------MÉTODO DA BISSEÇÃO------------
x=((41/K_positivo-1)/y)-1

a=1
b=1909992109
e=0.01

def f(y):
    return ((41/K_positivo-1)/y)-1

if f(a)*f(b)<0:
    while np.abs((b-a)/2)>e:
        xi=(a+b)/2
        if f(xi)==0:
            print ("A raiz é:", xi)
            break
        else:
            if f(a)*f(xi)<0:
                b=xi
            else: 
                a=xi             
    print ("O valor da raiz é:",xi)
else:
    print ("Não há raiz.")