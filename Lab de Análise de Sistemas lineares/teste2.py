#Relatório 01 - Laboratório de Análise de Sistema Lineares
#Modelagem de Aproximada de Processos Industriais-Modelagem Caixa Branca e Caixa Preta
#Matheus Alvim Santos Teixeira
#Weverton Carlos Souza

import numpy as np              #importando a biblioteca Numpy que trata as operações vetoriais
import matplotlib.pyplot as plt #importando a biblioteca que gera os gráficos.
import control as ct            #importando biblioteca
import statistics as sta        #importando a biblioteca Statistics 

plt.close('all') #fechando todas as janelas de gráficos

# Declaração das Constantes
V = 10        #volume do reator
C = 4500      #capacidade térmica da solução
F = 2         #vazao volumetrica
thetai = 20   #temperatura ambiente
h = 15        #coeficiente de convecção
Aext = 31.4   #superficie do tanque   
rho = 1000    #massa especifica da solução

theta_operacao = 50 #ponto de operação desejado
Q0 = -(rho*C*F*thetai)+(rho*C*F*theta_operacao)+(h*Aext*(thetai-theta_operacao)) #equação que calcula 
# a quantidade necessária para levar o sistema ao ponto de operação
#pre alocando vetores
T  = 0.1   #periodo de amostragem ou intervalo de integração (s)
tf = 300   #tempo final da simulação
vetor_tempo  = np.arange(0,tf,T) #criando vetor tempo, com passo de 1s
theta = np.empty(len(vetor_tempo)) #alocando o vetor para a temperatura do liquido
theta.fill(np.nan) #limpando o vetor e definindo-o como NAN
Q1 = Q0*np.ones(len(vetor_tempo)) #criando o vetor para a quantidade de calor inicial
ref = theta_operacao*np.ones(len(vetor_tempo)) #criando o vetor para o ponto de referencia

#Configurando a condição inicial do sistema
theta[0] = 20 # definindo a condição inicial do sistema, neste caso a temperatura inicial é a temperatura ambiente

#integração da equação diferencial
for k in range(len(vetor_tempo)-1): 
    theta[k+1] = theta[k]+T*((F*(thetai - theta[k])/V)+ (Q1[k]/(rho*V*C)) + ((h*Aext)*(thetai - theta[k])/(rho*V*C)))
    
#Plotando o gráfico da equação diferencial
plt.figure(1)
plt.subplot(211)
plt.title('Curva da EDO')
plt.plot(vetor_tempo,theta,'b',label='$\\theta$(t)')
plt.plot(vetor_tempo,ref,'--k',label='ref')
plt.grid()
plt.ylabel('$\\theta(t)  [°C]$')
plt.legend()
plt.ylim(15,55)
plt.xlim(0,35)
#Plotando gráfico do sinal de entrada
plt.subplot(212)
plt.plot(vetor_tempo,Q1,'b',label='$Q(t)$')
plt.grid()
plt.ylabel('Q(t) [J] ')
plt.xlabel('Tempo [s]')
plt.legend()
plt.ylim(0,1000000000)
plt.xlim(0,35)
plt.show()


#Função de transferência
s=ct.tf("s")
G=1/((s + (F/V+h*Aext/(V*rho*C)))*(rho*C*V)) #equação de transferência
t, ysp = ct.forced_response(G,T=vetor_tempo, U=Q1)
ysp = ysp + thetai #somando a temperatura ambiente

#plotando o gráfico da equação de transferencia e comparando com a curva da EDO
plt.figure(2)
plt.title('Comparação entre a Curva da EDO e a Curva da Função de Transfência')
plt.plot(vetor_tempo,theta,'b',label='$\\theta_1$(t)') 
plt.plot(t,ysp,'r',label='$\\theta_2$(t)')
plt.plot(vetor_tempo,ref,'--k',label='ref')
plt.grid()
plt.ylabel('$\\theta(t) [°C]$')
plt.legend()
plt.ylim(15,55)
plt.xlim(0,33)
plt.xlabel('Tempo [s]')

x0 = 9
xf = 9.8
y0 = 44.8
yf = 45.6



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
plt.plot(vetor_tempo,theta,'b',label='$\\theta_1$(t)') 
plt.plot(t,ysp,'r',label='$\\theta_2$(t)')
plt.plot(vetor_tempo,ref,'--k',label='ref')
plt.show()



#Inserindo degraus na entrada
Q2 = Q0*np.ones(len(vetor_tempo)) #criando o vetor para a quantidade de calor inicial

Q2[0:400]= Q0
Q2[400:800]= 1.10*Q0
Q2[800:1200]= Q0
Q2[1200:1600]= 0.90*Q0
Q2[1600:2000]= Q0
Q2[2000:2400]= 1.10*Q0
Q2[2400:2800]= Q0

for k in range(len(vetor_tempo)-1): 
    theta[k+1] = theta[k]+T*((F*(thetai - theta[k])/V)+ (Q2[k]/(rho*V*C)) + ((h*Aext)*(thetai - theta[k])/(rho*V*C)))

#Plotando as variáveis de interesse: theta e Q2 da equação diferencial
plt.figure(3)
plt.subplot(211)
plt.title('Comportamento da Curva com Degraus Positivos e Negativo')
plt.plot(vetor_tempo,theta,'b',label='$\\theta(t)$')
plt.plot(vetor_tempo,ref,'--k',label='ref')
plt.grid()
plt.ylabel('$\\theta(t) [°C]$')
plt.legend()
plt.ylim(0,100)
plt.xlim(0,tf)

x0 = 35
xf = 180
y0 = 45
yf = 55


plt.plot([x0,xf], [y0,y0], 'r--')
plt.plot([x0,xf], [yf,yf], 'r--')
plt.plot([x0,x0], [y0,yf], 'r--')
plt.plot([xf,xf], [y0,yf], 'r--')

plt.subplot(212)
plt.plot(vetor_tempo,Q2,'b',label='Q(t)')
plt.grid()
plt.ylabel('Q(J) ')
plt.xlabel('Tempo(s)')
plt.legend()
plt.ylim(100000000,1000000000)
plt.xlim(0,tf)

# ==============================
# PLOT NORMAL
# ==============================

plt.plot([x0,xf], [y0,y0], 'r--')
plt.plot([x0,xf], [yf,yf], 'r--')
plt.plot([x0,x0], [y0,yf], 'r--')
plt.plot([xf,xf], [y0,yf], 'r--')

# ==============================
# PLOT COM ZOOM
# ==============================

a = plt.axes([0.35, 0.77, 0.25, 0.10]) # left, right, width, height
plt.xlim(x0,xf)
plt.ylim(y0,yf)
plt.grid()
plt.plot(vetor_tempo,theta,'b',label='$\\theta(t)$')
plt.plot(vetor_tempo,ref,'--k',label='ref')
plt.show()



#Plotando o degrau positivo para analisar os parâmetros do método de 3 parâmetros
#Degrau positivo
m1 = (theta[402]-theta[401])/(vetor_tempo[402]- vetor_tempo[401]) # Coeficiente angular da curva tangente ao degrau de subida
eq1 = m1*vetor_tempo - m1*vetor_tempo[401] + theta[401] #equação da reta tangente
ref1 = 52.995596810296654*np.ones(len(vetor_tempo)) #referencia da resposta ao degrau 
ref2 = (0.63*(theta[799] - theta[400])+50)*np.ones(len(vetor_tempo)) #referencia 63% do gamho de Kpositivo, método de Miller

plt.figure(4) #plotando o degrau positivo
plt.title('Degrau Positivo aqui')
plt.plot(vetor_tempo,theta,'b',label='resposta Edo no degrau positivo')
plt.plot(vetor_tempo,ref,'--k',label='ref ponto de operação')
plt.plot(vetor_tempo,eq1,'r',label='reta tangente')
plt.plot(vetor_tempo,ref1,'--g',label='ref da resposta no degrau')
plt.plot(vetor_tempo,ref2,'--y',label='ref de 63% K positivo')
plt.grid()
plt.ylabel('$\\theta(t) [°C]$')
plt.xlabel('Tempo(s)')
plt.legend()
plt.ylim(49,54)
plt.xlim(38,70)
plt.show()

#Degrau Negativo
m2 = (theta[1202]-theta[1201])/(vetor_tempo[1202]- vetor_tempo[1201]) # Coeficiente angular da curva tangente ao degrau de subida
eq2 = m2*vetor_tempo - m2*vetor_tempo[1201] + theta[1201] #equação da reta tangente
ref3 = 46.998120884452426*np.ones(len(vetor_tempo)) 
ref4 = (0.63*(theta[1599] - theta[1200])+50)*np.ones(len(vetor_tempo)) #referencia 63% do gamho de Knegativo, método de Miller

plt.figure(5)
plt.title('Degrau Negativo')
plt.plot(vetor_tempo,theta,'b',label='Resposta EDO no degrau negativo')
plt.plot(vetor_tempo,ref,'--k',label='ref ponto de operação')
plt.plot(vetor_tempo,eq2,'r',label='reta tangente')
plt.plot(vetor_tempo,ref3,'--g',label='ref da resposta no degrau')
plt.plot(vetor_tempo,ref4,'--y',label='ref de 63% K negativo')
plt.grid()
plt.ylabel('$\\theta(t) [°C]$')
plt.xlabel('Tempo(s)')
plt.legend()
plt.ylim(46,51)
plt.xlim(118,150)
plt.show()


#Calculando a FT para o degrau positivo utilizando o segmento AC
Kp = (theta[799]-theta[400])/(Q2[799]-Q0)  #Calculando K método Ziegler Nichols p/ o degrau positivo
tau1 = 44.975-40 #calculando tau, segmento AC.
#como não existe atraso, a variável theta é igual a zero
Gp = Kp/(tau1*s+1) #modelo do degrau positivo
t, yp = ct.forced_response(Gp,T= vetor_tempo, U = Q2) #trazendo a equação de TF para o domínio do tempo
yp = yp + thetai #somando a temperatura inicial

#Calculando a FT para o degrau positivo utilizando o segmento AB, Método de Miller
tau2 = 3.231 #Segmento AB
Gpm = Kp/(tau2*s+1) #modelo do degrau positivo Miller
t, ypm = ct.forced_response(Gpm,T= vetor_tempo, U = Q2)
ypm = ypm + thetai

#Calculando a FT para o degrau negativo utilizando o segmento AC
Kn = (theta[1599]-theta[1200])/(Q2[1599]-Q0)  #Calculando K do método Ziegler Nichols p/ o degrau negativo
tau3 = 125.095-120 #calculando tau, segmento AC.
#como não existe atraso, a variável, theta é igual a zero
Gn = Kn/(tau1*s+1) #modelo do degrau negativo
t, yn= ct.forced_response(Gn,T= vetor_tempo, U = Q2) #trazendo a equação de TF para o domínio do tempo
yn = yn + thetai #somando a temperatura inicial

#Calculando a FT para o degrau negativo utilizando o segmento AB, Método de Miller
tau4 = 3.21 #SegmentoAB
Gnm = Kn/(tau4*s+1) #modelo do degrau negativo Miller
t, ynm = ct.forced_response(Gnm,T= vetor_tempo, U = Q2)
ynm = ynm + thetai

# Plotando todos os modelos - EDO, Método Ziegler Nichols Degrau Positivo e Negativo
# e Méodo de Miller Degrau positivo e negativo com degraus diferentes do que foram
# utilizados como base para a confecção dos modelos

Q3 = Q0*np.ones(len(vetor_tempo)) #criando o vetor para a quantidade de calor inicial
 
Q3[0:400]= Q0
Q3[400:800]= 1.13*Q0
Q3[800:1200]= Q0
Q3[1200:1600]= 0.87*Q0
Q3[1600:2000]= 1.13*Q0
Q3[2000:2400]= 0.87*Q0
Q3[2400:2800]= Q0

for k in range(len(vetor_tempo)-1): 
    theta[k+1] = theta[k]+T*((F*(thetai - theta[k])/V)+ (Q3[k]/(rho*V*C)) + ((h*Aext)*(thetai - theta[k])/(rho*V*C)))

t, yp = ct.forced_response(Gp,T= vetor_tempo, U = Q3)
yp = yp + thetai

t, ypm = ct.forced_response(Gpm,T= vetor_tempo, U = Q3)
ypm = ypm + thetai

t, yn = ct.forced_response(Gn,T= vetor_tempo, U = Q3)
yn = yn + thetai

t, ynm = ct.forced_response(Gnm,T= vetor_tempo, U = Q3)
ynm = ynm + thetai


plt.figure(6)
plt.title('Comparação dos Modelos')
plt.plot(vetor_tempo,theta,'b',label='Sistema Real') #EDO
plt.plot(t,yp,'g',label='Modelo 1 ') #Modelo degrau positivo método Ziegler Nichols
plt.plot(t,ypm,'r',label='Modelo 2') #Modelo degrau positivo método de Miller
plt.plot(t,yn,'c',label='Modelo 3')  #Modelo degrau negativo método Ziegler Nichols
plt.plot(t,ynm,'y',label='Modelo 4') #Modelo degrau negativo método de Miller
plt.plot(vetor_tempo,ref,'--k',label='ref')
plt.grid()
plt.ylabel('$\\theta(t) [°C]$')
plt.xlabel('Tempo(s)')
plt.legend()
plt.ylim(30,70)
plt.xlim(20,175)

x0 = 39.5
xf = 42
y0 = 49.5
yf = 52.5

x01 = 119.5
xf1 = 122
y01 = 48.5
yf1 = 51.5


# ==============================
# PLOT NORMAL
# ==============================

plt.plot([x0,xf], [y0,y0], 'b--')
plt.plot([x0,xf], [yf,yf], 'b--')
plt.plot([x0,x0], [y0,yf], 'b--')
plt.plot([xf,xf], [y0,yf], 'b--')

plt.plot([x01,xf1], [y01,y01], 'm--')
plt.plot([x01,xf1], [yf1,yf1], 'm--')
plt.plot([x01,x01], [y01,yf1], 'm--')
plt.plot([xf1,xf1], [y01,yf1], 'm--')


# ==============================
# PLOT COM ZOOM
# ==============================

a = plt.axes([0.25, 0.62, 0.20, 0.20]) # left, right, width, height
plt.xlim(x0,xf)
plt.ylim(y0,yf)
plt.plot(t,yp,'g',label='Modelo 1 ') #Modelo degrau positivo método Ziegler Nichols
plt.plot(t,ypm,'r',label='Modelo 2') #Modelo degrau positivo método de Miller
plt.plot(t,yn,'c',label='Modelo 3')  #Modelo degrau negativo método Ziegler Nichols
plt.plot(t,ynm,'y',label='Modelo 4') #Modelo degrau negativo método de Miller
plt.grid()
plt.show()

# PLOT COM ZOOM
# ==============================

a = plt.axes([0.50, 0.20, 0.20, 0.20]) # left, right, width, height
plt.xlim(x01,xf1)
plt.ylim(y01,yf1)
plt.plot(t,yp,'g',label='Modelo 1 ') #Modelo degrau positivo métodos do Ziegler Nichols
plt.plot(t,ypm,'r',label='Modelo 2') #Modelo degrau positivo método de Miller
plt.plot(t,yn,'c',label='Modelo 3 ') ##Modelo degrau negativo métodos do Ziegler Nichols
plt.plot(t,ynm,'y',label='Modelo 4') #Modelo degrau negativo método de Miller
plt.grid()



# Comparando os sistemas pelo índice RMSE
RMSE = np.zeros(12)

real = theta[400:800] #theta no intervalo dos degrais

Gp = yp[400:800] # Erro da curva parametrizada pelo degrau positivo pelo método Ziegler Nichols
RMSE[0] = np.sqrt(sta.mean(((Gp-real)**2)))

Gpm = ypm[400:800]  #Erro da curva parametrizada pelo degrau de positivo pelo método de Miller
RMSE[1] = np.sqrt(sta.mean(((Gpm-real)**2)))

Gn = yn[400:800] #Erro da curva parametrizada pelo degrau de negativo pelo método de Ziegler Nichols
RMSE[2] = np.sqrt(sta.mean(((Gn-real)**2)))

Gnm = ynm[400:800] # Erro da curva parametrizada pelo degrau negativo pelo método de Miller
RMSE[3] = np.sqrt(sta.mean(((Gnm-real)**2)))

real = theta[1200:1600] #theta no intervalo dos degrais

Gp = yp[1200:1600] # Erro da curva parametrizada pelo degrau positivo pelo método do Ziegler Nichols
RMSE[4] = np.sqrt(sta.mean(((Gp-real)**2)))

Gpm = ypm[1200:1600]  #Erro da curva parametrizada pelo degrau de positivo pelo método de Miller
RMSE[5] = np.sqrt(sta.mean(((Gpm-real)**2)))

Gn = yn[1200:1600] #Erro da curva parametrizada pelo degrau de negativo pelo método do Ziegler Nichols
RMSE[6] = np.sqrt(sta.mean(((Gn-real)**2)))

Gnm = ynm[1200:1600] # Erro da curva parametrizada pelo degrau negativo pelo método de Miller
RMSE[7] = np.sqrt(sta.mean(((Gnm-real)**2)))

real = theta[0:3000] #theta no intervalo dos degrais

Gp = yp[0:3000] # Erro da curva parametrizada pelo degrau positivo pelo método do Ziegler Nichols
RMSE[8] = np.sqrt(sta.mean(((Gp-real)**2)))

Gpm = ypm[0:3000]  #Erro da curva parametrizada pelo degrau de positivo pelo método de Miller
RMSE[9] = np.sqrt(sta.mean(((Gpm-real)**2)))

Gn = yn[0:3000] #Erro da curva parametrizada pelo degrau de negativo pelo método do Ziegler Nichols
RMSE[10] = np.sqrt(sta.mean(((Gn-real)**2)))

Gnm = ynm[0:3000] # Erro da curva parametrizada pelo degrau negativo pelo método de Miller
RMSE[11] = np.sqrt(sta.mean(((Gnm-real)**2)))

#Print dos índices 
#Degrau positivo
print('Degrau Positivo\n')
print('Índice RMSE para o modelo Ziegler Nichols degrau positivo\n', RMSE[0])
print('Índice RMSE para o modelo Miller degrau positivo\n', RMSE[1])
print('Índice RMSE para o modelo Ziegler Nichols degrau negativo\n', RMSE[2])
print('Índice RMSE para o modelo Miller degrau negativo \n', RMSE[3])
#Degrau negativo
print('---------------------------------------------------------------------')
print('Degrau Negativo\n')
print('Índice RMSE para o modelo Ziegler Nichols degrau positivo\n', RMSE[4])
print('Índice RMSE para o modelo Miller degrau positivo\n' , RMSE[5])
print('Índice RMSE para o modelo Ziegler Nichols degrau negativo\n', RMSE[6])
print('Índice RMSE para o modelo Miller degrau negativo\n', RMSE[7])
#Todo o tempo de simulação
print('---------------------------------------------------------------------')
print('Para todo tempo de simualação\n')
print('Índice RMSE para o modelo Ziegler Nichols degrau positivo\n', RMSE[8])
print('Índice RMSE para o modelo Miller degrau positivo\n' , RMSE[9])
print('Índice RMSE para o modelo Ziegler Nichols degrau negativo\n', RMSE[10])
print('Índice RMSE para o modelo Miller degrau negativo\n', RMSE[11])