#Laboratório de Análise de Sistemas Lineares
#@Autores: Thiago José da Silva e Luiza Gomes de Castro e Sá
#@Data: 18/01/2022

import numpy as np                          #Importando a biblioteca numpy
from matplotlib import pyplot as plt        #Importando a biblioteca matplotlib responsável por plotar gráficos
import control as ct                        #Importando a biblioteca control para realizar o controle do sistema
import statistics as sta                    #Importando a biblioteca statistic
import math

plt.close('all') #Fechando todas as abas de gráficos aberta
#----------------------------------------------------------------------------------------------------------------------------------------
" Inicio Guia 5"
#----------------------------------------------------------------------------------------------------------------------------------------

def model_update(t,x,u_array,params):
    u = u_array[0]
    u_massa = u_array[1]
    MT = m_t*(1 + u_massa)

    K_1 = (d_cm*rho_ar*C_a*L_a*L_1)/(2*MT*(((L_t**2)/12)+d_cm**2))
    K_2 = (g*d_cm)/ (((L_t**2)/12)+d_cm**2)
    K_3 = (atrito*(d_cm**2))/(MT*(((L_t**2)/12)+d_cm**2))

    # Variáveis de estado
    x1 = x[0] # Posição angular
    x2 = x[1] # Velocidade angular
    return np.array([x2, K_1*(np.cos(x1)**2)*u - (K_2*np.sin(x1) + K_3*x2)], dtype=object) # Retornando a EDO em forma de vetor


def model_output(t,x,u,params):
    # Para o caso em estudo, a saí­da do sistema, y, é o estado X[0], ou seja, y = theta.   
    return x[0]


def K_ganho(massa):
    K_1 = (d_cm*rho_ar*C_a*L_a*L_1)/(2*massa*(((L_t**2)/12)+d_cm**2))
    K_2 = (g*d_cm)/ (((L_t**2)/12)+d_cm**2)
    K_3 = (atrito*(d_cm**2))/(massa*(((L_t**2)/12)+d_cm**2))

    return K_1, K_2, K_3

def coeficiente_angular(valor_y, valor_x, j):
    """
    Função responsável por retornar o valor do coeficiente angular da reta
    """
    return (valor_y[j]-valor_y[j-1])/(valor_x[j]-valor_x[j-1])

def equacao(vetor1, vetor2, c_angular, k):
    """
    Retorna a equaçao da reta
    """
    return c_angular*vetor1 - c_angular*vetor1[k] + vetor2[k]

def valor_K(vetor1, vetor2, ponto_de_operacao, variavel_de_controle, k):
    return vetor1[k]-ponto_de_operacao, vetor2[k]-variavel_de_controle, (vetor1[k]-ponto_de_operacao)/(vetor2[k]-variavel_de_controle)

#----------------------------------------------------------------------------------------------------------------------------------------
" Definindo as varíaveis do sistema"
#----------------------------------------------------------------------------------------------------------------------------------------
L_a = 0.154       # Largura da placa móvel de alumínio
L_1 = 0.155       # Comprimento da placa abaixo do eixo de rotação.
L_t = 0.270       # Comprimento da placa abaixo do eixo de rotação.
d_cm = 0.020       # Distância do centro de massa da placa
rho_ar = 1.23   # Densidade do ar
m_t = 0.005
C_a = 2.05      # Coeficiente de arrasto
atrito = 5      # Coeficiented e atrito viscoso
g = 9.81        # Gravidade
ponto_de_operacao_grau = 38
ponto_radiano = math.radians(ponto_de_operacao_grau)
ponto_inicial = 0

#----------------------------------------------------------------------------------------------------------------------------------------
" Definindo as constantes K do sistema"
#----------------------------------------------------------------------------------------------------------------------------------------

K_1, K_2, K_3 = K_ganho(m_t)
# print(f"K1 = {K_1}, K2 = {K_2}, K3 = {K_3}, ")

tempo_final = 50   # Tempo final da análise
periodo = 0.001
time = np.arange(0, tempo_final, periodo) # Criando um lista que vai de 0 até 50000, com período 0.001
U_entrada = (K_2 * math.sin(ponto_radiano))/(K_1 * (math.cos(ponto_radiano)**2))
U0 = U_entrada * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao valor de entrada
ref = ponto_radiano * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao ponto de operação
u_massa = np.zeros(len(time))

print(U_entrada)
print(np.radians(3.64))
# print(f'Sinal de entrada = {U_entrada}')
#----------------------------------------------------------------------------------------------------------------------------------------
" Definindo os degraus do sistema"
#----------------------------------------------------------------------------------------------------------------------------------------
U_Degraus = U_entrada * np.ones(len(time))           #Criando um novo vetor para alocar os degraus
U_Degraus = np.array_split(U_Degraus, 4)      # Dividindo o vetor em 4 partes iguais

U_Degraus[0][:]=U0[0]                            # Definindo a primeira parte do vetor como U0
U_Degraus[1][:]= 1.20*U0[0]                       # Definindo a segunda parte do vetor como 110% de U0
U_Degraus[2][:]= U0[0]                           # Definindo a terceira parte do vetor como U0   
U_Degraus[3][:]= 0.80*U0[0]                           # Definindo a terceira parte do vetor como U0   
"""Com as linhas acima é definido os degraus"""

U_Degraus_conc = np.concatenate([U_Degraus[0], U_Degraus[1], U_Degraus[2], U_Degraus[3]]) # Unindo os vetores que haviam sido divididos

FanPlate = ct.NonlinearIOSystem(model_update, model_output , states=2, name='FanPlate', inputs = ('u','u_massa'), outputs = ('y'))

#----------------------------------------------------------------------------------------------------------------------------------------
" Sistema em Malha aberta:"
#----------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------
" Simulação do sistema sem aplicação de degraus na entrada:"
#----------------------------------------------------------------------------------------------------------------------------------------
t, y = ct.input_output_response(FanPlate, time, [U0, u_massa], ponto_inicial)
# print(y)

#----------------------------------------------------------------------------------------------------------------------------------------
" Plotando sistema sem aplicação de degraus na entrada:"
#----------------------------------------------------------------------------------------------------------------------------------------
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,ref,'--k', label='ref')
plt.plot(t,y[0],'b', label='$\\theta(t)$')
plt.ylabel('$\\theta(t)[rad]$')
plt.xlim(0,10)
plt.ylim(0,1.5)
plt.legend()
plt.title('Resposta temporal do sistema em malha aberta sem degrau')
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,U0,'b', label='U (t)')
plt.ylabel('U(t)')
plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,50)
plt.ylim(0,2.5)
plt.grid()
plt.show()

#----------------------------------------------------------------------------------------------------------------------------------------
" Plotando sistema com aplicação de degraus na entrada:"
#----------------------------------------------------------------------------------------------------------------------------------------
t_degrau, y_degrau = ct.input_output_response(FanPlate, time, [U_Degraus_conc, u_massa], ponto_inicial)
y_max_deg = y_degrau[0][24499]
# print(y_max_deg)
plt.figure(2)
plt.subplot(2,1,1)
plt.plot(t_degrau,ref,'--k', label='ref')
plt.plot(t_degrau,y_degrau[0],'b', label='$\\theta(t)$')
plt.ylabel('$\\theta(t)[rad]$')
plt.xlim(0,50)
plt.ylim(0,1.5)
plt.legend()
plt.title('Resposta temporal do sistema em malha aberta com degrau')
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,U_Degraus_conc,'b', label='Q (t)')
plt.ylabel('U(t)')
plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,50)
plt.ylim(0,2.5)
plt.grid()
plt.show()

#----------------------------------------------------------------------------------------------------------------------------------------
"Determinando os parâmetros de Ziegler-Nichols" 
#----------------------------------------------------------------------------------------------------------------------------------------
s = ct.TransferFunction.s # Cria um sistema de função de transferência.

"""Reta tangente ao degrau positivo"""
amplitude_degrau_pos = y_degrau[0][24999]                                                           # Definindo a amplitude do degrau positivo
coeficiente_angular_degrau_positivo = coeficiente_angular(y_degrau[0], time, 12600)                 # Definindo o coeficiente angular pela função
equacao_degrau_positivo = equacao(time, y_degrau[0], coeficiente_angular_degrau_positivo, 12600)    # Definindo a eqaução da reta
ref2 = amplitude_degrau_pos * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao ponto de operação
"""Determinando os parâmetros de Ziegler através da análise gráfica para o degrau positivo, o valor de theta é de 0.011, o que pode ser considerado aproximadamente nulo:"""
tau_positivo = 13.6383 - 12.5116
var_theta, var_u, K1 = valor_K(y_degrau[0], U_Degraus_conc, ponto_radiano, U_entrada, 24999)
# print(f'K = {K1}, Variaçao Saída = {var_theta}, Variação U = {var_u}, tau = {tau_positivo}')

#----------------------------------------------------------------------------------------------------------------------------------------
"Modelo de primeira ordem com os parâmetros de Ziegler-Nichols para o degrau positivo" 
#----------------------------------------------------------------------------------------------------------------------------------------

modelo_degrau_positivo = K1/(tau_positivo*s+1)                                                        # Modelo do degrau positivo
# print(modelo_degrau_positivo)
t_degrau_positivo, y_degrau_positivo = ct.forced_response(modelo_degrau_positivo,T= time, U = U_Degraus_conc - U_entrada) # Simulando a saída de um sistema linear

"""Reta tangente ao degrau negativo"""
amplitude_degrau_negativo = y_degrau[0][49999]                                                           # Definindo a amplitude do degrau positivo
coeficiente_angular_degrau_negativo = coeficiente_angular(y_degrau[0], time, 37600)                 # Definindo o coeficiente angular pela função
equacao_degrau_negativo = equacao(time, y_degrau[0], coeficiente_angular_degrau_negativo, 37600)    # Definindo a eqaução da reta
ref3 = amplitude_degrau_negativo * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao ponto de operação
"""Determinando os parâmetros de Ziegler através da análise gráfica para o degrau negativo, o valor de theta é de 0.012, o que pode ser considerado aproximadamente nulo:"""
tau_negativo =  38.8910 - 37.5130
var_theta_n, var_u_n, K2 = valor_K(y_degrau[0], U_Degraus_conc, ponto_radiano, U_entrada, 49999)
# print(f'K-negativo = {K2}, Variaçao Saída = {var_theta_n}, Variação U = {var_u_n}, tau = {tau_negativo}')

#----------------------------------------------------------------------------------------------------------------------------------------
"Modelo de primeira ordem com os parâmetros de Ziegler-Nichols para o degrau negativo" 
#----------------------------------------------------------------------------------------------------------------------------------------

modelo_degrau_negativo = K2/(tau_negativo*s+1)                                                          # Modelo do degrau positivo
# print(modelo_degrau_negativo)
t_degrau_negativo, y_degrau_negativo = ct.forced_response(modelo_degrau_negativo,T= time, U = U_Degraus_conc - U_entrada) # Simulando a saída de um sistema linear

#----------------------------------------------------------------------------------------------------------------------------------------
"Gráficos método de Ziegler-Nichols" 
#----------------------------------------------------------------------------------------------------------------------------------------
plt.figure(3)
plt.title("Gráfico Método de 3 parâmetros - degrau positivo")
plt.plot(time,y_degrau[0],'b',label='$\\theta$')
plt.plot(time,y_degrau_positivo + ponto_radiano,'g',label='$\\theta_{Ziegler-Nichols}$')
plt.plot(time,equacao_degrau_positivo,'r',label='reta tangente')
plt.plot(time,ref,'--k',label='ref')
plt.plot(time,ref2,'--m',label='ref_degrau')
plt.grid()
plt.ylabel('$\\theta[rad]$')
plt.legend()
plt.ylim(0.650,0.75)
plt.xlim(11,25)

x0 = 12.50
xf = 12.52
y0 = 0.663
yf = 0.664

x0_ = 12
xf_ = 13
y0_ = 0.666
yf_ = 0.66

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

a = plt.axes([0.35, 0.5, 0.15, 0.15]) # left, right, width, height
plt.xlim(x0,xf)
plt.ylim(y0,yf)
plt.grid()
plt.plot(time,y_degrau[0],'b',label='$g(s)$')
plt.plot(time,equacao_degrau_positivo,'r',label='ref degrau')
plt.plot(t_degrau_negativo,y_degrau_positivo + ponto_radiano,'g',label='$\\theta_{3 parâmetros}$')
plt.plot(time,ref,'--k',label='ref')
plt.show()

"""Degrau negativo"""
plt.figure(4)
plt.title("Gráfico Método de 3 parâmetros - degrau negativo")
plt.plot(time,y_degrau[0],'b',label='$\\theta$')
plt.plot(t_degrau_negativo,y_degrau_negativo + ponto_radiano,'g',label='$\\theta_{3 parâmetros}$')
plt.plot(time,equacao_degrau_negativo,'r',label='reta tangente')
plt.plot(time,ref,'--k',label='ref')
plt.plot(time,ref3,'--c',label='ref_degrau')
plt.grid()
plt.ylabel('$\\theta[^{\\circ}C]$')
plt.legend()
plt.ylim(0.57,0.68)
plt.xlim(36,50)

x0_n = 37.495
xf_n = 37.52
y0_n = 0.662
yf_n = 0.6635

x0_ = 37
xf_ = 38
y0_ = 0.66
yf_ = 0.67

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

a = plt.axes([0.35, 0.5, 0.15, 0.15]) # left, right, width, height
plt.xlim(x0_n,xf_n)
plt.ylim(y0_n,yf_n)
plt.grid()
plt.plot(time,y_degrau[0],'b',label='$g(s)$')
plt.plot(time,equacao_degrau_negativo,'r',label='ref degrau')
plt.plot(t_degrau_negativo,y_degrau_negativo + ponto_radiano,'g',label='$\\theta_{Ziegler-Nichols}$')
plt.plot(time,ref,'--k',label='ref')
plt.show()

"""Realizando o atraso por Padé"""
#1ª ordem
n1 = 1 
N1 = ct.pade(0.15,n1) #aproximação de padé
Gd1 = ct.TransferFunction(np.array(N1[0]),np.array(N1[1])) #construção da função de transferência.
Gr1 = modelo_degrau_positivo*Gd1
print(Gr1)
t1, y1 = ct.forced_response(Gr1,T= time, U = U0-U_entrada) 
y1 = y1 + ponto_radiano

#3ª ordem
n3 = 3 
N3 = ct.pade(0.15,n3) #aproximação de padé
Gd3 = ct.TransferFunction(np.array(N3[0]),np.array(N3[1])) #construção da função de transferência.
Gr3 = modelo_degrau_positivo*Gd3
t3, y3 = ct.forced_response(Gr1,T= time, U = U0-U_entrada) 
y3 = y3 + ponto_radiano

#5ª ordem
n5 = 5 
N5 = ct.pade(0.15,n5) #aproximação de padé
Gd5 = ct.TransferFunction(np.array(N5[0]),np.array(N5[1])) #construção da função de transferência.
Gr5 = modelo_degrau_positivo*Gd5
print(Gr5)
t5, y5 = ct.forced_response(Gr5,T= time, U = U0-U_entrada) 
y5 = y5 + ponto_radiano

#9ª ordem
n9 = 9 
N9 = ct.pade(0.15,n9) #aproximação de padé
Gd9 = ct.TransferFunction(np.array(N9[0]),np.array(N9[1])) #construção da função de transferência.
Gr9 = modelo_degrau_positivo*Gd9
t9, y9 = ct.forced_response(Gr9,T= time, U = U0-U_entrada) 
y9 = y9 + ponto_radiano

#conversão das funções de transferência do atraso para um sistema em espaço de estado
atraso1 = ct.tf2io(Gd1, name='atraso1', inputs='u', outputs='y1')#1ª ordem
atraso3 = ct.tf2io(Gd3, name='atraso3', inputs='u', outputs='y3')#3ª ordem
atraso5 = ct.tf2io(Gd5, name='atraso5', inputs='u', outputs='y5')#5ª ordem
atraso9 = ct.tf2io(Gd9, name='atraso9', inputs='u', outputs='y9')#9ª ordem

#-------------------------------------------------------------------------
#construção da malha aberta

#1ª ordem------------------------------------------------------
malha_aberta1 = ct.InterconnectedSystem(
    (atraso1, FanPlate), name='malha_aberta1',
    connections = (('atraso1.u',),('FanPlate.u','atraso1.y1')),
    inplist = ('atraso1.u'),
    outlist = ('FanPlate.y'))

X1 = np.zeros(n1+2) #condições iniciais nulas
Ta, y1 = ct.input_output_response(malha_aberta1, time, U_Degraus_conc, X1)
#--------------------------------------------------------------

#3ª ordem------------------------------------------------------
malha_aberta3 = ct.InterconnectedSystem(
    (atraso3, FanPlate), name='malha_aberta3',
    connections = (('atraso3.u',),('FanPlate.u','atraso3.y3')),
    inplist = ('atraso3.u'),
    outlist = ('FanPlate.y'))

X3 = np.zeros(n3+2) #condições iniciais nulas
Ta, y3 = ct.input_output_response(malha_aberta3, time, U_Degraus_conc, X3)
#--------------------------------------------------------------

#5ª ordem------------------------------------------------------
malha_aberta5 = ct.InterconnectedSystem(
    (atraso5, FanPlate), name='malha_aberta5',
    connections = (('atraso5.u',),('FanPlate.u','atraso5.y5')),
    inplist = ('atraso5.u'),
    outlist = ('FanPlate.y'))

X5 = np.zeros(n5+2) #condições iniciais nulas
Ta, y5 = ct.input_output_response(malha_aberta5, time, U_Degraus_conc, X5)
#--------------------------------------------------------------

#9ª ordem------------------------------------------------------
malha_aberta9 = ct.InterconnectedSystem(
    (atraso9, FanPlate), name='malha_aberta9',
    connections = (('atraso9.u',),('FanPlate.u','atraso9.y9')),
    inplist = ('atraso9.u'),
    outlist = ('FanPlate.y'))

X9 = np.zeros(n9+2) #condições iniciais nulas
Ta, y9 = ct.input_output_response(malha_aberta9, time, U_Degraus_conc, X9)
#--------------------------------------------------------------

'''
--------------------------------------------------------------------
Plotagem da resposta temporal do sistema com atraso
--------------------------------------------------------------------
'''
plt.figure(5)
plt.title("Gráfico aproximação por padé com atraso")
plt.plot(Ta,ref,'--k',label='ref')
plt.plot(Ta, y_degrau_positivo + ponto_radiano, 'black', label = 'Modelo - 3P') #3P
plt.plot(Ta, y1, 'purple', label='Modelo - 1ª ordem') #atraso 1 
plt.plot(Ta, y3, 'g', label='Modelo - 2ª ordem') #atraso 3 
plt.plot(Ta, y5, 'r', label='Modelo - 5ª ordem') #atraso 5 
plt.plot(Ta, y9, 'b', label='Modelo - 9ª ordem') #atraso 9 
plt.legend()
plt.ylabel('$\\theta[rad]$')
plt.xlabel('Tempo')
plt.ylim(0.66,0.74)
plt.xlim(10,25)
plt.grid()

x0_n = 12.45
xf_n = 12.70
y0_n = 0.6615
yf_n = 0.6645

x0_ = 12
xf_ = 13
y0_ = 0.66
yf_ = 0.67

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

a = plt.axes([0.4, 0.5, 0.15, 0.15]) # left, right, width, height
plt.xlim(x0_n,xf_n)
plt.ylim(y0_n,yf_n)
plt.grid()
plt.plot(Ta,ref,'--k',label='ref')
plt.plot(Ta, y_degrau_positivo + ponto_radiano, 'black', label = 'Modelo - 3P') #3P
plt.plot(Ta, y1, 'purple', label='Modelo - 1ª ordem') #atraso 1 
plt.plot(Ta, y3, 'g', label='Modelo - 2ª ordem') #atraso 3 
plt.plot(Ta, y5, 'r', label='Modelo - 5ª ordem') #atraso 5 
plt.plot(Ta, y9, 'b', label='Modelo - 9ª ordem') #atraso 9 
plt.show()

"""Comparando o degrau positivo de Ziegler com a 5 aproximação de Padé"""
plt.figure(6)
plt.title("Comparando o degrau positivo de Ziegler com a 5 aproximação de Padé")
plt.plot(Ta,ref,'--k',label='ref')
plt.plot(Ta, y_degrau_positivo + ponto_radiano, 'black', label = 'Modelo - 3P') #3P
plt.plot(Ta, y5, 'r', label='Modelo - 5ª ordem') #atraso 5 
plt.legend()
plt.ylabel('$\\theta[rad]$')
plt.xlabel('Tempo')
plt.ylim(0.66,0.74)
plt.xlim(10,25)
plt.grid()
x0_n = 12.45
xf_n = 12.70
y0_n = 0.6615
yf_n = 0.6645

x0_ = 12
xf_ = 13
y0_ = 0.66
yf_ = 0.67

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

a = plt.axes([0.4, 0.5, 0.15, 0.15]) # left, right, width, height
plt.xlim(x0_n,xf_n)
plt.ylim(y0_n,yf_n)
plt.grid()
plt.plot(Ta,ref,'--k',label='ref')
plt.plot(Ta, y_degrau_positivo + ponto_radiano, 'black', label = 'Modelo - 3P') #3P
plt.plot(Ta, y5, 'r', label='Modelo - 5ª ordem') #atraso 5 
plt.show()

"""
Projetando os controladores P e PI
"""
"""Método da curva de Ziegler - Nichols"""
theta_atraso = 0.15

num, dem = ct.pade(theta_atraso, 5) # Definindo as aproximações do atraso por pade de 5ª ordem
Gd7 = ct.tf(num, dem) # Função de transferência da aproximação do atraso de 5ª ordem
Atraso = ct.tf2io(Gd7,name ='atraso',inputs='u', outputs='y')

#Determinanddo os valores do controlador Proporcional
Kc_P_ZN = tau_positivo / (K1 * theta_atraso)
num_P_ZN = [Kc_P_ZN] 
den_P_ZN = 1
G_P_ZN = ct.tf(num_P_ZN,den_P_ZN)

#Determinanddo os valores do controlador Proporcional Integral
Kc_PI_ZN = (0.9*tau_positivo)/(K1*theta_atraso) #Kc tabela 8.7
T_PI_ZN=(theta_atraso*10)/3 #Ti tabela 8.6
num_PI_ZN = [Kc_PI_ZN * T_PI_ZN, Kc_PI_ZN]
den_PI_ZN = [T_PI_ZN, 0]
G_PI_ZN = ct.tf(num_PI_ZN,den_PI_ZN)

"""Método da curva de CHR"""
#Determinanddo os valores do controlador Proporcional
Kc_P_CHR = (0.3*tau_positivo) / (K1 * theta_atraso)
num_P_CHR = [Kc_P_CHR] 
den_P_CHR = 1
G_P_CHR = ct.tf(num_P_CHR,den_P_CHR)

#Determinanddo os valores do controlador Proporcional Integral
Kc_PI_CHR = (0.6*tau_positivo)/(K1*theta_atraso) #Kc tabela 8.7
T_PI_CHR= 4*theta_atraso #Ti tabela 8.6
num_PI_CHR = [Kc_PI_CHR * T_PI_CHR, Kc_PI_CHR]
den_PI_CHR = [T_PI_CHR, 0]
G_PI_CHR = ct.tf(num_PI_CHR,den_PI_CHR)

time, y_degrau_G_P_ZN = ct.forced_response(G_P_ZN,T= time, U = U_Degraus_conc - U_entrada) # Simulando a saída de um sistema linear
print(f"proporcional Ziegler Nichls: {G_P_ZN}")
print(f"proporcional-integral Ziegler Nichls: {G_PI_ZN}")
print(f"proporcional CHR: {G_P_CHR}")
print(f"proporcional-integral CHR: {G_PI_CHR}")


"""Atividade 10-a"""

tempo_final = 240 # Duração da simulação (s)
T = 0.001         # Período de amostragem (s)
time = np.arange(0, tempo_final, T) # Vetor do tempo de simulação espaçados de 0.01s
u_massa = np.zeros(len(time))
U0 = U_entrada * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao valor de entrada
ref_4 = np.zeros(len(time))

controlador_P_ZN = ct.tf2io(G_P_ZN, name='controlador_P_ZN', inputs='u', outputs='y') #controlador P da Curva de Reação Ziegler-Nichols
controlador_PI_ZN = ct.tf2io(G_PI_ZN, name='controlador_PI_ZN', inputs='u', outputs='y') #controlador PI da Curva de Reação Ziegler-Nichols
controlador_P_CHR = ct.tf2io(G_P_CHR, name='controlador_P_CHR', inputs='u', outputs='y') #controlador P do CHR
controlador_PI_CHR = ct.tf2io(G_PI_CHR, name='controlador_PI_CHR', inputs='u', outputs='y') #controlador PI do CHR

# print(controlador_P_ZN)
# print(controlador_PI_ZN)
# print(controlador_P_CHR)
# print(controlador_PI_CHR)

# saturacao = ct.NonlinearIOSystem(0, saturation_output, states=0, name='saturacao', inputs=('u'), outputs=('y5'))

#----------------------------------------------------------------------------------------------------------------------------------------
" Definindo os degraus do sistema para malha fechada"
#----------------------------------------------------------------------------------------------------------------------------------------

for k in range(len(time)):
    if time[k] < 30: # Estabilizando no ponto de equilíbrio
        ref_4[k] = ponto_radiano
    elif time[k] < 60: # Degrau unitário superior
        ref_4[k] = ponto_radiano + np.radians(3.64)
    elif time[k] < 90: # Condição de equilíbrio
        ref_4[k] = ponto_radiano
    elif time[k] < 120: # Degrau unitário inferior
        ref_4[k] = ponto_radiano - np.radians(3.64)
    elif time[k] < 150: # Metade do degrau unitário superior
        ref_4[k] = ponto_radiano + np.radians(1.82)
    elif time[k] < 180: # Degrau unitário superior
        ref_4[k] = ponto_radiano + np.radians(3.64)
    elif time[k] < 210: # Metade do degrau unitário inferior
        ref_4[k] = ponto_radiano - np.radians(1.82)
    else: # Condição de equilíbrio
        ref_4[k] = ponto_radiano

#------------------------------------------------------------------------------
#Interconexão entre o controlador Proporcional (Zigler-Nichols), Sistema não linear e Atraso
#------------------------------------------------------------------------------

FanPlate_P_ZN = ct.InterconnectedSystem(
    (controlador_P_ZN, Atraso, FanPlate), name='FanPlate_P_ZN',
    connections = (('controlador_P_ZN.u', '-FanPlate.y'),
                   ('atraso.u','controlador_P_ZN.y'),
                   ('FanPlate.u','atraso.y')),
    inplist = ('controlador_P_ZN.u','atraso.u', 'FanPlate.u_massa'),
    inputs = ('theta_ref_P_ZN', 'U0', 'u_massa'),
    outlist = ('FanPlate.y', 'atraso.u'),
    outputs = ('y','u'))

FanPlate_PI_ZN = ct.InterconnectedSystem(
    (controlador_PI_ZN, Atraso, FanPlate), name='FanPlate_PI_ZN',
    connections = (('controlador_PI_ZN.u', '-FanPlate.y'),
                   ('atraso.u','controlador_PI_ZN.y'),
                   ('FanPlate.u','atraso.y')),
    inplist = ('controlador_PI_ZN.u','atraso.u', 'FanPlate.u_massa'),
    inputs = ('theta_ref_PI_ZN', 'U0', 'u_massa'),
    outlist = ('FanPlate.y', 'atraso.u'),
    outputs = ('y','u'))

#controladorador P - método CHR
FanPlate_P_CHR = ct.InterconnectedSystem(
    (controlador_P_CHR, Atraso, FanPlate), name='FanPlate_P_CHR',
    connections = (('controlador_P_CHR.u', '-FanPlate.y'),
                   ('atraso.u','controlador_P_CHR.y'),
                   ('FanPlate.u','atraso.y')),
    inplist = ('controlador_P_CHR.u','atraso.u', 'FanPlate.u_massa'),
    inputs = ('theta_ref_P_CHR', 'U0', 'u_massa'),
    outlist = ('FanPlate.y', 'atraso.u'),
    outputs = ('y','u'))

#controladorador PI - método CHR
FanPlate_PI_CHR = ct.InterconnectedSystem(
    (controlador_PI_CHR, Atraso, FanPlate), name='FanPlate_PI_CHR',
    connections = (('controlador_PI_CHR.u', '-FanPlate.y'),
                   ('atraso.u','controlador_PI_CHR.y'),
                   ('FanPlate.u','atraso.y')),
    inplist = ('controlador_PI_CHR.u','atraso.u', 'FanPlate.u_massa'),
    inputs = ('theta_ref_PI_CHR', 'U0', 'u_massa'),
    outlist = ('FanPlate.y', 'atraso.u'),
    outputs = ('y','u'))

posicao_inicial = np.radians(35) #Condição inicial do sistema (posição).
velocidade_inicial = 0 #Condição inicial do sistema (velocidade).
u_massa = np.zeros(len(time))

t_P_PI, y_P_ZN = ct.input_output_response(FanPlate_P_ZN, time, [ref_4, U0, u_massa], [0,0,0,0,0,0.9*posicao_inicial,velocidade_inicial]) #controlador P do método da Curva de Reação Ziegler-Nichols
t_P_PI, y_PI_ZN = ct.input_output_response(FanPlate_PI_ZN, time, [ref_4, U0, u_massa], [0,0,0,0,0,0,0.9*posicao_inicial,velocidade_inicial]) #controlador PI do método da Curva de Reação Ziegler-Nichols
t_P_PI, y_P_CHR = ct.input_output_response(FanPlate_P_CHR, time, [ref_4, U0, u_massa], [0,0,0,0,0,0.9*posicao_inicial,velocidade_inicial]) #controlador P do método CHR
t_P_PI, y_PI_CHR = ct.input_output_response(FanPlate_PI_CHR, time, [ref_4, U0, u_massa], [0,0,0,0,0,0,0.9*posicao_inicial,velocidade_inicial]) #controlador PI do método CHR

print(y_P_ZN[0])


plt.figure(7)
plt.subplot(2,1,1)
plt.title("Controladores P e PI pelo método da curva de Ziegler-Nichols e CHR")
plt.plot(time, ref_4,'--k', label='ref')
plt.plot(time, y_P_ZN[0],'b', label='P-ZN') #controlador P do método da Curva de Reação Ziegler-Nichols
plt.plot(time, y_PI_ZN[0],'r', label='PI-ZN') #controlador PI do método da Curva de Reação Ziegler-Nichols
plt.plot(time, y_P_CHR[0],'g', label='P-CHR') #controlador P do método CHR
plt.plot(time, y_PI_CHR[0],'purple', label='PI-CHR') #controlador PI do método CHR
plt.ylabel('$\\theta[rad]$')
plt.xlim(0,max(time))
plt.ylim(0.525, 0.8)
plt.legend(loc='upper right')
plt.grid()

x0_n = 28.000
xf_n = 34.000
y0_n = 0.65
yf_n = 0.78

# ==============================
# PLOT NORMAL
# ==============================

plt.plot([x0_n,xf_n], [y0_n,y0_n], 'c--')
plt.plot([x0_n,xf_n], [yf_n,yf_n], 'c--')
plt.plot([x0_n,x0_n], [y0_n,yf_n], 'c--')
plt.plot([xf_n,xf_n], [y0_n,yf_n], 'c--')
# ==============================
# PLOT COM ZOOM
# ==============================

a = plt.axes([0.2, 0.57, 0.10, 0.10]) # left, right, width, height
plt.xlim(x0_n,xf_n)
plt.ylim(y0_n,yf_n)
plt.grid()
plt.plot(time, ref_4,'--k', label='ref')
plt.plot(time, y_P_ZN[0],'b', label='theta(t)') #controlador P do método da Curva de Reação Ziegler-Nichols
plt.plot(time, y_PI_ZN[0],'r', label='theta(t)') #controlador PI do método da Curva de Reação Ziegler-Nichols
plt.plot(time, y_P_CHR[0],'g', label='theta(t)') #controlador P do método CHR
plt.plot(time, y_PI_CHR[0],'purple', label='theta(t)') #controlador PI do método CHR

plt.subplot(2,1,2)
plt.plot(time, y_P_ZN[1],'b', label='P-ZN') #controlador P do método da Curva de Reação Ziegler-Nichols
plt.plot(time, y_PI_ZN[1],'r', label='PI-ZN') #controlador PI do método da Curva de Reação Ziegler-Nichols
plt.plot(time, y_P_CHR[1],'g', label='P-CHR') #controlador P do método CHR
plt.plot(time, y_PI_CHR[1],'purple', label='PI-CHR') #controlador PI do método CHR
plt.ylabel('U')
plt.xlabel('Tempo')
plt.xlim(0,max(time))
plt.legend(loc='best')
plt.grid()
plt.show()

for k in range(len(time)): # Aplica a variação na massa do sistema entre 60s e 160s
    if time[k] >= 60 and time[k] < 160:
        u_massa[k] = 0.2

ref_10 = ponto_radiano * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao ponto de operação

t1_, y1_ = ct.input_output_response(FanPlate_P_ZN, time, [ref_10,U0,u_massa], [0,0,0,0,0,0.9*posicao_inicial,velocidade_inicial]) #controlador P do método da Curva de Reação Ziegler-Nichols
t2_, y2_ = ct.input_output_response(FanPlate_PI_ZN, time, [ref_10,U0,u_massa], [0,0,0,0,0,0,posicao_inicial,velocidade_inicial]) #controlador PI do método da Curva de Reação Ziegler-Nichols
t3_, y3_ = ct.input_output_response(FanPlate_P_CHR, time, [ref_10,U0,u_massa], [0,0,0,0,0,0.9*posicao_inicial,velocidade_inicial]) #controlador P do método da Curva de CHR
t4_, y4_ = ct.input_output_response(FanPlate_PI_CHR, time, [ref_10,U0,u_massa], [0,0,0,0,0,0,posicao_inicial,velocidade_inicial]) #controlador PI do método da Curva de CHR

plt.figure(8)
plt.subplot(2,1,1)
plt.title("Controladores P e PI pelo método da curva de Ziegler-Nichols e CHR - Rejeição a perturbação")
plt.plot(t1_,ref_10,'--k',label='ref')
plt.plot(t1_,y1_[0,:],'mediumvioletred',label='P-ZN')
plt.plot(t2_,y2_[0,:],'goldenrod',label='PI-ZN')
plt.plot(t3_,y3_[0,:],"-.b",label='P-CHR')
plt.plot(t4_,y4_[0,:],"-.g",label='PI-CHR')
plt.ylabel('$\\theta(t) [rad]$')
plt.xlim(0,tempo_final)
plt.legend()
plt.legend(prop={'size':10})
plt.grid()

plt.subplot(2,1,2)
plt.plot(t1_,y1_[1,:],'mediumvioletred',label='P-ZN')
plt.plot(t2_,y2_[1,:],'goldenrod',label='PI-ZN')
plt.plot(t3_,y3_[1,:],"-.b",label='P-CHR')
plt.plot(t4_,y4_[1,:],"-.g",label='PI-CHR')
plt.ylabel('U')
plt.xlabel('Tempo')
plt.xlim(0,tempo_final)
plt.legend()
plt.legend(prop={'size':10})
plt.grid()
plt.show()

"""Atividade 11"""
"""
    Normalizando o sistema
    ======================
"""
# Normalizando a curva de resposta com o controlador 1 (P_zn)
yn_P_zn = (y_P_ZN[0] - np.radians(38))/np.radians(3.64) # Normalização

# Normalizando a curva de resposta com o controlador 2 (PI_zn)
yn_PI_zn = (y_PI_ZN[0] - np.radians(38))/np.radians(3.64) # Normalização

# Normalizando a curva de resposta com o controlador 3 (P_chr)
yn_P_chr = (y_P_CHR[0] - np.radians(38))/np.radians(3.64) # Normalização

# Normalizando a curva de resposta com o controlador 4 (PI_chr)
yn_PI_chr = (y_PI_CHR[0] - np.radians(38))/np.radians(3.64) # Normalização

# Vetores de referência para os gráficos normalizados
refn = np.ones(time.shape)              # Referência do ponto de equilíbrio
reftp = (1+ 0.02)*np.ones(time.shape)   # Referência da margem de +2%
reftn = (1- 0.02)*np.ones(time.shape)   # Referência da margem de -2%

plt.figure(9)
plt.subplot(2,2,1)
plt.title("Respostas dos controladores P")
plt.plot(time,yn_P_zn,label='P-ZN')
plt.plot(time, refn*yn_P_zn[49999], "--r", label = 'Ref')
plt.plot(time, reftp*yn_P_zn[49999], "--k", label = '$\pm 2\%$ ref')
plt.plot(time, reftn*yn_P_zn[49999], "--k")
plt.ylabel('$\\theta[rad]$')
plt.xlim(29.4,36)
plt.ylim(-0.05,1.4)
plt.legend()
plt.grid()

# Controlador 2 (PI_zn)----------------------------------------------------------
plt.subplot(2,2,2)
plt.title("Respostas dos controladores PI")
plt.plot(time,yn_PI_zn,label='PI-ZN')
plt.plot(time, refn*yn_PI_zn[49999], "--r", label = 'Ref')
plt.plot(time, reftp*yn_PI_zn[49999], "--k", label = '$\pm 2\%$ ref')
plt.plot(time, reftn*yn_PI_zn[49999], "--k")
plt.xlim(29.4,36)
plt.ylim(-0.05,1.75)
plt.legend()
plt.grid()

# Controlador 3 (P_chr)-----------------------------------------------------------
plt.subplot(2,2,3)
plt.plot(time,yn_P_chr,label='P-CHR')
plt.plot(time, refn*yn_P_chr[49999], "--r", label = 'Ref')
plt.plot(time, reftp*yn_P_chr[49999], "--k", label = '$\pm 2\%$ ref')
plt.plot(time, reftn*yn_P_chr[49999], "--k")
plt.ylabel('$\\theta[rad]$')
plt.xlim(29.4,36)
plt.ylim(-0.05,0.8)
plt.legend()
plt.grid()

# Controlador 4 (PI_chr)--------------------------------------------------------
plt.subplot(2,2,4)
plt.plot(time,yn_PI_chr,label='PI-CHR')
plt.plot(time, refn*yn_PI_chr[49999], "--r", label = 'Ref')
plt.plot(time, reftp*yn_PI_chr[49999], "--k", label = '$\pm 2\%$ ref')
plt.plot(time, reftn*yn_PI_chr[49999], "--k")
plt.xlabel('Tempo')
plt.xlim(29.4,36)
plt.ylim(-0.05,1.4)
plt.legend()
plt.grid()
plt.show()

# =============================================================================
# Parâmetros obtidos
# =============================================================================

# P_ZN (Controlador 1)
a1 = 0.883
b1 = 1.374 - a1
tp1 = 30.45      # tempo de pico
tr1 = 30.30      # tempo de subida
ts1 = 32.27     # tempo de acomodação
mp1 = b1/a1     # sobressinal máximo

# PI_ZN (Controlador 2)
a2 = 1
b2 = 1.719 - a2
tp2 = 30.51      # tempo de pico
tr2 = 30.31      # tempo de subida
ts2 = 32.74     # tempo de acomodação
mp2 = b2/a2     # sobressinal máximo

# P_CHR (Controlador 3)
"""Se trata de um sistema criticamente amortecido, desse modo, não é possível obter o parâmetros referentes ao guia 4"""

# PI_CHR (Controlador 4)
a4 = 1
b4 = 1.338 - a4
tp4 = 30.0196     # tempo de pico
tr4 = 30.397      # tempo de subida
ts4 = 31.089     # tempo de acomodação
mp4 = b4/a4     # sobressinal máximo



"""
    Avaliando a performance por meio de índices de desempenho
    =========================================================
"""
#------------------------------------------------------------------------------
#==============================================================================
#                                 Questão 12                                  =
#==============================================================================
# Para todos os métdos abaixo, os 4 primeiros valores dos vetores estão relacionados à sequência de degraus.
# Já os 4 últimos estão relacionados à rejeição de perturbação.
# Calculando o erro em relação à referência ----------------------------------------------------------------------------
e0 = np.abs(y_P_ZN[0] - ref_4) # Calculo do erro do Controlador Proporcional por Zigler-Nichols para a sequencia de degraus
e1 = np.abs(y_PI_ZN[0] - ref_4) # Calculo do erro do Controlador Proporcional Integral por Zigler-Nichols para a sequencia de degraus
e2 = np.abs(y_P_CHR[0] - ref_4) # Calculo do erro do Controlador Proporcional por CHR para a sequencia de degraus
e3 = np.abs(y_PI_CHR[0] - ref_4) # Calculo do erro do Controlador Proporcional Integral por CHR para a sequencia de degraus

e4 = np.abs(y1_[0] - ref_10) # Calculo do erro do Controlador Proporcional por Zigler-Nichols para a rejeição à perturbação
e5 = np.abs(y2_[0] - ref_10) # Calculo do erro do Controlador Proporcional Integral por Zigler-Nichols para a rejeição à perturbação
e6 = np.abs(y3_[0] - ref_10) # Calculo do erro do Controlador Proporcional por CHR para a rejeição à perturbação
e7 = np.abs(y4_[0] - ref_10) # Calculo do erro do Controlador Proporcional Integral por CHR para a rejeição à perturbação

# Cálculo pelo método IAE ----------------------------------------------------------------
IAE = np.zeros(8) # Alocando um vetor de 8 posições para o erro pelo método de IAE

IAE[0] = T*np.sum(e0) # Calculo do IAE Controlador Proporcional por Zigler-Nichols para a sequencia de degraus
IAE[1] = T*np.sum(e1) # Calculo do IAE Controlador Proporcional Integral por Zigler-Nichols para a sequencia de degraus
IAE[2] = T*np.sum(e2) # Calculo do IAE Controlador Proporcional por CHR para a sequencia de degraus
IAE[3] = T*np.sum(e3) # Calculo do IAE Controlador Proporcional Integral por CHR para a sequencia de degraus

IAE[4] = T*np.sum(e4) # Calculo do IAE Controlador Proporcional por Zigler-Nichols para a rejeição à perturbação
IAE[5] = T*np.sum(e5) # Calculo do IAE Controlador Proporcional Integral por Zigler-Nichols para a rejeição à perturbação
IAE[6] = T*np.sum(e6) # Calculo do IAE Controlador Proporcional por CHR para a rejeição à perturbação
IAE[7] = T*np.sum(e7) # Calculo do IAE Controlador Proporcional Integral por CHR para a rejeição à perturbação


print(f"IAE ZN P é:  {IAE[0]}")
print(f"IAE ZN PI é: {IAE[1]}")
print(f"IAE CHR P é: {IAE[2]}")
print(f"IAE CHR PI é: {IAE[3]}")

print(f"Rejeição a perturbação - IAE ZN P  é: {IAE[4]}")
print(f"Rejeição a perturbação - IAE ZN PI  é: {IAE[5]}")
print(f"Rejeição a perturbação - IAE CHR P  é: {IAE[6]}")
print(f"Rejeição a perturbação - IAE CHR PI  é: {IAE[7]}")


# Cálculo pelo método ITAE para degrais seccionados ----------------------------------------------------------------

ITAE = np.zeros(8) # Alocando um vetor de 8 posições para o erro pelo método de ITAE

ITAE[0] = T*np.sum(time[0:30000]*(e0[0:30000] + e0[30000:60000] + e0[60000:90000] + e0[90000:120000] + e0[120000:150000] + e0[150000:180000] + e0[180000:210000] + e0[210000:240000])/8) # Calculo do ITAE Controlador Proporcional por Zigler-Nichols para a sequencia de degraus
ITAE[1] = T*np.sum(time[0:30000]*(e1[0:30000] + e1[30000:60000] + e1[60000:90000] + e1[90000:120000] + e1[120000:150000] + e1[150000:180000] + e1[180000:210000] + e1[210000:240000])/8) # Calculo do ITAE Controlador Proporcional Integral por Zigler-Nichols para a sequencia de degraus
ITAE[2] = T*np.sum(time[0:30000]*(e2[0:30000] + e2[30000:60000] + e2[60000:90000] + e2[90000:120000] + e2[120000:150000] + e2[150000:180000] + e2[180000:210000] + e2[210000:240000])/8) # Calculo do ITAE Controlador Proporcional por CHR para a sequencia de degraus
ITAE[3] = T*np.sum(time[0:30000]*(e3[0:30000] + e3[30000:60000] + e3[60000:90000] + e3[90000:120000] + e3[120000:150000] + e3[150000:180000] + e3[180000:210000] + e3[210000:240000])/8) # Calculo do ITAE Controlador Proporcional Integral por CHR para a sequencia de degraus

ITAE[4] = T*np.sum(time[0:30000]*(e4[0:30000] + e4[30000:60000] + e4[60000:90000] + e4[90000:120000] + e4[120000:150000] + e4[150000:180000] + e4[180000:210000] + e4[210000:240000])/8) # Calculo do ITAE Controlador Proporcional por Zigler-Nichols para a rejeição à perturbação;
ITAE[5] = T*np.sum(time[0:30000]*(e5[0:30000] + e5[30000:60000] + e5[60000:90000] + e5[90000:120000] + e5[120000:150000] + e5[150000:180000] + e5[180000:210000] + e6[210000:240000])/8) # Calculo do ITAE Controlador Proporcional Integral por Zigler-Nichols para a rejeição à perturbação;
ITAE[6] = T*np.sum(time[0:30000]*(e6[0:30000] + e6[30000:60000] + e6[60000:90000] + e6[90000:120000] + e6[120000:150000] + e6[150000:180000] + e6[180000:210000] + e7[210000:240000])/8) # Calculo do ITAE Controlador Proporcional por CHR para a rejeição à perturbação;
ITAE[7] = T*np.sum(time[0:30000]*(e7[0:30000] + e7[30000:60000] + e7[60000:90000] + e7[90000:120000] + e7[120000:150000] + e7[150000:180000] + e7[180000:210000] + e7[210000:240000])/8) # Calculo do ITAE Controlador Proporcional Integral por CHR para a rejeição à perturbação;

print(f"ITAE ZN P é: {ITAE[0]}")
print(f"ITAE ZN PI é: {ITAE[1]}")
print(f"ITAE CHR P é: {ITAE[2]}")
print(f"ITAE CHR PI é: {ITAE[3]}")

print(f"Rejeição a perturbação - ITAE ZN P  é: {ITAE[4]}")
print(f"Rejeição a perturbação - ITAE ZN PI  é: {ITAE[5]}")
print(f"Rejeição a perturbação - ITAE CHR P  é: {ITAE[6]}")
print(f"Rejeição a perturbação - ITAE CHR PI  é: {ITAE[7]}")

# Cálculo pelo método RMSE ----------------------------------------------------------------

RMSE = np.zeros(8) # Alocando um vetor de 8 posições para o erro pelo método de RMSE

RMSE[0] = np.sqrt((e0**2).mean()) # Calculo do RMSE Controlador Proporcional por Zigler-Nichols para a sequencia de degraus
RMSE[1] = np.sqrt((e1**2).mean()) # Calculo do RMSE Controlador Proporcional Integral por Zigler-Nichols para a sequencia de degraus
RMSE[2] = np.sqrt((e2**2).mean()) # Calculo do RMSE Controlador Proporcional por CHR para a sequencia de degraus
RMSE[3] = np.sqrt((e3**2).mean()) # Calculo do RMSE Controlador Proporcional Integral por CHR para a sequencia de degraus

RMSE[4] = np.sqrt((e4**2).mean()) # Calculo do RMSE Controlador Proporcional por Zigler-Nichols para a rejeição à perturbação
RMSE[5] = np.sqrt((e5**2).mean()) # Calculo do RMSE Controlador Proporcional Integral por Zigler-Nichols para a rejeição à perturbação
RMSE[6] = np.sqrt((e6**2).mean()) # Calculo do RMSE Controlador Proporcional por CHR para a rejeição à perturbação
RMSE[7] = np.sqrt((e7**2).mean()) # Calculo do RMSE Controlador Proporcional Integral por CHR para a rejeição à perturbação


print(f"RMSE ZN P é: {RMSE[0]}")
print(f"RMSE ZN PI é: {RMSE[1]}")
print(f"RMSE CHR P é: {RMSE[2]}")
print(f"RMSE CHR PI é: {RMSE[3]}")

print(f"Rejeição a perturbação - RMSE ZN P  é: {RMSE[4]}")
print(f"Rejeição a perturbação - RMSE ZN PI  é: {RMSE[5]}")
print(f"Rejeição a perturbação - RMSE CHR P  é: {RMSE[6]}")
print(f"Rejeição a perturbação - RMSE CHR PI  é: {RMSE[7]}")
 
