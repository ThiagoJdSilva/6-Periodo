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
" Inicio Guia 4"
#----------------------------------------------------------------------------------------------------------------------------------------

def model_update(t,x,u,params):
    # Variáveis de estado
    x1 = x[0] # Posição angular
    x2 = x[1] # Velocidade angular
    return np.array([x2, K_1*(np.cos(x1)**2)*u - (K_2*np.sin(x1) + K_3*x2)]) # Retornando a EDO em forma de vetor


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
m_t = 0.100
C_a = 2.05      # Coeficiente de arrasto
atrito = 5      # Coeficiented e atrito viscoso
g = 9.81        # Gravidade
ponto_de_operacao = 50
ponto_radiano = math.radians(ponto_de_operacao)
ponto_inicial = 0

#----------------------------------------------------------------------------------------------------------------------------------------
" Definindo as constantes K do sistema"
#----------------------------------------------------------------------------------------------------------------------------------------

K_1, K_2, K_3 = K_ganho(m_t)
print(K_1)
print(K_2)
print(K_3)

tempo_final = 20   # Tempo final da análise
periodo = 0.001
time = np.arange(0, tempo_final, periodo) # Criando um lista que vai de 0 até 99, com período 1
U_entrada = (K_2 * math.sin(ponto_radiano))/(K_1 * (math.cos(ponto_radiano)**2))
U0 = U_entrada * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao valor de entrada
ref = ponto_radiano * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao ponto de operação
print(f"valor de entrada {U_entrada}")
print(f"valor de degrau {1.20*U_entrada}")
#----------------------------------------------------------------------------------------------------------------------------------------
" Definindo os degraus do sistema"
#----------------------------------------------------------------------------------------------------------------------------------------
U_Degraus = U_entrada * np.ones(len(time))           #Criando um novo vetor para alocar os degraus
U_Degraus = np.array_split(U_Degraus, 3)      # Dividindo o vetor em 4 partes iguais

U_Degraus[0][:]=U0[0]                            # Definindo a primeira parte do vetor como U0
U_Degraus[1][:]= 1.20*U0[0]                       # Definindo a segunda parte do vetor como 110% de U0
U_Degraus[2][:]= U0[0]                           # Definindo a terceira parte do vetor como U0   
"""Com as linhas acima é definido os degraus"""

U_Degraus_conc = np.concatenate([U_Degraus[0], U_Degraus[1], U_Degraus[2]]) # Unindo os vetores que haviam sido divididos

FanPlate = ct.NonlinearIOSystem(model_update, model_output , states=2, name='FanPlate', inputs = ('u'), outputs = ('y'))

#----------------------------------------------------------------------------------------------------------------------------------------
" Sistema em Malha aberta:"
#----------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------
" Simulação do sistema sem aplicação de degraus na entrada:"
#----------------------------------------------------------------------------------------------------------------------------------------
t, y = ct.input_output_response(FanPlate, time, U0, ponto_inicial)

#----------------------------------------------------------------------------------------------------------------------------------------
" Plotando sistema sem aplicação de degraus na entrada:"
#----------------------------------------------------------------------------------------------------------------------------------------
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,ref,'--k', label='ref')
plt.plot(t,y,'b', label='$\\theta$(t)')
# plt.ylabel('T(t)[C]')
plt.xlim(0,10)
# plt.yticks(range(0, 100, 10))
plt.ylim(0,1.5)
plt.legend()
plt.title('Resposta temporal do sistema em malha aberta sem degrau')
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,U0,'b', label='U (t)')
# plt.ylabel('Q (t)')
# plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,10)
plt.ylim(60,75)
plt.grid()
plt.show()

#----------------------------------------------------------------------------------------------------------------------------------------
" Plotando sistema com aplicação de degraus na entrada:"
#----------------------------------------------------------------------------------------------------------------------------------------
t_degrau, y_degrau = ct.input_output_response(FanPlate, time, U_Degraus_conc, ponto_inicial)
plt.figure(2)
plt.subplot(2,1,1)
plt.plot(t_degrau,ref,'--k', label='ref')
plt.plot(t_degrau,y_degrau,'b', label='$\\theta$(t)')
# plt.ylabel('T(t)[C]')
plt.xlim(0,20)
# plt.yticks(range(0, 100, 10))
plt.ylim(0,1.5)
plt.legend()
plt.title('Resposta temporal do sistema em malha aberta com degrau')
plt.grid()
plt.subplot(2,1,2)
plt.plot(time,U_Degraus_conc,'b', label='U(t)')
# plt.ylabel('Q (t)')
# plt.xlabel('Tempo[s]')
plt.legend()
plt.xlim(0,20)
plt.ylim(60,75)
plt.grid()
plt.show()

#----------------------------------------------------------------------------------------------------------------------------------------
" Normalizando a curva do degrau positivo analisado:"
#----------------------------------------------------------------------------------------------------------------------------------------
yn = y_degrau[6666:13000]
tn = t_degrau[0:len(yn)] 

thetan = (yn - min(yn))/(y_degrau[9999] - min(yn))

#----------------------------------------------------------------------------------------------------------------------------------------
" Referências para acomadação mais ou menos 2% ou mais ou menos 5%:"
#----------------------------------------------------------------------------------------------------------------------------------------
refn = np.ones(tn.shape)          # Referência do ponto de equilíbrio
reftsp = 1.02*np.ones(tn.shape)   # Referência da margem de +2%
reftsn = 0.98*np.ones(tn.shape)   # Referência da margem de -2%
reftsp5 = 1.05*np.ones(tn.shape)  # Referência da margem de +5%
reftsn5 = 0.95*np.ones(tn.shape)  # Referência da margem de -5%

#----------------------------------------------------------------------------------------------------------------------------------------
" Estabelecendo os principais parâmetros que caracterizam a resposta transitória, verificando o último gráfico:"
#----------------------------------------------------------------------------------------------------------------------------------------
tp = 0.348  # Instante de pico
b = 0.5877   # y(tp)
Mp = b/1     # Sobressinal máximo (overshoot)
tr = 0.1909  # Tempo de subida
ts = 2.498  # Tempo de acomodação (2%)
zeta = 0.17465
wn = 9.16848

print(f"valor de zeta:{zeta}")
print(f"valor de wn:{wn}")

#----------------------------------------------------------------------------------------------------------------------------------------
" Determinando as curvas envoltórias:"
#----------------------------------------------------------------------------------------------------------------------------------------
Ep = 1 + (np.e**((-zeta)*wn*tn))/(np.sqrt(1-(zeta**2)))  # Curva envoltória superior
En = 1 - (np.e**((-zeta)*wn*tn))/(np.sqrt(1-(zeta**2)))  # Curva envoltória inferior

#----------------------------------------------------------------------------------------------------------------------------------------
" Plotando Gráfico da resposta e das curvas envoltórias:"
#----------------------------------------------------------------------------------------------------------------------------------------
plt.figure(3)
plt.plot([0.7,0.7], [1,1.5895], ':m', label = 'linhas de cota')
plt.plot([0.3459,0.7], [1.5895,1.5895], ':k')
plt.text(0.75,1.4,"$M_p$")
plt.plot([0,0.3459], [1.5895,1.5895], ':m|')
plt.text(0.115,1.619,"$t_p$")
plt.plot([0.1888,0.1888], [1,1.5], ':k')
plt.plot([0,0.1888], [1.4,1.4], ':m|')
plt.plot([0.1888,0.1888], [1.4,1.4], '--k')
plt.text(0.035,1.43,"$t_r$")
plt.plot([2.5059,2.5059], [1.02,0.77], ':k')
plt.plot([0,2.5059], [0.85,0.85], ':m|')
plt.text(2.1,0.78,"$t_{s,2\%}$")
plt.plot([1.8311,1.8311], [1.0496,0.50], ':k')
plt.plot([0,1.8311], [0.6,0.6], ':m|')
plt.text(1.35,0.535,"$t_{s,5\%}$")
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
plt.legend()
plt.show()

#----------------------------------------------------------------------------------------------------------------------------------------
" Definindo novos degraus para sistema"
#----------------------------------------------------------------------------------------------------------------------------------------
U_Degraus_eq_transf = U_entrada * np.ones(len(time))            #Criando um novo vetor para alocar os degraus
U_Degraus_eq_transf = np.array_split(U_Degraus_eq_transf, 6)    # Dividindo o vetor em 4 partes iguais

U_Degraus_eq_transf[0][:]=U_entrada                             # Definindo a primeira parte do vetor como U0
U_Degraus_eq_transf[1][:]=1.20*U_entrada                        # Definindo a segunda parte do vetor como 110% de U0
U_Degraus_eq_transf[2][:]= U_entrada                            # Definindo a terceira parte do vetor como U0 
U_Degraus_eq_transf[3][:]=0.80*U_entrada                        # Definindo a primeira parte do vetor como U0
U_Degraus_eq_transf[4][:]= U_entrada                            # Definindo a segunda parte do vetor como 110% de U0
U_Degraus_eq_transf[5][:]= 1.30*U_entrada
"""Com as linhas acima é definido os degraus"""

U_Degraus_conc_eq_transf = np.concatenate([U_Degraus_eq_transf[0], U_Degraus_eq_transf[1], U_Degraus_eq_transf[2], U_Degraus_eq_transf[3], U_Degraus_eq_transf[4], U_Degraus_eq_transf[5]]) # Unindo os vetores que haviam sido divididos

t_degrau2, y_degrau2 = ct.input_output_response(FanPlate, time, U_Degraus_conc_eq_transf, 0)

#----------------------------------------------------------------------------------------------------------------------------------------
" Equação após transformada de Laplace - Modelo Linear"
#----------------------------------------------------------------------------------------------------------------------------------------
s=ct.tf("s") # Frequência

#----------------------------------------------------------------------------------------------------------------------------------------
" Determinando o valor do ganho estático"
#----------------------------------------------------------------------------------------------------------------------------------------
K = abs((math.radians(50) - y_degrau[9999])/(1.20*U_entrada- U_entrada))
print(abs(math.radians(50) - y_degrau[9999]))
print(abs(1.20*U_entrada- U_entrada))
print(K)

#----------------------------------------------------------------------------------------------------------------------------------------
" Equação geral para sistemas subamortecidos"
#----------------------------------------------------------------------------------------------------------------------------------------
eqTransf = K*(wn**2)/((s**2)+(2*s*wn*zeta)+(wn**2)) #Equação de transferência
print(eqTransf)
timeTransfer, tempOut = ct.forced_response(eqTransf,T=time, U=U_Degraus_conc_eq_transf-U_entrada) 

timeTransfer, thetanonL = ct.input_output_response(FanPlate, time, U_Degraus_conc_eq_transf, [ponto_radiano,0]) # Sistema Real

#----------------------------------------------------------------------------------------------------------------------------------------
" Plotando Gráfico para comparação do sistema real e do modelo obtido:"
#----------------------------------------------------------------------------------------------------------------------------------------
plt.figure(4)
plt.subplot(211) 
#plt.title('Equação de Transfêrencia')
plt.plot(timeTransfer,tempOut+(np.radians(50)),'b',label='Modelo')
plt.plot(t_degrau2,thetanonL,'y', label='Real')
plt.plot(time,ref,'--k',label='ref(0.873rad)')
plt.grid()
plt.ylabel('$\\theta$(t)[rad]')
plt.legend()
plt.xlim(0,20)
plt.ylim(0.3,1.5)

plt.subplot(212) 
plt.plot(timeTransfer,U_Degraus_conc_eq_transf,'b', label = '$U(t)$')
plt.grid()
plt.ylabel('${u(t)}[m^2/s^2]$')
plt.xlabel('Tempo [s]')
plt.legend()
plt.xlim(0,20)
plt.ylim(45,90)
plt.show()

#----------------------------------------------------------------------------------------------------------------------------------------
" Fim guia 4"
#----------------------------------------------------------------------------------------------------------------------------------------
