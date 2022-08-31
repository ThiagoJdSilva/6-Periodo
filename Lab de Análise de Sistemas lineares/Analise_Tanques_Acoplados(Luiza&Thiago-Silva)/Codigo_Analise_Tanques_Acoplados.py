#Laboratório de Análise de Sistemas Lineares
#@Autores: Thiago José da Silva e Luiza Gomes de Castro e Sá
#@Data: 28/11/2021

import numpy as np                          #Importando a biblioteca numpy
from matplotlib import pyplot as plt        #Importando a biblioteca matplotlib responsável por plotar gráficos
import control as ct                        #Importando a biblioteca control para realizar o controle do sistema
import statistics as sta                    #Importando a biblioteca statistic

plt.close('all') #Fechando todas as abas de gráficos abertas

def model_uptade(t,T,Q,params):
    """
    Função resposável por retornar o valor da EDO
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


def miller(tau):
    """
    Retorna o valor de Tau pelo método de Miller
    """
    return 0.63*tau


SYSTEM = ct.NonlinearIOSystem(model_uptade, model_output,
states=1, name='SYSTEM', inputs=('u'), outputs=('y'))

#----------------------------------------------------------------------------------------------------------------------------------------
" Definindo as varíaveis do sistema"
#----------------------------------------------------------------------------------------------------------------------------------------
V = 10         # Volume do reator
C = 4500       # Capacidade térmica da solução
F = 3          # Vazao volumetrica
h = 15         # Coeficiente de convecção
Aext = 31.4    # Superficie do tanque   
rho = 1000     # Massa especifica da solução
T_inicial = 20 # Temperatura inicial
A_ext = 31.4   # Área do tanque
ponto = 60     # Ponto de operação
T_inicial = 20 # Temperatura ambiente

tempo_final = 99    # Tempo final da análise
time = np.linspace(0, tempo_final, tempo_final+1) # Criando um lista que vai de 0 até 99, com período 1
Q_entrada = h*A_ext*(T_inicial-ponto)- (rho*F*C)*(T_inicial-ponto) # Determina o valor de calor de entrada do sistema
Q0 = Q_entrada * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao valor de entrada
ref = ponto * np.ones(time.shape) # Cria uma lista do tamanho da lista de tempo com todos os valores iguais ao ponto de operação
T0 = 20 # Representa a temperatura ambiente

T,y = ct.input_output_response(SYSTEM, time, Q0, T0, method='Radau')
"""
Simula um sistema dinâmico, EDO de tanques acoplados, com a entrada Q0 e retorna seus valores de tempo de saída e a resposta do sistema.
"""

#----------------------------------------------------------------------------------------------------------------------------------------
" Plotando gráfico da EDO:"
#----------------------------------------------------------------------------------------------------------------------------------------
plt.figure(1)                               # Definindo como figura 1
plt.subplot(2,1,1)                          # Definindo a posição da imagem, primeira coluna e segunda linha
plt.title('Resposta temporal de tanques acoplados em malha aberta') # Nome do Gráfico
plt.plot(T,ref,'--k', label='ref')          # Plotando a primeira curva e sua legenda
plt.plot(T,y,'b', label='$\\theta(t)$')        # Plotando a segunda curva e sua legenda
plt.ylabel('$\\theta(t)[^{\\circ}C]$')      # Definindo o nome do eixo y
plt.xlim(0,(tempo_final/4))                 # Definindo a escala do eixo x
plt.legend()                                # Plotando a legenda na imagem
plt.ylim(0,100)                             # Definindo a escala do eixo y
plt.grid()                                  # Define a presença de grades na imagem

plt.subplot(2,1,2)                          # Definindo a posição da imagem, primeira coluna e segunda linha
plt.plot(time,Q0,'b', label='Q0')           # Plotando a primeira curva 
plt.ylabel('Q(t)')                          # Definindo o nome do eixo y
plt.legend()                                # Plotando as legenda
plt.xlabel('Tempo[s]')                      # Definindo o nome do eixo x
plt.xlim(0,(tempo_final/4))                 # Definindo a escala do eixo x
plt.ylim(0,1000000000)                      # Definindo a escala do eixo y
plt.grid()                                  # Definindo a presença de grades na imagem
plt.show()                                  # Plotando a figura
#----------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------
" Expansão em série de Taylor"
#----------------------------------------------------------------------------------------------------------------------------------------

temperatura_aproximacao = np.empty(len(time)) # Alocando o vetor para a temperatura do liquido
temperatura_aproximacao.fill(np.nan)          # Limpando o vetor e definindo-o como NAN
temperatura_aproximacao[0] = 20               # Definindo a posição 0 do vetor sendo 20

for k in range(len(time)-1):
    """Realiza a aproximação por série de Taylor"""
    temperatura_aproximacao[k+1]= temperatura_aproximacao[k] - 0.3000020933*(temperatura_aproximacao[k] - ponto)

plt.figure(2)
plt.title("Aproximação linear")                                                      
plt.plot(T,ref,'--k', label='ref')
plt.plot(T,temperatura_aproximacao,'b', label='Aprox. Linear')
plt.plot(T,y,'g',label='$\\theta(t)$')
plt.ylabel('$\\theta(t)[^{\\circ}C]$')
plt.xlim(0,(tempo_final/4))
plt.legend()
plt.ylim(0,100)
plt.grid()
plt.show()
#----------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------
" Equação após transformada de Laplace - Modelo Linear"
#----------------------------------------------------------------------------------------------------------------------------------------
s = ct.tf("s") # Cria um sistema de função de transferência.

transformada_de_laplace = 1/(V*rho*C*(s+(F/V)+((h*Aext)/(V*rho*C)))) # Cria a equação da transformada de laplace

t_laplace, y_laplace = ct.forced_response(transformada_de_laplace,T=time, U=Q0)  # Simulando a saída de um sistema linear, são parâmetros a equação linear, o tempo e a varíavel de controle

y_laplace += T_inicial # Somando a temperatura inicial a sáida do sistema

#----------------------------------------------------------------------------------------------------------------------------------------
" Plotando Gráfico da Função de Transferência:"
#----------------------------------------------------------------------------------------------------------------------------------------
plt.figure(3)
plt.subplot(211)
plt.title('Função de transferência') 
plt.plot(t_laplace,y_laplace,'b',label='$\\theta$')
plt.plot(time,ref,'--k',label='ref')
plt.grid()
plt.ylabel('$\\theta(s)[^{\\circ}C]$')
plt.legend()
plt.ylim(0,100)
plt.xlim(0,(tempo_final/4))
plt.subplot(212) 
plt.plot(t_laplace,Q0,'b', label = '$Q0$')
plt.grid()
plt.ylabel('Q (J)')
plt.xlabel('Tempo (s)')
plt.legend()
plt.xlim(0,(tempo_final/4))
plt.ylim(0,1000000000)
plt.show()
#----------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------
"Geração de degraus positivos e negativos em torno do ponto de operação"
#----------------------------------------------------------------------------------------------------------------------------------------
temperatura_degrau = np.empty(len(time))            #Alocando o vetor para a temperatura do liquido
temperatura_degrau.fill(np.nan)                     #Limpando o vetor e definindo-o como NAN

CalorDegraus = Q0[0] * np.ones(len(time))           #Criando um novo vetor para alocar os degraus
CalorDegraus = np.array_split(CalorDegraus, 4)      # Dividindo o vetor em 4 partes iguais

CalorDegraus[0][:]=Q0[0]                            # Definindo a primeira parte do vetor como Q0
CalorDegraus[1][:]=1.10*Q0[0]                       # Definindo a segunda parte do vetor como 110% de Q0
CalorDegraus[2][:]= Q0[0]                           # Definindo a terceira parte do vetor como Q0
CalorDegraus[3][:]=0.90*Q0[0]                       # Definindo a quarta parte do vetor como 90% de Q0
"""Com as linhas acima é definido os degraus"""

CalorDegraus = np.concatenate([CalorDegraus[0], CalorDegraus[1], CalorDegraus[2], CalorDegraus[3]]) # Unindo os vetores que haviam sido divididos

temperatura_degrau[0] = T_inicial                   # Definindo a temperatura ambiente como a posição 0 do novo vetor

for k in range(len(time)-1): 
    temperatura_degrau[k+1] = temperatura_degrau[k]+ ((F*(T_inicial - temperatura_degrau[k])/V)+ (CalorDegraus[k]/(rho*V*C)) + ((h*A_ext)*(T_inicial - temperatura_degrau[k])/(rho*V*C)))

#----------------------------------------------------------------------------------------------------------------------------------------
" Plotando o gráfico dos degraus"
#----------------------------------------------------------------------------------------------------------------------------------------
plt.figure(4)
plt.subplot(211)
plt.title('Degraus')
plt.plot(time, temperatura_degrau,'b',label='$\\theta$')
plt.plot(time,ref,'--k',label='ref')
plt.grid()
plt.ylabel('$\\theta(t)^{\\circ}C$')
plt.legend()
plt.ylim(10,70)
plt.xlim(0,tempo_final)

plt.subplot(212)
plt.plot(time,CalorDegraus,'b',label='$Q_{\\degraus}$')
plt.grid()
plt.ylabel('Q (J) ')
plt.xlabel('Tempo (s)')
plt.legend()
plt.ylim(400000000,700000000)
plt.xlim(0,tempo_final)
plt.show()
#----------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------
" Aplicando o método de Ziegler-Nichols "
#----------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------
#Aplicando para degrau positivo:
#----------------------------------------------------------------------------------------------------------------------------------------

amplitude_degrau_pos = temperatura_degrau[49]                                                           # Definindo a amplitude do degrau positivo
coeficiente_angular_degrau_positivo = coeficiente_angular(temperatura_degrau, time, 26)                 # Definindo o coeficiente angular pela função
equacao_degrau_positivo = equacao(time, temperatura_degrau, coeficiente_angular_degrau_positivo, 26)    # Definindo a eqaução da reta
valores_amplitude_d_pos = amplitude_degrau_pos*np.ones(len(time))                                       # Criando um novo vetor com o tamanho da lista tempo e com os valores iguais o da amplitude
tau_positivo = 28.33 - 25   
                                                                          # Definindo o valor de tau

modelo_degrau_positivo = K1/(tau_positivo*s+1)                                                          # Modelo do degrau positivo

t_degrau_positivo, y_degrau_positivo = ct.forced_response(modelo_degrau_positivo,T= time, U = CalorDegraus) # Simulando a saída de um sistema linear
y_degrau_positivo += T_inicial                                                                              # Somando a temperatura inicial a sáida do sistema

#----------------------------------------------------------------------------------------------------------------------------------------
#Aplicando para degrau negativo
#----------------------------------------------------------------------------------------------------------------------------------------
amplitude_degrau_neg = temperatura_degrau[99]                                                           # Definindo a amplitude do degrau negativo
coeficiente_angular_degrau_negativo = coeficiente_angular(temperatura_degrau, time, 76)                 # Definindo o coeficiente angular pela função
equacao_degrau_negativo = equacao(time, temperatura_degrau, coeficiente_angular_degrau_negativo, 76)    # Definindo a eqaução da reta
valores_amplitude_d_neg = amplitude_degrau_neg*np.ones(len(time))                                       # Criando um novo vetor com o tamanho da lista tempo e com os valores iguais o da amplitude
var_t, var_u2, K2 = valor_K(temperatura_degrau, CalorDegraus, ponto, Q0[0], 49)                                        # Definindo o valor de 
tau_negativo = tau_positivo                                                                             # Definindo o valor de tau negativo

modelo_degrau_negativo = K2/(tau_negativo*s+1)                                                           # Modelo do degrau negativo

t_degrau_negativo, y_degrau_negativo = ct.forced_response(modelo_degrau_positivo,T= time, U = CalorDegraus) # Simulando a saída de um sistema linear
y_degrau_negativo += T_inicial                                                                              # Somando a temperatura inicial a sáida do sistema

#----------------------------------------------------------------------------------------------------------------------------------------
"Gráficos método de Ziegler-Nichols" 
#----------------------------------------------------------------------------------------------------------------------------------------
plt.figure(5)
plt.subplot(211)
plt.title("Gráfico Método de Ziegler-Nichols - degrau positivo")
plt.plot(time,temperatura_degrau,'b',label='$\\theta$')
plt.plot(t_degrau_positivo,y_degrau_positivo,'g',label='$\\theta_{Ziegler-Nichols}$')
plt.plot(time,equacao_degrau_positivo,'r',label='')
plt.plot(time,ref,'--k',label='ref')
plt.plot(time,valores_amplitude_d_pos,'--y',label='K')
plt.grid()
plt.ylabel('$\\theta[^{\\circ}C]$')
plt.legend()
plt.ylim(56,66)
plt.xlim(20,45)

plt.subplot(212)
plt.title("Gráfico Método de Ziegler-Nichols - degrau negativo")
plt.plot(time,temperatura_degrau,'b',label='$\\theta$')
plt.plot(t_degrau_negativo,y_degrau_negativo,'g',label='$\\theta_{Ziegler-Nichols}$')
plt.plot(time,equacao_degrau_negativo,'r',label='')
plt.plot(time,ref,'--k',label='ref')
plt.plot(time,valores_amplitude_d_neg,'--y',label='K')
plt.grid()
plt.ylabel('$\\theta[^{\\circ}C]$')
plt.xlabel('Tempo (s)')
plt.legend(loc='right')
plt.ylim(46,66)
plt.xlim(70,95)
plt.show()

#----------------------------------------------------------------------------------------------------------------------------------------
"Utilizando método de Miller"
#----------------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------------
# Aplicando para Miller Degrau positivo
#----------------------------------------------------------------------------------------------------------------------------------------

tau_miller_positivo = miller(tau_positivo)                                                           # Calculando o valor de tau de Miller pela fórmula na funçao                                                   
modelo_miller_positivo = K1/(tau_miller_positivo*s+1)                                                # Modelo do degrau positivo para Miller
t_miller_pos, y_miller_pos = ct.forced_response(modelo_miller_positivo,T= time, U = CalorDegraus)    # Simulando a saída de um sistema linear
y_miller_pos += T_inicial                                                                            # Somando a temperatura inicial a sáida do sistema

#----------------------------------------------------------------------------------------------------------------------------------------
# Aplicando para Miller Degrau negativo
#----------------------------------------------------------------------------------------------------------------------------------------

tau_miller_negativo = tau_miller_positivo                                                           # Definindo o valor para tau negativo
modelo_miller_negativo = K2/(tau_miller_negativo*s+1)                                               # Modelo do degrau negativo para Miller
t_miller_neg, y_miller_neg = ct.forced_response(modelo_miller_negativo,T= time, U = CalorDegraus)   # Simulando a saída de um sistema linear
y_miller_neg += T_inicial                                                                           # Somando a temperatura inicial a sáida do sistema

#----------------------------------------------------------------------------------------------------------------------------------------
"Gráficos método de Miller" 
#----------------------------------------------------------------------------------------------------------------------------------------
plt.figure(6)
plt.subplot(211)
plt.title("Gráfico Método de Miller - degrau positivo")
plt.plot(time, temperatura_degrau,'b',label='$\\theta_{EDO}$')
plt.plot(t_miller_pos, y_miller_pos,'g',label='$\\theta_{Miller}$')
plt.plot(time,equacao_degrau_positivo,'r',label='')
plt.plot(time,ref,'--k',label='ref')
plt.plot(time,(0.63*valores_amplitude_d_pos),'--y',label='K')
plt.grid()
plt.ylabel('$\\theta[^{\\circ}C]$')
plt.legend()
plt.ylim(56,66)
plt.xlim(20,45)

plt.subplot(212)
plt.title("Gráfico Método de Miller - degrau negativo")
plt.plot(time, temperatura_degrau,'b',label='$\\theta$')
plt.plot(t_miller_neg,y_miller_neg,'g',label='$\\theta_{Miller}$')
plt.plot(time,equacao_degrau_negativo,'r',label='tt')
plt.plot(time,ref,'--k',label='ref')
plt.plot(time,valores_amplitude_d_neg,'--y',label='K')
plt.grid()
plt.ylabel('$\\theta[^{\\circ}C]$')
plt.legend()
plt.ylim(46,66)
plt.xlim(70,95)
plt.show()

#----------------------------------------------------------------------------------------------------------------------------------------
"Plotando todas as curvas"
#----------------------------------------------------------------------------------------------------------------------------------------

plt.figure(7)
plt.title("Curvas de todos os métodos")
plt.plot((time),temperatura_degrau,'b',label='$\\theta_{EDO}$')                                       # Modelo EDO
plt.plot(t_degrau_positivo,y_degrau_positivo,'r',label='$\\theta_{Ziegler-Nichols - positivo}$')    # Modelo degrau positivo método Ziegler-Nichols
plt.plot(t_miller_pos,y_miller_pos,'g',label='$\\theta_{Miller - positivo}$')                       # Modelo degrau positivo método de Miller
plt.plot(t_degrau_negativo,y_degrau_negativo,'y',label='$\\theta_{Ziegler-Nichols - negativo}$')   # Modelo degrau negativo Ziegler-Nicholas
plt.plot(t_miller_neg,y_miller_neg,'c',label='$\\theta_{Miller - positivo}$')                       # Modelo degrau negativo método de Miller
plt.plot(time,ref,'--k',label='ref')
plt.grid()
plt.ylabel('$\\theta[^{\\circ}C]$')
plt.legend()
plt.ylim(10,80)
plt.xlim(0,tempo_final)
plt.show()

#----------------------------------------------------------------------------------------------------------------------------------------
"Comparando os sistemas pelo índice RMSE"
#----------------------------------------------------------------------------------------------------------------------------------------
RMSE = np.zeros(4)

temperatura_nos_degrais = temperatura_degrau[:] # Temperatura no instante dos degrais

erro_degrau_subida_tres_parametros = y_degrau_positivo[:] # Erro da curva parametrizada pelo degrau de subida pelo método dos três parâmetros
RMSE[0] = np.sqrt(sta.mean(((erro_degrau_subida_tres_parametros - temperatura_nos_degrais)**2)))

erro_degrau_descida_tres_parametros = y_degrau_negativo[:] # Erro da curva parametrizada pelo degrau de subida pelo método dos três parâmetros
RMSE[1] = np.sqrt(sta.mean(((erro_degrau_descida_tres_parametros - temperatura_nos_degrais)**2)))

erro_degrau_subida_miller = y_miller_pos[:] # Erro da curva parametrizada pelo degrau de subida pelo método dos três parâmetros
RMSE[2] = np.sqrt(sta.mean(((erro_degrau_subida_miller - temperatura_nos_degrais)**2)))

erro_degrau_descida_miller = y_miller_neg[:] # Erro da curva parametrizada pelo degrau de subida pelo método dos três parâmetros
RMSE[3] = np.sqrt(sta.mean(((erro_degrau_descida_miller - temperatura_nos_degrais)**2)))

print('Índice RMSE para o modelo de Ziegler-Nichols degrau positivo: ', RMSE[0])
print('Índice RMSE para o modelo de Ziegler-Nichols degrau negativo: ', RMSE[1])
print('Índice RMSE para o modelo de Miller degrau positivo:          ', RMSE[2])
print('Índice RMSE para o modelo de Miller degrau negativo:          ', RMSE[3])