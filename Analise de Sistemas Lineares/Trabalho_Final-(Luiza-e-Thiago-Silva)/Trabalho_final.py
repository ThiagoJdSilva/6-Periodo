#@Análise de Sistemas Lineares e Laboratório de Análise de Sistemas Lineares
#@Autores: Thiago José da Silva e Luiza Gomes de Castro e Sá
#@Data: 07/02/2022

import numpy as np                          # Importando a biblioteca numpy
from matplotlib import pyplot as plt        # Importando a biblioteca matplotlib responsável por plotar gráficos
import control as ct                        # Importando a biblioteca control para realizar o controle do sistema
import statistics as sta                    # Importando a biblioteca statistic
import math                                 # Importando a biblioteca math

s = ct.tf('s')          #Criando s como transferência
"""
----------------------------------------------------------------------------------------------------------------------------------------------
Questão 1 à 4
----------------------------------------------------------------------------------------------------------------------------------------------
"""

tempo_final = 100   # Tempo final da análise
periodo = 0.1
time = np.arange(0, tempo_final, periodo) # Criando um lista que vai de 0 até 10, com período 0.1

G_s1 = (-0.125*(s+0.435))/((s**2 + 0.226*s + 0.0169)*(s+1.23)) # Equação de transferência de 3 ordem que representa o sistema
G_s2 = (-0.0442)/(s**2 + 0.226*s + 0.0169) # Equação aproximada de Segunda Ordem para represetar o sistema

t, y = ct.step_response(-G_s1, time)        # Resposta ao degrau negativo para a função de 3 ordem
t, y2 = ct.step_response(-G_s2, time)       # Resposta ao degrau negativo para a função de 2 ordem

informacao1 = ct.step_info(G_s1)            # Devolve os parâmetros de Segunda ordem para a função de 3 ordem
informacao2 = ct.step_info(G_s2)            # Devolve os parâmetros de Segunda ordem para a função de 2 ordem

print(informacao1)
print(informacao2)

ref =  y[999] * np.ones(time.shape)

"""
Plotando a resposta temporal da Função de Transferência de 3 ordem
"""
plt.figure(1)
plt.plot(t, y, 'm', label='Resposta Sistema Terceira Ordem')
plt.plot(t, ref, 'k', label='Ref estabilidade')
plt.plot(t, 0.98*ref, 'b', label='Acomodação $\pm 2%$')
plt.plot(t, 1.02*ref, 'b')
plt.grid()
plt.xlim(0, 100)
plt.ylim(0, 2.8)
plt.ylabel('$\\theta$')
plt.xlabel('Tempo[s]')
plt.legend()
plt.show()

"""
Plotando a resposta temporal da Função de Transferência de 3 ordem e da Função de Transferência de 2 ordem
"""
plt.figure(2)
plt.plot(t, y, 'm', label='Resposta Sistema Terceira Ordem')
plt.plot(t, y2, 'c', label='Resposta Sistema Aproximado Segunda Ordem')
plt.plot(t, ref, 'k', label='Ref estabilidade')
plt.plot(t, 0.98*ref, 'b', label='Acomodação $\pm 2%$')
plt.plot(t, 1.02*ref, 'b')
plt.grid()
plt.xlim(0, 100)
plt.ylim(0, 2.8)
plt.ylabel('$\\theta$')
plt.xlabel('Tempo[s]')
plt.legend()
plt.show()

"""
----------------------------------------------------------------------------------------------------------------------------------------------
Questão 7
----------------------------------------------------------------------------------------------------------------------------------------------
"""
K1 = np.arange(-0.38229, 25.907, 0.1)       # Definindo os valores de K1 que mantém o sistema estável
valores_erro = []                           # Criando uma lista para alocar os valores de erro calculados
A_s = 2/(s+2)                               # Definindo a Equação A(s)
M_s = -s/(0.01*s+1)                         # Definindo a Equação M(s)

Mf1 = (G_s1*A_s*M_s)/(1+ G_s1*A_s*M_s)      # Definindo a malha fechada interna


for k in K1:
    Mf2 = (k*(-Mf1))/(1 + k*(-Mf1))         # Definindo a malha fechada externa - total
    t, y = ct.step_response(Mf2)            # Realizando uma entrada degrau para a malha fechada externa 
    e = 1/(1 +(-k*(-0.1087/0.04157)))       # Calculando o erro estacionário pelo limite s*E(s) com s tendendo 0
    valores_erro.append(e)                  # Alocando o erro na lista

print(min(valores_erro))                    # Printando o menor valor obtido para o erro
print(f"O K1 que apresenta o menor erro é: {K1[valores_erro.index(min(valores_erro))]}")    # Plotando o menor valor de K1 que garante estabilidade a malha 

"""
----------------------------------------------------------------------------------------------------------------------------------------------
Questão 8
----------------------------------------------------------------------------------------------------------------------------------------------
"""
A_s = -2 /(s+2)                             # Definindo a Equação A(s)

FT_ma = G_s1*A_s                            # Definindo a FT em malha aberta

"""
Plotando o Lugar Geométrico das Raízes
"""
plt.figure(3)
rlist, klist = ct.root_locus(FT_ma, print_gain=True, grid=True)
out_st_magn = ct.stability_margins(FT_ma)
Gain_margin = out_st_magn
plt.show()

K1 = 2.828                                 # Valor obtido de K1 para o coeficiente de amortecimento igual a 0.4559

Malha_Fechada_2 =(FT_ma*K1)/(1+(FT_ma*K1))      # Obtendo a equação da malha fechada com K1 = 2.828
t, y1 = ct.step_response(Malha_Fechada_2)       # Realizando uma entrada degrau para a malha fechada 
print(ct.step_info(Malha_Fechada_2))            # Printando os parâmetros de Segunda ordem para a função da malha fechada

ref = y1[999]*np.ones(t.shape)                  # Definindo um vetor referencia 

"""
Plotando a resposta temporal do sistema para K1 = 2.828
"""
plt.figure(4)
plt.plot(t, y1, 'c', label='K = 2.828')
plt.plot(t, ref, 'k', label='Referência')
plt.grid()
plt.legend()
plt.xlim(0, 60)
plt.ylim(0, 1.2)
plt.ylabel('$\\theta$')
plt.xlabel('Tempo[s]')
plt.show()

"""
----------------------------------------------------------------------------------------------------------------------------------------------
Questão 9
----------------------------------------------------------------------------------------------------------------------------------------------
"""

tempo_final = 100   # Tempo final da análise
periodo = 0.1
time = np.arange(0, tempo_final, periodo) # Criando um lista que vai de 0 até 100, com período 0.1

G = -0.125*(s+0.435)/((s+1.23)*(s**2+0.226*s+0.0169))   #função de transferência G
A = 2/(s+2)                                             #função de transferência A

FT_ma = (-0.25*s**2+0.14125*s+0.10875)/(s**4+3.456*s**3+3.20688*s**2+0.610547*s+0.041574) 

plt.figure(5)
rlist, klist = ct.root_locus(FT_ma, print_gain=True, grid=True)
out_st_magn = ct.stability_margins(FT_ma)
Gain_margin = out_st_magn
plt.show()

K1 = 1.049                                              # Ganho encontrado pelo lugar das raízes

C = K1*(-1)                                             # Definindo o controlador para o ganho K1 encontrado
M = K1*s                                                # Definindo M para o ganho K1 encontrado

MF = ct.feedback(G*A,M)                                 # Fechando a malha interna
MF2 = ct.feedback(MF*C,1)                               # Fechando a malha externa

t, y1 = ct.step_response(MF2, time)                     # Realizando uma entrada degrau para a malha fechada 

print(ct.step_info(MF2))                                # Printando os parâmetros de Segunda ordem para a função da malha fechada

ref = y1[999]*np.ones(t.shape)                          # Definindo um vetor referencia

"""
Plotando a resposta temporal do sistema para K1 = 1.049
"""
plt.figure(6)
plt.plot(t, y1, 'c', label='K = 1.049')
plt.plot(t, ref, 'k', label='Referência')
plt.grid()
plt.legend()
plt.xlim(0, 60)
plt.ylim(0, 1.2)
plt.ylabel('$\\theta$')
plt.xlabel('Tempo[s]')
plt.show()

"""
----------------------------------------------------------------------------------------------------------------------------------------------
Questão 10
----------------------------------------------------------------------------------------------------------------------------------------------
"""

valor = list(np.arange(-100, 100,0.1))          # Criando uma lista de valores para plotar a parabola

Ki = 0 * np.ones(len(valor))                    # Valor mínimo de Ki

Kp = []                                         # Array para alocar os valores de Kp

intersecao_x = []                               # Array área de estabilidade em x
intersecao_y = []                               # Array área de estabilidade em y

for kp in valor:
    ki = -0.1621*kp**2 + 4.3077*kp + 4.4698
    if kp > -1 and kp < 0:
        intersecao_x.append(kp + 0.1)           # Valor para Kp estavavel
        intersecao_y.append(ki + 0.1)           # Valor para Ki estavel
    elif kp > 0 and kp < 13:
        intersecao_x.append(kp + 0.1)           # Valor para Kp estavavel
        intersecao_y.append(ki - 0.1)           # Valor para Ki estavel
    elif kp > 13 and kp < 27.57:
        intersecao_x.append(kp - 0.1)           # Valor para Kp estavavel
        intersecao_y.append(ki - 0.1)           # Valor para Ki estavel
    Kp.append(ki)

"""
Plotando a regiãop estavel de Kp x Ki
"""
plt.figure(7)
plt.title("Região estável")
plt.plot(valor,Ki,'--k', label='Minimo de $K_i$')
plt.plot(valor,Kp,'b', label='Valores $K_p$')
plt.plot(valor)
plt.xlim(-50,60)
plt.ylim(-50,50)
for n in list(range(-1,0)):
    plt.vlines(x=intersecao_x[n], ymin=0, ymax=intersecao_y[n], color='m', label='Região estavel')
for n in list(range(0,len(intersecao_y))):
    plt.vlines(x=intersecao_x[n], ymin=0, ymax=intersecao_y[n], color='m')
plt.xlabel("$K_p$")
plt.ylabel("$K_i$")
plt.legend()
plt.grid()
plt.show()

tempo_final = 2                                         # Tempo final da análise
periodo = 0.1                                           # Período de análise
time = np.arange(0, tempo_final, periodo)               # Criando um array que vai de 0 até 10, com período 0.1

G_s = 1 / ((s+1)*(0.2*s+1)*(0.06*s+1))                  # Definindo a função G(s)

valores_y = []                                          # Array para alocar as respostas ao degrau para cada valor de Kp e Ki
par = []                                                # Array para pares de Kp x Ki

ref = 1.2 * np.ones(len(time))                          # Criando uma referência para entrada degrau

"""
A função for abaixo pega para cada valor de Kp, ja demilitado, todos os valores no intervalo de Ki com periodo de 0.1, realizando uma 
análise completa por toda a região de estabilidade.
"""
for n in intersecao_x:
    ki = -0.1621*n**2 + 4.3077*n + 4.4698
    for x in np.arange(0.1,ki-0.1,0.1):
        C_s = n + (x/s)
        malha_fechada = (G_s*C_s)/(1+(G_s*C_s))
        t, y = ct.forced_response(malha_fechada, time, ref)
        valores_y.append(y)
        par.append([n, x])
        print(f'{n} e {x}')

valores_ITAE = []                                       # Vetor para alocar os valores dos erros ITAE
"""
O for abaixo calcula os erros ITAE para todo par de Kp,Ki testado.
"""
for n in valores_y:
    e0 = np.abs(n - ref)                                # Calculo do erro 
    ITAE = periodo * np.sum(time * e0)                  # Calculo erro ITAE
    valores_ITAE.append(ITAE)

print(min(valores_ITAE))                                # Retorna o menor erro obtido
print(f"O par de (Kp, Ki) que apresentaram o melhor desempenho ITAE é: {par[valores_ITAE.index(min(valores_ITAE))]}") # Retorna o par Kp,Ki que produziu o menor erro obtido
