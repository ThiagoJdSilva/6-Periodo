import numpy as np
from matplotlib import pyplot as plt
import control as ct

plt.close('all')

# Variéveis do sistema

V = 10       # Volume do reator químico
c_p = 4500      # Capacidade térmica da solução
F = 3        # Vazao volumétrica
h = 15        # Coeficiente de conveção natural
A_ext = 31.4   # Superficie do tanque   
rho = 1000    # Massa especifica da solução
A_ext = 31.4  # Superfície do tanque para troca de calor por convecção
operationPoint = 45 # Ponto de operação
timeFirst = 20 # Temperatura ambiente

timeEnd = 49 # Tempo de execução do sistema
time = np.linspace(0, timeEnd, timeEnd+1)  # Array do tempo de execução do sistema
heatToPoint = int(-(F*rho*c_p*(timeFirst - operationPoint)) - (h*A_ext*(timeFirst - operationPoint))) # Valor do calor para a temperatura estabilizar no ponto de operação
arrayHeatToPoint = heatToPoint * np.ones(time.shape) # Array de heatToPoint
ref = operationPoint * np.ones(time.shape) # Array de referência para o ponto de operação

# Plot resposta temporal do sistema
def resp_temp():
    
    def model_uptade(t,T,Q,params):
        dT = (Q/(rho*V*c_p))+(h*A_ext*(timeFirst-T)/(rho*V*c_p))+((F*(timeFirst-T)/V))
        return dT

    def model_output(t,T,Q,params):
        return T

    SYSTEM = ct.NonlinearIOSystem(model_uptade, model_output,
    states=1, name='SYSTEM', inputs=('u'), outputs=('y'))

    T,y = ct.input_output_response(SYSTEM, time, arrayHeatToPoint, timeFirst)

    plt.figure(1)
    plt.subplot(2,1,1)
    plt.plot(T,ref,'--k', label='ref')
    plt.plot(T,y,'b', label='T(t)')
    plt.ylabel('T(t)[C]')
    plt.xlim(0,30)
    plt.yticks(range(0, 100, 10))
    plt.legend()
    plt.ylim(0,60)
    plt.title('Resposta temporal de tanques acoplados em malha aberta')
    plt.grid()
    
    plt.subplot(2,1,2)
    plt.plot(time,arrayHeatToPoint,'b', label='Q0')
    plt.ylabel('Q(t)')
    plt.legend()
    plt.xlabel('Tempo[s]')
    plt.xlim(0,30)
    plt.ylim(0,1000000000)
    plt.grid()
    plt.show()

# Plot equação de tranferência
def eq_tranf():
    
    s=ct.tf("s") # Frequência
    eqTransf = 1/(V*rho*c_p*(s+(F/V)+((h*A_ext)/(V*rho*c_p)))) #Equação de transferência
    timeTransfer, tempOut = ct.forced_response(eqTransf,T=time, U=arrayHeatToPoint) #tempOut é a saída da temperatura e timeTransfer é otempo de integraçào

    tempOut += timeFirst # Iniciar a temperatura de saída com a de ambiente
    
    plt.figure(2)
    plt.subplot(211) 
    plt.title('Equação de Transfêrencia')
    plt.plot(timeTransfer,tempOut,'b',label='$\\theta_0$')
    plt.plot(time,ref,'--k',label='ref')
    plt.grid()
    plt.ylabel('$\\theta (°)$')
    plt.legend()
    plt.ylim(15,55)
    plt.xlim(0,30)
    
    plt.subplot(212) 
    plt.plot(timeTransfer,arrayHeatToPoint,'b', label = '$Q1$')
    plt.grid()
    plt.ylabel('Q (J)')
    plt.xlabel('Tempo (s)')
    plt.legend()
    plt.xlim(0,30)
    plt.ylim(0,1000000000)
    plt.show()

# Plot degraus
def degr():

    valueAllocated = np.empty(len(time)) # Separando local na memoria para alocar os valores de temperatura do liquido
    valueAllocated.fill(np.nan) # Preechendo array com Not a Number

    heatSteps = heatToPoint * np.ones(len(time)) # Foi criado um novo array com os degraus 
    heatSteps = np.array_split(heatSteps, 4) # O array foi separado em 4 novos arrays

    # Foi adicionado os degraus para cada intervalo de temperatura
    heatSteps[0].fill(heatToPoint)
    heatSteps[1].fill(1.15 * heatToPoint)
    heatSteps[2].fill(heatToPoint)
    heatSteps[3].fill(0.85 * heatToPoint)     
        
    heatSteps = np.concatenate([heatSteps([0], heatSteps[1], heatSteps[2], heatSteps[3])]) # Faço o concat dos arrays dos degraus
    valueAllocated[0] = timeFirst # Iniciar com um valor previo, pos um Not a Number não pode ser somado na integração abaixo

    # Integrando pelo método de Taylor
    for k in range(len(time)-1): 
        valueAllocated[k+1] = valueAllocated[k]+((F*(timeFirst - valueAllocated[k])/V)+ (heatSteps[k]/(rho*V*c_p)) + ((h*A_ext)*(timeFirst - valueAllocated[k])/(rho*V*c_p)))

    plt.figure(3)
    plt.subplot(211)
    plt.title('Equação Edo')
    plt.plot(time,valueAllocated,'b',label='$\\t$')
    plt.plot(time,ref,'--k',label='ref')
    plt.grid()
    plt.ylabel('$\\t (°)$')
    plt.legend()
    plt.ylim(15,55)
    plt.xlim(0,timeEnd)
    
    plt.subplot(212)
    plt.plot(time,heatSteps,'b',label='Q1')
    plt.grid()
    plt.ylabel('Q (J) ')
    plt.xlabel('Tempo (s)')
    plt.legend()
    plt.ylim(0,1000000000)
    plt.xlim(0,timeEnd)
    plt.show()
    
resp_temp()
eq_tranf()
degr()