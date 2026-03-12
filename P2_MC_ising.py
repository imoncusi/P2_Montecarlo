
import numpy as np
from numpy.random import rand  
import matplotlib.pyplot as plt
import os


#----------------------------------------------------------------------
##  BLOCK OF FUNCTIONS USED IN THE MAIN CODE
#----------------------------------------------------------------------

#Generation of a random initial state for NxN spins
def initialstate(N):
    ''' generates a random spin configuration for initial condition'''
    state = 2*np.random.randint(2, size=(N,N))-1 
    return state


H_ext=[0.0,0.2,0.5,0.8,1.0,1.5,2.0]  #definim una llista de valors, els que es vulguin

Tpic_C_llista=[] #un cop acabat els bucles per cada h in H_ext, es fan servir les dades d'aquestes 2 llistes
Tpic_X_llista=[]
r=plt.figure(dpi=100)


def configPlot(ax,config,step,N,T_val): #definim aquesta funció per finalment obtenir una figura amb subplots
   X, Y = np.meshgrid(range(N), range(N))
   ax.clear()
   ax.pcolormesh(X, Y, config, vmin=-1.0, vmax=1.0, cmap='RdBu_r')
   ax.set_title(f'T={T_val:.3f}, step={step}')
   ax.axis("off")
  
for h in H_ext:
  def mcmove(config, beta, H_ext): 
    '''Monte Carlo move using Metropolis algorithm '''
    for i in range(N):  #es fan 2 bulces per tenir un total de N+N iteracions
        for j in range(N):
                #select random spin from NxN system
                a = np.random.randint(0, N) 
                b = np.random.randint(0, N)
                s =  config[a, b]  
                #calculate energy cost of this new configuration (the % is for calculation of periodic boundary condition)
                nb = config[(a+1)%N,b] + config[a,(b+1)%N] + config[(a-1)%N,b] + config[a,(b-1)%N] 
                cost = 2*s*nb + 2*s*h #AQUI ES ON ES POSA h
                if cost < 0:
                    s *= -1  
                elif rand() < np.exp(-cost*beta): 
                    s *= -1                      
                config[a, b] = s               
    return config
  
#This function calculates the energy of a given configuration for the plots of Energy as a function of T
  def calcEnergy(config):
    '''Energy of a given configuration'''
    energy = 0
    for i in range(len(config)): #i=filas
        for j in range(len(config)):  
            S = config[i,j] 
            nb = config[(i+1)%N, j] + config[i,(j+1)%N] + config[(i-1)%N, j] + config[i,(j-1)%N]
            energy += -nb*S -h*S
    return energy/4  

#This function calculates the magnetization of a given configuration
  def calcMag(config):
    '''Magnetization of a given configuration'''
    mag = np.sum(config)
    return mag

#
# MAIN PROGRAM
#
# Initial parameters for calculation
## change the parameter below if you want to simulate a smaller system
  nt      = 2**6        # number of temperature points
  N       = 2**4        # size of the lattice, N x N
  eqSteps = 2**10       # number of MC sweeps for equilibration
  mcSteps = 2**10       # number of MC sweeps for calculation
  
  n1, n2  = 1.0/(mcSteps*N*N), 1.0/(mcSteps*mcSteps*N*N) #constants de normalitzacio
#Generate a random distribution of temperatures to make an exploration
  tm = 2.269;    T=np.random.normal(tm, .64, nt)  #tm es la tc del model teoric, surt d'una formula
  T  = T[(T>1.2) & (T<6.0)];    nt = np.size(T)
#la funcio np.random.normla(centre de la distribucio,desviacio std,num de valors a generar)
  Energy       = np.zeros(nt);   Magnetization  = np.zeros(nt)
  SpecificHeat = np.zeros(nt);   Susceptibility = np.zeros(nt)
#es creen arrays q van guardant les dades Enery   Magne   Speci   Susce
#                                         Enery   Magne   Speci   Susce
#----------------------------------------------------------------------
#  SIMULATION LOOP
#----------------------------------------------------------------------
  print('Starting Simulations at ',len(T),' different temperatures.')
  for m in range(len(T)):  #index de la temperatura escollida, T1 T2...
    E1 = M1 = E2 = M2 = 0
    config = initialstate(N) #llama a la funcion q genera red inicial: linea 11
    iT=1.0/T[m] #inversa de la temperatura actual en unitats reduides
    iT2=iT*iT   #per calcular calor especific
    print('Simulation ',m+1,' of',len(T),' at reduced temperature T=',T[m])

    for i in range(eqSteps):         # equilibrate
        mcmove(config, iT, h)           # Monte Carlo moves

    for i in range(mcSteps):
        mcmove(config, iT, h)
        Ene = calcEnergy(config)     # calculate the energy
        Mag = calcMag(config)        # calculate the magnetisation

        E1 = E1 + Ene
        M1 = M1 + Mag
        M2 = M2 + Mag*Mag
        E2 = E2 + Ene*Ene

        Energy[m]         = n1*E1
        Magnetization[m]  = n1*M1
        SpecificHeat[m]   = (n1*E2 - n2*E1*E1)*iT2
        Susceptibility[m] = (n1*M2 - n2*M1*M1)*iT
        
  
# Plot everything
# 
  plot_folder="./simu" + str(h)
  plots_finals="C:\\Users\\User\\Desktop\\SNN\\Pmontecarlo"
  os.makedirs(plot_folder, exist_ok=True) 

  f = plt.figure(figsize=(20, 10)); # plot the calculated values

  sp =  f.add_subplot(2, 2, 1 );
  plt.plot(T, Energy, 'o', color="#A60628");
  plt.xlabel("Temperature (T)", fontsize=20);
  plt.ylabel("Energy ", fontsize=20);

  sp =  f.add_subplot(2, 2, 2 );
  plt.plot(T, abs(Magnetization), 'o', color="#348ABD");
  plt.xlabel("Temperature (T)", fontsize=20);
  plt.ylabel("Magnetization ", fontsize=20);

  sp =  f.add_subplot(2, 2, 3 );
  plt.plot(T, SpecificHeat, 'o', color="#A60628");
  plt.xlabel("Temperature (T)", fontsize=20);
  plt.ylabel("Specific Heat ", fontsize=20);

  sp =  f.add_subplot(2, 2, 4 );
  plt.plot(T, Susceptibility, 'o', color="#348ABD");
  plt.xlabel("Temperature (T)", fontsize=20);
  plt.ylabel("Susceptibility", fontsize=20);
  plt.savefig(os.path.join(plot_folder, "plotgrande"+str(h)+".png"))
  plt.close()

  idxC=np.argmax(SpecificHeat) #no et torna el valor, sino la posicio de l'array
  idxX=np.argmax(Susceptibility)
  
  Tpic_C_llista.append(T[idxC]) #això et dona un número, el qual es guarda en Tpic. T es una llista, [idx] denota index de la llista
  Tpic_X_llista.append(T[idxX]) #la idea és fer una llista per despres poder fer un grafic
                         # i veure el desplaçament del pic
                         #Per això abans del loop es fa un Tpic_llista
  print("picos hechos")

#Plotejar el desplaçament dels pics
plt.plot(H_ext, Tpic_C_llista,'o', 'g' );
plt.xlabel("H_ext", fontsize=20);
plt.ylabel("Temperature (T) of " \
"max. Specific Heat", fontsize=20);
plt.savefig(os.path.join(plots_finals, "plotC.png"))
plt.close()

plt.plot(H_ext, Tpic_X_llista,'o', 'g' );
plt.xlabel("H_ext", fontsize=20);
plt.ylabel("Temperature (T) of " \
"max. Susceptibility", fontsize=20);
plt.savefig(os.path.join(plots_finals, "plotX.png"))
plt.close()
#
#
#
#
print(Tpic_C_llista)
iterations=100
plt.ion()


for q in Tpic_C_llista:
    T_snapshots=[0.9*q, 1.0*q, 1.1*q] 
    fig,axes=plt.subplots(1,3,figsize=(12,4)) #crear figura amb 3 subplts, un per a cada valor de la llista Tpicllista
    configs=[] #crear configs independents;;;; de manera que per cada valor es facin 3 snapshots en una mateixa figura
    betas=[]

    for T_val in T_snapshots:
      configs.append(initialstate(N)) #per a cada una de les 3 temperatures per a cada Tpic
      betas.append(1.0/T_val)         #es construeis una configuraio inicial d spins independent
    
    #montecarlo en si

    for step in range(iterations):
       for k in range(3):   #k representa l'index per recorrer aquestes 3 temperatures per a cada Tpic per a cada valor de h
          mcmove(configs[k], N, betas[k])
     
          if step%5 == 0:
               configPlot(axes[k], configs[k], step, N, T_snapshots[k])
  
       plt.pause(0.1)   
    plt.suptitle("From left to right: 0.9Tmax, Tmax, 1.1Tmax")
    plt.savefig(f"snapshots_Tpeak_{q:.3f}.png")   
    plt.draw()
 
print("Simulation ended")