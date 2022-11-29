from itertools import cycle
import random
import sys
import subprocess
import numpy as np
from keras.models import Sequential
from keras.layers import Dense, Activation
import matplotlib.pylab as plt
import grafica as gr

tmax=10000
dt=0.5
N=int(tmax/dt)


load_saved_pool = False
save_brains = True
brains = []
fitness_t=[] #it stores average fitness for each generation
total_brains = 50
update=True
generation = 0



def save_pool():
    for i in range(total_brains):
        brains[i].save_weights("Modelos/" + str(i) + ".keras")
    print('se guardaron los modelos')

def create_model():
    model = Sequential() #build it one by one

    model.add(Dense(3, input_shape=(3,)))       # we need phi and theta start and end plus time
    model.add(Activation('relu'))
    model.add(Dense(25, input_shape=(3,)))
    model.add(Activation('relu'))
    model.add(Dense(25, input_shape=(25,)))
    model.add(Activation('relu'))
    model.add(Dense(2, input_shape=(25,)))      #one controls the power of the motors and the other the angle
    model.add(Activation('tanh'))               #constraint the output to the range (-1,1)

    model.compile(loss='mse', optimizer='adam')  #it needs to be compiled even if we wont use gradient based optimization

    return model

def predict(theta_start,theta_target, brain_number):
    global brains 

    time=np.linspace(0,tmax,N)
    time=np.reshape(time,(N,1))
    data=np.ones((N,3))

    data[:,0]=data[:,0]*theta_start
    data[:,1]=data[:,1]*theta_target
    data[:,2]=time.reshape(1,-1)

    input=data
    output = brains[brain_number].predict(input, 200,verbose=0)

    return output #array of values predicted by the network

def model_crossover(parent1, parent2):
  global brains

  weight1 = brains[parent1].get_weights()  #it brings weights and biases
                                           
  weight2 = brains[parent2].get_weights()  

  new_weight1 = weight1
  new_weight2 = weight2
  
  gene = random.randint(0,len(new_weight1)-1) #we change a random weight
  #print(gene,'gene------')
  new_weight1[gene] = weight2[gene]
  new_weight2[gene] = weight1[gene]
  q=np.asarray([new_weight1,new_weight2],dtype=object)
  #print(type(q))
  return q 

def model_mutate(weights,var):
    for i in range(len(weights)):
        for j in range(len(weights[i])):
            if( random.uniform(0,1) < 0.2): #learing rate of 20%
                change = np.random.uniform(-var,var,weights[i][j].shape)
                weights[i][j] += change
            
    return weights

def start():# creates first model generation
    for i in range(total_brains):
        model = create_model()
        brains.append(model) 

    if load_saved_pool:
      for i in range(total_brains):
          brains[i].load_weights("Modelos/" + str(i) + ".keras")

def total_prediction(theta1,theta2):
  resultados=np.zeros((N,2))
  for i in range(total_brains):
    h=predict(theta1,theta2,i)
    resultados=np.append(resultados,h,axis=1)
    
  resultados=resultados.T
  resultados=np.delete(resultados, np.s_[0:2], axis=0)
  return resultados



def evolve(best_fit1,best_fit2):
  global update
  global generation
  global best_brain
  global best_brain2

  #print(best_fit1,best_fit2)
  mutations=[]
  for i in range(total_brains//2):
    cross_weights=model_crossover(best_fit1,best_fit2)
    mutation1=model_mutate(cross_weights[0],0.5)
    mutation2=model_mutate(cross_weights[1],0.5)

    mutations.append(mutation1)
    mutations.append(mutation2)

  i=0
  for brain in brains:
    brain.set_weights(mutations[i])
    i+=1
  generation+=1




def find_best_fit():
  global update
  fitness=np.loadtxt("fitness.txt")
  print(f"fitness promedio {np.mean(fitness)} en la generacion {generation}")
  print(f"fitness max es {np.max(fitness)} en la generacion {generation} ")
  fitness_t.append(np.mean(fitness))
  maxfit1=np.max(fitness)
  best_fit1=np.where(fitness==maxfit1)[0]
  fitness[best_fit1]=0
  maxfit2=np.max(fitness)
  best_fit2=np.where(fitness==maxfit2)[0]

  if len(best_fit1)>1: 
    best_fit1=best_fit1[0]
  if len(best_fit2)>1:
    best_fit2=best_fit2[0]

  return int(best_fit1),int(best_fit2)


def test_no_io():

	process = subprocess.Popen("a.exe", shell=False)
	
	## Get exit codes
	out, err = process.communicate()
	errcode = process.returncode
	#print(errcode)

	process.kill() 
	process.terminate()




start() 


total_generations=100


for i in range(total_generations):
  theta1=np.random.rand()*2*np.pi
  theta2=np.random.rand()*2*np.pi
  resultados=total_prediction(theta1,theta2)
  np.savetxt(f"datos_gen{1}.txt",resultados) #it stores all the predictions from all the networks                                                    # 2 rows per brain, one for delta_theta and another for thrust 
  np.savetxt(f"condiciones_gen{1}.txt",np.asarray([[theta1,theta2]]))
  test_no_io()
  gr.graficar(i)
  bf1,bf2=find_best_fit()
  evolve(bf1,bf2)
  print(bf1,bf2)
  if i%10==0:
    save_pool()


