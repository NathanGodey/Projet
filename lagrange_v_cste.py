import scipy as sc
import numpy as np
from scipy import linalg as LS
from numpy import linalg as LA
import matplotlib.pyplot as plt
import time

# ------------ CODE DE CALCUL POUR LE CHAMP DE VITESSE UNIFORME ----------------------

dt=1/256
epsilon=0.000001
nb_steps=256
nbParticles=8192
center=[0.25,0.25]
cov=1/50
dx=dy=1/2**8

vitesse=0.5;




def norme(v):
    return np.sqrt(v[0]**2+v[1]**2)
 

def v(theta):
    return np.array([vitesse*np.cos(theta),vitesse*np.sin(theta)])

class particle:
    x=0
    x0=0
    y=0
    y0=0
    vx=0
    vy=0
    def __init__(self):
        self.x=self.y=0
        

Particles=[]
for i in range(nbParticles):
    Particles.append(particle())
    # Les particules sont initialement disposées selon une gaussienne cetrée en (0.25 , 0.25) de variance 0.02
    Particles[i].x0=Particles[i].x=np.random.normal(center[0],cov)
    Particles[i].y=Particles[i].y0=np.random.normal(center[1],cov)

def point_fixe(A,k, theta):
    f_k = A + dt*v(theta)/2
    t = dt*(k+1)
    
    B = f_k + dt*v(theta)/2

    while norme( B-A) > epsilon :
        A,B = B, f_k + dt*v(theta)/2
    return B
    

def resoudre(X0,n,theta):
    X = np.array([ X0 for k in range(n)] )
    for k in range(n-1):
        X[k+1] = point_fixe(X[k],k,theta)
    
    return X
#-------------------------------------------------------------------------------

X,Y=np.meshgrid(np.linspace(0,1,20),np.linspace(0,1,20))

# Valeurs possibles de theta
THETA = np.array([ np.pi*k/32 for k in range(16)]);

#Le calcul se fait ici pour k = 8
k=8;
theta = THETA[k]

#Calcul des trajectoires des particules
for i_p in range(len(Particles)):
    p=Particles[i_p]
    allXp[i_p]=resoudre([p.x0,p.y0],nb_steps,theta)

for i_t in range(nb_steps):
    # Affichage du champs de vitesse
    Vx,Vy = [v(theta)],[v(theta)]

    plt.quiver(np.linspace(0,1,20),np.linspace(0,1,20),Vx,Vy)
    # Affichage de la trajectoire de chaque particule
    for i_p in range(len(Particles)):
        p=Particles[i_p]
        plt.plot([p.x],[p.y],marker='o',color='red')
        
    
        p.x,p.y=allXp[i_p][i_t]
    plt.show()
    
    plt.pause(0.02)
    plt.clf()
    
# Calcul de la matrice des snapshots pour une valeur de theta
def snap(theta):
    Snapshot = np.zeros([nbParticles*2, nb_steps])
    res=[]
    for i_p in range(nbParticles):
        p=Particles[i_p]
        res=resoudre([p.x0,p.y0],nb_steps,theta)
        for t in range(nb_steps):
            Snapshot[i_p,t] = res[t,0]
            Snapshot[i_p + nbParticles,t] = res[t,1]
    return Snapshot
M=snap(theta);
pos0 = np.copy(M[:,0])
for i in range(nb_steps):
    M[:,i] = M[:,i] -M[:,0]

        

    
#On calcule la décomposition SVD de la matrice M
print("Calcul de la décomposition SVD de la matrice M")
Ut,St,Vt = LS.svd(M,False) 

#On sauve la matrice M obtenue et les vecteurs singuliers à gauche de M
np.save('M_mat', M)
np.save('Ut_mat', Ut)


#On trace les valeurs singulières de M en fonction de leur indice
Sigma = np.diag(St)
sig = St[0:20]

#print(sig)

svec = np.zeros(20) 
for i in range(0,20):
    svec[i] = i 

plt.plot(svec, np.log(sig)/np.log(10)) 
plt.xlabel('i')
plt.ylabel('log 10 ieme valeur singuliere') 
plt.show()

#Visualiser le premier mode

svec2 = np.zeros(2*nbParticles)

plt.plot( np.linspace(np.linspace(0,1,20),np.linspace(0,1,20));
plt.plot( Ut[0,:nbParticles], Ut[0,nbParticles:2*nbParticles])

Vx,Vy=-Ut[:nbParticles,0], -Ut[nbParticles:2*nbParticles,0]


plt.quiver(pos0[:nbParticles],pos0[nbParticles:2*nbParticles],Vx,Vy)

#Calcul de la matrice des snapshots pour les 16 valeurs de theta
A=snap(0)
for k in range(1,8):
    print(k)
    B=snap(THETA[k])
    A = np.concatenate((A,B),axis=1);
np.save('Theta_mat',A)


#------------------------- INTERPOLATION-----------------------------------------


def Lower_left(x,y):
    i = int(np.floor((x+1)/dx))
    j = int(np.floor((1+y)/dy))
    return np.array([i,j])

def phi(k,x,y):
    i,j= Lower_left(x,y)[0],Lower_left(x,y)[1]
    x_tilde = x/dx - i
    y_tilde = y/dy - j
    if k==1:
        return (1-x_tilde)*(1-y_tilde)
    if k==2:
        return x_tilde*(1-y_tilde)
    if k==3:
        return x_tilde*y_tilde
    else :
        return y_tilde*(1-x_tilde)