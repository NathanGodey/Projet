import scipy as sc
import numpy as np
from scipy import linalg as LS
from numpy import linalg as LA
import matplotlib.pyplot as plt
import time


# ------------ CODE DE CALCUL POUR LE CHAMP DE VITESSE D'UN ECOULEMENT CELLULAIRE ----------------------
    
dt=0.002
epsilon=0.00001
nb_steps=256
nbParticles=100
center=[0.25,0.25]
cov=0.02
dx=dy=0.01

# theta est le vecteur contenant les paramètres de cette simulation
theta = [0.2, 3.12, 2.69] 

def norme(v):
    return np.sqrt(v[0]**2+v[1]**2)

#--------------------------- Champ de vitesse dérivant d'un écoulement cellulaire-------------------------------

def Vx(x,y,theta):
    res = 2*np.pi*np.sin(2*np.pi*x)*np.cos(2*np.pi*y) - theta[0]*2*np.pi*theta[2]*np.cos(2*np.pi*theta[1]*x)*np.sin(2*np.pi*theta[2]*y)
    return res

def Vy(x,y,theta):
    res = 2*np.pi*np.sin(2*np.pi*y)*np.cos(2*np.pi*x) - theta[0]*2*np.pi*theta[1]*np.cos(2*np.pi*theta[1]*y)*np.sin(2*np.pi*theta[2]*x)
    return -res
    
def vitesse_cell(x,y,theta):
    return np.array([ Vx(x,y,theta), Vy(x,y,theta)])
#------------------------------------------------------------------------------------------------------------
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
# Les particules sont initialement disposées selon une gaussienne cetrée en (0.25 , 0.25) de variance 0.02
for i in range(nbParticles):
    Particles.append(particle())
    Particles[i].x0=Particles[i].x=np.random.normal(center[0],cov)
    Particles[i].y=Particles[i].y0=np.random.normal(center[1],cov)

def point_fixe(A,theta):
    f_k = A + dt*vitesse_cell(A[0],A[1],theta)/2
    B = f_k + dt*vitesse_cell(A[0],A[1],theta)/2
    while norme( B-A) > epsilon :
        A,B = B, f_k + dt*vitesse_cell(A[0],A[1],theta)/2
    return B
    

def resoudre(X0,n,theta):
    X = np.array([ X0 for k in range(n)] )
    for k in range(n-1):
        X[k+1] = point_fixe(X[k],theta)
    return X
#-------------------------------------------------------------------------------

X,Y=np.meshgrid(np.linspace(0,1,100),np.linspace(0,1,100))


allXp=nbParticles*[1]


#Calcul des trajectoires des particules
for i_p in range(len(Particles)):
    p=Particles[i_p]
    allXp[i_p]=resoudre([p.x0,p.y0],nb_steps,theta)
    
for i_t in range(nb_steps):
    # Affichage du champs de vitesse
    Vx1,Vy1 = vitesse_cell(X,Y,theta)

    plt.quiver(np.linspace(0,1,100),np.linspace(0,1,100),Vx1,Vy1)
    
    # Affichage de la trajectoire de chaque particule
    for i_p in range(len(Particles)):
        p=Particles[i_p]
        plt.plot([p.x],[p.y],marker='o',color='red')
        
    
        p.x,p.y=allXp[i_p][i_t]
    plt.show()
    
    plt.pause(0.02)
    plt.clf()
    
# Calcul de la matrice des snapshots pour une valeur du vecteur theta
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
#On calcule la décomposition SVD de la matrice M
print("Calcul de la décomposition SVD de la matrice M")
Ut,St,Vt = LS.svd(M,False)
Sigma = np.diag(St)
sig = St[0:40]

print(sig)

svec = np.zeros(40) 
for i in range(0,40):
    svec[i] = i 

# Affichage du log des valeurs singulières
plt.plot(svec, np.log(sig)/np.log(10)) 
plt.xlabel('i')
plt.ylabel('log 10 ieme valeur singuliere') 
plt.show()