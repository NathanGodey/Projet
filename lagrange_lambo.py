import scipy as sc
import numpy as np
from scipy import linalg as LS
from numpy import linalg as LA
import matplotlib.pyplot as plt
import time

# ------------ CODE DE CALCUL POUR LE CHAMP DE VITESSE DE LAMB OSEEN ----------------------  

dt=0.002
epsilon=0.00001
nb_steps=256
nbParticles=100
center=[0.5,0.5]
cov=0.02
dx=dy=0.01

gamma=10
nu=0.5
rc_0=0.7





def norme(v):
    return np.sqrt(v[0]**2+v[1]**2)
    
def u_theta(x,y):
    
    r=norme([x,y])

    return [-y/r,x/r]


def v_LO(x,y,t):

    r=norme([x,y])

    rc=np.sqrt(4*nu*t+rc_0**2)
    V=gamma/(2*np.pi*r)*(1-np.exp(-(r/rc)**2))
    return np.array([V*u_theta(x,y)[0],V*u_theta(x,y)[1]])

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

def point_fixe(A, k):
    f_k = A + dt*v_LO(A[0],A[1],k*dt)/2
    t = dt*(k+1)
    
    B = f_k + dt*v_LO(A[0],A[1],t)/2

    while norme( B-A) > epsilon :
        A,B = B, f_k + dt*v_LO(A[0],A[1],k*dt)/2
    return B
    

def resoudre(X0,n):
    X = np.array([ X0 for k in range(n)] )
    for k in range(n-1):
        X[k+1] = point_fixe(X[k],k)
    
    return X
#-------------------------------------------------------------------------------

X,Y=np.meshgrid(np.linspace(0,1,20),np.linspace(0,1,20))


allXp=nbParticles*[1]
#Calcul des trajectoires des particules
for i_p in range(len(Particles)):
    p=Particles[i_p]
    allXp[i_p]=resoudre([p.x0,p.y0],nb_steps)
    
for i_t in range(nb_steps):
    # Affichage du champs de vitesse
    Vx,Vy=v_LO(X,Y,i_t*dt)

    plt.quiver(np.linspace(0,1,20),np.linspace(0,1,20),Vx,Vy)
    
    # Affichage de la trajectoire de chaque particule
    for i_p in range(len(Particles)):
        p=Particles[i_p]
        plt.plot([p.x],[p.y],marker='o',color='red')
        
    
        p.x,p.y=allXp[i_p][i_t]
    plt.show()
    
    plt.pause(0.02)
    plt.clf()

# Calcul de la matrice des snapshots
def snap():
    Snapshot = np.zeros([nbParticles*2, nb_steps])
    res=[]
    for i_p in range(nbParticles):
        p=Particles[i_p]
        res=resoudre([p.x0,p.y0],nb_steps)
        for t in range(nb_steps):
            Snapshot[i_p,t] = res[t,0]
            Snapshot[i_p + nbParticles,t] = res[t,1]
    return Snapshot
    
M=snap();
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
#------------------------------------------------------------------------

    

    
    