import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0,0.5,100)
x = np.linspace(0,1,256)

dx = (1-0)/(256-1)
dt = (0.5-0)/(100-1)

x_0 = 0.25
M = np.zeros((len(x),len(t)))
sigma = 0.05
vitesse = 1
nb_modes = 10

for i in range(len(x)):
    for j in range(len(t)):
        M[i][j] = np.exp(-(x[i]-x_0-vitesse*t[j])**2/(2*sigma**2))

u, s, vh = np.linalg.svd(M,False)

##Trace des colonnes de M en fonction du temps (sur la base u)

#plt.plot((np.diag(s) @ vh)[0,:])
#plt.plot((np.diag(s) @ vh)[1,:])
#plt.plot((np.diag(s) @ vh)[2,:])
#plt.plot((np.diag(s) @ vh)[3,:])

## Portrait de phase

#plt.plot((np.diag(s) @ vh)[0,:],((np.diag(s) @ vh)[1,:]))
#plt.plot((np.diag(s) @ vh)[0,:],((np.diag(s) @ vh)[2,:]))
#plt.plot((np.diag(s) @ vh)[0,:],((np.diag(s) @ vh)[3,:]))

M_app=(u @ np.diag(s)) @ vh



## Modes
#plt.plot(u[:,0])
#plt.plot(u[:,1])
#plt.plot(u[:,2])


##Trace d'une gaussienne quelconque grace aux modes
target=np.zeros(len(x))
for i in range(len(x)):
    target[i]=np.exp(-(x[i]-0.57)**2/(2*sigma**2))
    
#plt.plot(target)

target_app = u[:,0:2] @ np.transpose(u[:,0:2]) @ target

#plt.plot(target_app)



## Plot des courbes approchees et analytiques
#for i in range(5):
    #plt.plot(x,M_app[:,i*10],'b')
    
    #plt.plot(x,M[:,50],'g')
    
## plot des valeurs singulieres normalisees en echelle semi-log      
    
#plt.semilogy([i for i in range(40)],s[:40]/s[0])

plt.show()

##Discretisation de l'equation de transport (terme derive)

Dx_phi=np.zeros((len(x),nb_modes))

for i in range(nb_modes):
    tempder = np.zeros(len(x))
    for j in range(len(x)-1):
        tempder[j]=(u[j+1,i]-u[j,i])/dx
    tempder[len(x)-1]=(u[len(x)-1,i]-u[len(x)-2,i])/dx
    Dx_phi[:,i]=tempder

D=np.dot(np.transpose(u[:,:nb_modes]),Dx_phi)


##Condition initiale

c_0=np.zeros(len(x))
for i in range(len(x)):
    c_0[i]=np.exp(-(x[i]-0.27)**2/(2*sigma**2))
    
a_0=np.transpose(u[:,:nb_modes]) @ c_0
a=a_0
c=c_0
plt.plot(c)
         
##Euler implicite
for j in range (75):
    a = np.linalg.solve(np.eye(nb_modes)+dt*D,a)
    c = u[:,:nb_modes] @ a
    plt.pause(0.1)
    plt.clf()
    plt.plot(c)
plt.show()
    
