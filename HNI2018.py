import math

#-------------------------------------------------------------------------
#------------------------------ Parameters -------------------------------
#-------------------------------------------------------------------------

## vicosity=0

dimsys=1

myfilename="translatinghill"

#----- time range --------------------------------------------------------

Tinitial = 0.

Tfinal = 0.1

Nsteps = 1e4

CFLhard = True

tstpcst = False ## keep same time step as finest grid at each iteration ?

#----- grid --------------------------------------------------------------

mygeofile = "geo_square.py"

Lx = 1.
Ly = 1.

def mymask(xx,yy):
    return False

## levelmax = 7; levelmin = 4;    

levelmax = 5
levelmin = levelmax

errorcomp = True

massshow = True

deltalevel = 1 ## << to compute speed of convergence

myhist = False  ## history of levelmax: to compute order method

#-------------------------------------------------------------------------

onedx = False; onedy=False;
BC=["periodic","periodic","periodic","periodic"]

def fux(t,x,y):
    return 1.
    
def fuy(t,x,y):
    return 1.

## flot analytique : x(t,x0)=x0+t
def piedx(t,x,y):
    return x-t

## flot analytique : x(t,x0)=x0+t    
def piedy(t,x,y):
    return y-t    
    
myx0 = 0.25
myy0 = 0.25


## // x and y must have the same length
def gana(t,x,y):
    xx = piedx(t,x,y)
    yy = piedy(t,x,y)
    mywdt = 0.076
    u = math.exp(-((not onedy)*(xx-myx0)**2+(not onedx)*(yy-myy0)**2)/mywdt**2)*(((xx-myx0)**2+(yy-myy0)**2) < (0.25**2-1e-14) ) 
    return u

mymin=[0]
mymax=[1]
Tinitial = 0.
Tfinal = math.sqrt(2)/4
    
def myentropy1(state):
    return 0
    
def state2plot(i,qin):
    return qin
    
#-------------------------------------------------------------------------
    
localtest = False ## flag to perform or not an entropy test within the Riemann solver

localtest2 = False ## flag to perform or not an entropy test within the Riemann solver

globaltest = False ## flag to perform or not an entropy test after projection step

import riemann_scal

## neuman0 everywhere
def myBC(i,statein,cellcoord):
    xx = cellcoord[0]
    yy = cellcoord[1]

#-------------------------------------------------------------------------

plotevery = 1. ## << postprocessing and plot options

period = Tfinal/1 ## << convert -delay 10  img*.gif anim.gif

myplot = [True,False] ## << to plot or not : as many as dimensions dimsys+1

leg = ["c"] ## <<

sty = []
sty.append([2,-1,3]) ## <<

mygif = False ## << to animate gifs :>  convert -delay 3 img*.gif anim.gif

mypdf = False ## << to flip/flop : convert .gif -flip .jpg

mycontour = False ## << works only if ~onedx ~onedy

mycontourgif = False ## << works only if mycontour

mycontourpdf = False ## << works only if mycontour

mysleep = 10

withini = False

#-------------------------------------------------------------------------
#------------------------------ E0 Parameters -------------------------------
#-------------------------------------------------------------------------

MASS = []
L0L0 = []
L0L1 = []
L0L2 = []
L2L1 = []
L2L2 = []
for i in range(dimsys):
    MASS.append([])
    L0L0.append([])
    L0L1.append([])
    L0L2.append([])
    L2L1.append([])
    L2L2.append([])
    
levelref = levelmax+deltalevel
 
## CONVERGENCE STUDY : START BY FINEST
    
    
    
    
    
    
    