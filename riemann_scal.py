def RIEMANN(lambda1,state1,lambda2,state2,NX,NY,bndryflag):
    lambda0=max(abs(lambda1),abs(lambda2)) ## maximal wave speed
    riemann=(lambda1*state1-lambda2*state2+lambda0*(state1+state2))/2. #??
    rightf=[(lambda1*state1+lambda2*state2+lambda0*(state1-state2))/2.,0]
    leftf=[-rightf[0],0] ## conservative, entropy flux=0
    return [leftf,rightf,lambda0] ## end of Lax-Friedrichs Riemann solver
    