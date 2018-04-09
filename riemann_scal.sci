
function [leftf,rightf,lambda] = RIEMANN(lambda1,state1,lambda2,state2,NX,NY,bndryflag)
  
  lambda = max(abs(lambda1),abs(lambda2)); // maximal wave speed
  
  riemann = ( lambda1*state1 - lambda2*state2 + lambda*(state1+state2) )/2.;
  
  rightf = [( lambda1*state1 + lambda2*state2 + lambda*(state1-state2) )/2., 0];
    
  leftf = -rightf; // conservative , entropy flux = 0

endfunction // -- end of Lax-Friedrichs Riemann solver
