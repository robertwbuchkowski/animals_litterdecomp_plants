# Model with a separate detritus pool:Not used in the MS  -------------------------------------

fullmodel_R_detritus <- function(pars, io=inoutfull){
  
  fullmodel <-function(t, y,pars){
    
    with(as.list(c(pars,y)),{
      
      mfr = (Vlm*L*M/(Klm + M) + Vlm*D*M/(Klm + M))/((Vsm*S*M/(Ksm + M)) + Vlm*L*M/(Klm + M) + Vlm*D*M/(Klm + M))
      
      dL = -Vlm*L*M/(Klm + M) - Vlw*L*W/(Klw + L)
      dM = SUE*(Vlm*L*M/(Klm + M) + Vlm*D*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - Vlw*mfr*M*W/(Klw + L) - tm*M
      dW = SUEW*(Vlw*L*W/(Klw + L) + Vlw*mfr*M*W/(Klw + L) + Vlw*D*W/(Klw + D)) - tw*W
      dN = IN - q*N - fi*N + fo*S + (1-SUE)*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M) + Vlm*D*M/(Klm + M)) - Vpf*N*P/(Kpf + N)
      
      dS = tm*M - Vsm*S*M/(Ksm + M) + fi*N - fo*S
      dP = Vpf*N*P/(Kpf + N) - tp*P
      dD = tp*P+ tw*W - Vlm*D*M/(Klm + M) - Vlw*D*W/(Klw + D)
      
      dL15 = -Vlm*L15*M/(Klm + M) - Vlw*L15*W/(Klw + L)
      dM15 = SUE*(Vlm*L15*M/(Klm + M) + Vlm*D15*M/(Klm + M) + Vsm*S15*M/(Ksm + M)) - Vlw*mfr*M15*W/(Klw + L) - tm*M15
      dW15 = SUEW*(Vlw*L15*W/(Klw + L) + Vlw*mfr*M15*W/(Klw + L) + Vlw*D15*W/(Klw + D)) - tw*W15
      dN15 = IN*Rin -q*N15 +  (1-SUE)*(Vlm*L15*M/(Klm + M) + Vlm*D15*M/(Klm + M) + Vsm*S15*M/(Ksm + M)) - Vpf*N15*P/(Kpf + N)
      
      dS15 = tm*M15 - Vsm*S15*M/(Ksm + M)
      dP15 = Vpf*N15*P/(Kpf + N) - tp*P15
      dD15 = tp*P15 + tw*W15 - Vlm*D15*M/(Klm + M) - Vlw*D15*W/(Klw + D)
      
      return(list(c(dL, dM, dW, dN, dS, dP, dD,
                    dL15, dM15, dW15, dN15, dS15, dP15, dD15)))
      
    }
    )
  }
  
  reps=40
  
  output = array(NA, dim=c(117,15,reps))
  
  for(rep in 1:reps){
    
    y = c(L = io$L0[rep],
          M = io$M0[rep],
          W = io$W0[rep],
          N = io$N0[rep],
          S = io$S0[rep],
          P = io$P0[rep],
          D = io$D0[rep],
          L15 = io$L150[rep],
          M15 = io$M150[rep],
          W15 = io$W150[rep],
          N15 = io$N150[rep],
          S15 = io$S150[rep],
          P15 = io$P150[rep],
          D15 = io$D150[rep])
    
    output[,,rep] = ode(y=y, times = 1:117, func = fullmodel, parms = pars,
                        method="euler")
    
  }
  
  output
}

