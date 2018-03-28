# Dynamical Systems model for Paper on Herbivory, Detritivory, Litter Decomposition, and Plant Growth

library(FME)

library(readr)
inoutfull <- read_csv("~/Dropbox (Personal)/Yale PhD/Dissertation/Chapter 2/animals_litterdecomp_plants/experimentaldataformodel.csv")
head(inoutfull)

# Get the plotting functions need for displaying the results
source("plotfiguresfunction.R")


# Load the dynamical systems model and run using inital parameter  --------

fullmodel_R <- function(pars, io=inoutfull){
  
  fullmodel <-function(t, y,pars){
    
    with(as.list(c(pars,y)),{
      
      mfr = (Vlm*L*M/(Klm + M))/((Vsm*S*M/(Ksm + M)) + Vlm*L*M/(Klm + M))
      
      dL = -Vlm*L*M/(Klm + M) - Vlw*L*W/(Klw + L)
      dM = SUE*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - Vlw*mfr*M*W/(Klw + L) - tm*M
      dW = SUEW*(Vlw*L*W/(Klw + L) + Vlw*mfr*M*W/(Klw + L)) - tw*W
      dN = IN - q*N - fi*N + fo*S + (1-SUE)*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - Vpf*N*P/(Kpf + N) 
      
      dS = tm*M - Vsm*S*M/(Ksm + M) + fi*N - fo*S + cv1*(tp*P+ tw*W)
      dP = Vpf*N*P/(Kpf + N) - tp*P
      dD = cv2*(tp*P+ tw*W)
      
      dL15 = -Vlm*L15*M/(Klm + M) - Vlw*L15*W/(Klw + L)
      dM15 = SUE*(Vlm*L15*M/(Klm + M) + Vsm*S15*M/(Ksm + M)) - Vlw*mfr*M15*W/(Klw + L) - tm*M15
      dW15 = SUEW*(Vlw*L15*W/(Klw + L) + Vlw*mfr*M15*W/(Klw + L)) - tw*W15
      dN15 = IN*Rin -q*N15 +  (1-SUE)*(Vlm*L15*M/(Klm + M) + Vsm*S15*M/(Ksm + M)) - Vpf*N15*P/(Kpf + N)
      
      dS15 = tm*M15 - Vsm*S15*M/(Ksm + M) + cv1*(tp*P15 + tw*W15)
      dP15 = Vpf*N15*P/(Kpf + N) - tp*P15
      dD15 = cv2*(tp*P15 + tw*W15)
      
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
    
    output[,,rep] = ode(y=y, times = 1:117, func = fullmodel, parms = pars) #, method="euler")
    
  }
  
  output
}

# Initial parameters used to begin model selection
params<- c(Vlm = 0.025,
           Klm = 0.0005,
           Vlw = 0.011,
           Klw = 0.3,
           Vsm = 0.1,
           Ksm = 3000,
           Vpf = 0.1,
           Kpf = 10,
           SUE = 0.50,
           SUEW = 0.8,
           q = 0.01,
           IN=3,
           Rin=0.00367738555411437,
           tm = 0.1,
           tw = 0.01,
           tp = 0.03,
           fi=0.5,
           fo=0.002,
           cv1=1,
           cv2=0)

out = fullmodel_R(pars=params)


# Run model comparison to fielddata ---------------------------------------

fielddata = as.data.frame(cbind(inoutfull$L,inoutfull$M,inoutfull$W,inoutfull$N,inoutfull$S,inoutfull$P,
             inoutfull$L15,inoutfull$M15,inoutfull$W15,inoutfull$N15,inoutfull$S15,inoutfull$P15))

names(fielddata) = c("L","M","W","N","S","P",
                     "L15","M15","W15","N15","S15","P15")

plotcompare(model=out, data=fielddata, number=6)


# Cost function separating empirical data from model predictions
fullmodelcost = function(pars){
  
  
  # Get model output minus time, detritus, and detritus 15N
  modelout = fullmodel_R(pars=pars)[117,c(-1,-8,-15),]
  
  if(min(modelout) < 0) return(1e6)
  
  
  #Objective function, each variable rescaled by mean value
  sum(t(matrix(rep(c(rep(1,6), rep(1, 6)),40),byrow=T,nrow=40,ncol=12))*((t(fielddata) - modelout)/t(matrix(rep(apply(fielddata,2,sd),40),byrow=T,nrow=40,ncol=12)))^2)
  
}

fullmodelcost(pars=params)

#Updated cost function separating empirical data from model predictions for parameter sense
fullmodel_R_sense <- function(pars){
  
  output = fullmodel_R(pars=pars)
  
  t(output[117,c(-1,-8,-15),])
}


# See sensFun help file for explanation of the collinearity check
Sfun = sensFun(func=fullmodel_R_sense, params)

Sfun2 = sensFun(func=fullmodel_R_sense, params, 
                senspar= c(1,3,5,7,11,12,14,15,16))

summary(Sfun)

plot(Sfun, xlab="time", lwd = 2)

pairs(Sfun)

ident <- collin(Sfun2)
plot(ident, log = "y")
abline(h=15, col="red")
abline(h=10, col="orange")
head(ident, 20)
#First try
subset(collin(Sfun2, N=3), collinearity <10)

subset(collin(Sfun, N=11), collinearity <10)

# Run model WITH SOM AND INPUT/LOSS: Diffuse and indexed as 2 -------------------------------------------------

#start with the same steps as above
params2<- c(Vlm = 0.025,
            Klm = 0.005,
           Vlw = 0.015,
           Klw = 30,
           Vsm = 0.1,
           Ksm = 2500,
           Vpf = 0.15,
           Kpf = 10,
           SUE = 0.50,
           SUEW = 0.8,
           q = 0.2,
           IN=3,
           Rin=0.00367738555411437,
           tm = 0.1,
           tw = 0.01,
           tp = 0.03,
           fi=0.6,
           fo=0.002,
           cv1=1,
           cv2=0)

out2 = fullmodel_R(pars=params2)
range(out2)
range(out2[117,2:8,])
range(out2[117,9:15,])
fullmodelcost(pars=params2)

plotcompare(model=out2, data=fielddata, number=6, gt="trend")
plotcompare(model=out2, data=fielddata, number=6, gt="compare")
plotcompareF15(model=out2, data=fielddata, number=6)
fullmodelcost(pars=params2)

lpars2 = log(c(Vlm = 0.025, Vlw = 0.015, Vsm = 0.1, Vpf = 0.15,tm = 0.1,tp = 0.03))

fullmodelcost2 <- function(lpars2){
  
  pars2 = c(exp(lpars2),
            Klm = 0.005,
            Klw = 30,
            Ksm = 2500,
            Kpf = 10,
            SUE = 0.50,
            SUEW = 0.8,
            q = 0.2,
            IN=3,
            Rin=0.00367738555411437,
            tw = 0.01,
            fi=0.6,
            fo=0.002,
            cv1=1,
            cv2=0)
  
  fullmodelcost(pars=pars2)
  
}

fullmodelcost2(lpars2)


# optimize the fit using the optim function
optim2_1 = optim(lpars2, fullmodelcost2)
write.csv(optim2_1$par, file="optim2_1_par.csv")

out2_fit = fullmodel_R(pars=c(exp(optim2_1$par),Klm = 0.005,
                          Klw = 30,
                          Ksm = 2500,
                          Kpf = 10,
                          SUE = 0.50,
                          SUEW = 0.8,
                          q = 0.2,
                          IN=3,
                          Rin=0.00367738555411437,
                          tw = 0.01,
                          fi=0.6,
                          fo=0.002,
                          cv1=1,
                          cv2=0))

plotcompare(model=out2_fit, data=fielddata, number=6, gt="trend")
plotcompare(model=out2_fit, data=fielddata, number=6, gt="compare")
plotcompareF15(model=out2_fit, data=fielddata, number=6)

params2_fit=c(exp(optim2_1$par),Klm = 0.005,
             Klw = 30,
             Ksm = 2500,
             Kpf = 10,
             SUE = 0.50,
             SUEW = 0.8,
             q = 0.2,
             IN=3,
             Rin=0.00367738555411437,
             tw = 0.01,
             fi=0.6,
             fo=0.002,
             cv1=1,
             cv2=0)

# Run the model WITHOUT SOM sorbtion or decomposition: Tight and indexed as 3 -----------------------------------

params3<- c(Vlm = 0.025,
            Klm = 0.005,
            Vlw = 0.015,
            Klw = 30,
            Vsm = 0,
            Ksm = 0,
            Vpf = 0.15,
            Kpf = 10,
            SUE = 0.50,
            SUEW = 0.8,
            q = 0,
            IN= 0,
            Rin=0.00367738555411437,
            tm = 0.02,
            tw = 0.013,
            tp = 0.0005,
            fi=0,
            fo=0,
            cv1=1,
            cv2=0)

out3 = fullmodel_R(pars=params3)
range(out3)
fullmodelcost(pars=params3)

plotcompare(model=out3, data=fielddata, number=6, gt="trend")
plotcompare(model=out3, data=fielddata, number=6, gt="compare")
plotcompareF15(model=out3, data=fielddata, number=6)

Sfun3 = sensFun(func=fullmodel_R_sense, params3,
                senspar= c(1,3,5,7,14,15,16))

pairs(Sfun3)
ident <- collin(Sfun3)
plot(ident, log = "y")
abline(h=15, col="red")
abline(h=10, col="orange")

lpars3 = log(c(Vlm = 0.025, Vlw = 0.015, Vpf = 0.15,tm = 0.02,tp = 0.0005))

fullmodelcost3 <- function(lpars3){
  
  pars3 = c(exp(lpars3),
            Klm = 0.005,
            Klw = 30,
            Vsm = 0,
            Ksm = 0,
            Kpf = 10,
            SUE = 0.50,
            SUEW = 0.8,
            q = 0,
            IN=0,
            Rin=0.00367738555411437,
            tw = 0.01,
            fi=0,
            fo=0,
            cv1=1,
            cv2=0)
  
  fullmodelcost(pars=pars3)
  
}

fullmodelcost3(lpars3)

#Optimize the results as above
optim3_1 = optim(lpars3, fullmodelcost3)
write.csv(optim3_1$par, file="optim3_1_par.csv")
out3_fit = fullmodel_R(pars=c(exp(optim3_1$par),Klm = 0.005,
                          Klw = 30,
                          Vsm = 0,
                          Ksm = 0,
                          Kpf = 10,
                          SUE = 0.50,
                          SUEW = 0.8,
                          q = 0,
                          IN=0,
                          Rin=0.00367738555411437,
                          tw = 0.01,
                          fi=0,
                          fo=0,
                          cv1=1,
                          cv2=0))

plotcompare(model=out3_fit, data=fielddata, number=6, gt="trend")
plotcompare(model=out3_fit, data=fielddata, number=6, gt="compare")
plotcompareF15(model=out3_fit, data=fielddata, number=6)

params3_fit = c(exp(optim3_1$par),Klm = 0.005,
               Klw = 30,
               Vsm = 0,
               Ksm = 0,
               Kpf = 10,
               SUE = 0.50,
               SUEW = 0.8,
               q = 0,
               IN=0,
               Rin=0.00367738555411437,
               tw = 0.01,
               fi=0,
               fo=0,
               cv1=1,
               cv2=0)

# Plot the results of model manual and automatic fits ----------------------

fielddata2 = cbind(fielddata, inoutfull$NS_Oct,inoutfull$Label,inoutfull$Herbivory, inoutfull$L0, 
                   inoutfull$H0, inoutfull$P0, inoutfull$M0)
names(fielddata2)[13:19] = c("NS_Oct", "Label", "Herbivory", "L0", "H0", "P0", "M0")
fielddata2["Label1"] = ifelse(fielddata2$Label==1, "Fert./Lablel","Unfert./No Label")
fielddata2["Herbivory1"] = ifelse(fielddata2$Herbivory==1, "Herbivory","No Herbivory")

plot_figureLP2(data=fielddata2, model1=out3, model2 = out2)
plot_figureLP2_scaled(data=fielddata2, model1=out3, model2 = out2)

plot_figure(data=fielddata2,modeltight=out3, modelreal = out2)
plot_figureISO2(data=fielddata2, model1=out3, model2 = out2)

# Manual Fit
plot_figureLP2_scaled(data=fielddata2, model1=out3, model2 = out2)
plot_figureLP2_plantunscaled(data=fielddata2, model1=out3, model2 = out2)

# Optimization Fit
plot_figureLP2_scaled(data=fielddata2, model1=out3_fit, model2 = out2_fit, fname = paste0("model_optimfit_scaled_",as.character(Sys.Date()),"LP.pdf"))
plot_figureLP2_plantunscaled(data=fielddata2, model1=out3_fit, model2 = out2_fit, fname = paste0("model_optimfit_plantunscaled",as.character(Sys.Date()),"LP.pdf"))


# Plot trends in soil nitrogen between tight and realistic cycling --------

pdf("figureS2_INtrace_gradientfit.pdf", width=8, height=5)
    par(mfrow=c(1,2), oma=c(1,0,0,0))
    
    plot(out3[,5,1]~out3[,1,1], type="l",
         lwd=2, col=1, main="Tight Cycling", xlab="", ylab=expression(Inorganic~N~(mg[N])),
         ylim=c(min(min(out3[,5,1:8]),min(fielddata[1:8,4])),
                max(max(out3[,5,1:8]),max(fielddata[1:8,4]))))
    points(117, fielddata[1,4], col=1, pch=19)
    
    for(u in 2:8){
      points(out3[,5,u]~out2[,1,u], type="l",
             lwd=2, col=u)
      points(117, fielddata[u,4], col=u, pch=19)
    }
    text(0,60, label="a")
    
    plot(out2[,5,1]~out2[,1,1], type="l",
         lwd=2, col=1, main="Diffuse Cycling", xlab="", ylab="",
         ylim=c(0, 60))
    points(117, fielddata[1,4], col=1, pch=19)
    
    for(u in 2:8){
      points(out2[,5,u]~out2[,1,u], type="l",
             lwd=2, col=u)
      points(117, fielddata[u,4], col=u, pch=19)
    }
    text(0,60, label="b")
    
    mtext(text="Days", side=1, line=-1, outer=TRUE)
dev.off()

# Run the model with stratified inputs ------------------------------------

stratinput = inoutfull

stratinput$P0 = mean(inoutfull$P0)
stratinput$M0 = mean(inoutfull$M0)
stratinput$N0 = mean(inoutfull$N0)
stratinput$S0 = mean(inoutfull$S0)
stratinput$P150 = mean(inoutfull$P150)
stratinput$M150 = mean(inoutfull$M150)
stratinput$N150 = mean(inoutfull$N150)
stratinput$S150 = mean(inoutfull$S150)

colMeans(inoutfull)

#With SOM and INputs
outstd2 = fullmodel_R(pars=params2_fit, io=stratinput)

#Without SOM or INputs
outstd3 = fullmodel_R(pars=params3_fit, io=stratinput)

plot_figureISO_realcompare(data=fielddata2, model1=out2, model2 = outstd2, fname="figureS1_13Mar2018.pdf")
plot_figureISO_realcompare(data=fielddata2, model1=out3, model2 = outstd3, fname="figureS1_tight_13Mar2018.pdf")

# Plot the results of a gradient between tight cycling and diffuse cycling --------

params3_mod = params3_fit
params3_mod["Ksm"] = 2500

eachlength = 4

INgrad = seq(0,3, length=eachlength)
qgrad = seq(0,0.2, length=eachlength)

figrad = seq(0, 0.6, length=eachlength)
fograd = seq(0, 0.002, length=eachlength)

Vsmgrad = seq(0, 0.01138388, length=eachlength)
selmatrix = expand.grid(INgrad, figrad, Vsmgrad)

pgfinalresults = matrix(NA, nrow= dim(selmatrix)[1], ncol=9)

for(i in 1:dim(selmatrix)[1]){
      
      params3_mod["IN"] = selmatrix[i, 1]
      params3_mod["fi"] = selmatrix[i,2]
      params3_mod["Vsm"] = selmatrix[i,3]
      params3_mod["q"] = qgrad[match(params3_mod["IN"], INgrad)]
      params3_mod["fo"] = fograd[match(params3_mod["fi"], figrad)]
      
      pgfinalresults[i,5]=params3_mod["IN"]
      pgfinalresults[i,6]=params3_mod["q"]
      pgfinalresults[i,7]=params3_mod["Vsm"]
      pgfinalresults[i,8]=params3_mod["fi"]
      pgfinalresults[i,9]=params3_mod["fo"]
      
      paramgradient = matrix(NA, nrow = 40, ncol=10)
      
      model_paramgradient = fullmodel_R(pars=params3_mod)
      
      paramgradient[,1] = (model_paramgradient[1,2,] + model_paramgradient[1,9,]- model_paramgradient[117,2,] - model_paramgradient[117,9,])/117
      paramgradient[,2] = (model_paramgradient[117,7,] + model_paramgradient[117,14,])
      paramgradient[,3] = (model_paramgradient[117,7,] + model_paramgradient[117,14,])/colSums(model_paramgradient[117,,])
      
      pgfinalresults[i,1:2] = range(model_paramgradient)
      pgfinalresults[i,3] = summary(lm(paramgradient[,2]~paramgradient[,1]))$r.squared
      pgfinalresults[i,4] = summary(lm(paramgradient[,3]~paramgradient[,1]))$r.squared
}

pdf("gradient_of_open.pdf")
par(mfrow=c(1,3))
plot(pgfinalresults[,3]~pgfinalresults[,5], xlab="IN", ylab=expression(R^2), pch=19)
points(pgfinalresults[,4]~pgfinalresults[,5], xlab="IN", ylab=expression(R^2), pch=19, col="red")

plot(pgfinalresults[,3]~jitter(pgfinalresults[,7]), xlab="Vsm", ylab=expression(R^2), pch=19)
points(pgfinalresults[,4]~jitter(pgfinalresults[,7]), xlab="Vsm", ylab=expression(R^2), pch=19, col="red")

plot(pgfinalresults[,3]~pgfinalresults[,8], xlab="fi", ylab=expression(R^2), pch=19)
points(pgfinalresults[,4]~pgfinalresults[,8], xlab="fi", ylab=expression(R^2), pch=19, col="red")
dev.off()

fisymbol = pgfinalresults[,8]

fisymbol[which(fisymbol==0)] = 1
fisymbol[which(fisymbol<0.3 & fisymbol>0.1)] = 3
fisymbol[which(fisymbol<0.5 & fisymbol>0.3)] = 10
fisymbol[which(fisymbol==0.6)] = 19

fisymbol = round(fisymbol, 0)

pdf("gradient_of_open_2_gradientfit.pdf", height=6, width=11)
par(mfrow=c(1,2), oma=c(1,0,0,0), mar=c(5,4,3,1))
plot(pgfinalresults[,3]~jitter(pgfinalresults[,7]), xlab="", ylab=expression(R^2), pch=fisymbol, 
     col=as.factor(pgfinalresults[,5]), ylim=c(0,1), main="Plant")
abline(h=0.03, lty=2)
legend("topright", legend=c("0,0", "1, 7e-2", "2, 1.3", "3, 0.2"), col=as.factor(pgfinalresults[,5]), pch=19, bty="n", title = expression(paste(I[N],",",~q)))
text(0,1, label="a")

plot(pgfinalresults[,4]~jitter(pgfinalresults[,7]), xlab="", ylab="", pch=fisymbol, 
     col=as.factor(pgfinalresults[,5]), ylim=c(0,1), main="Plant Proportion")

legend("topright", legend=c("0,0", "0.2, 7e-4", "0.4, 1.3e-3", "0.6, 2e-3"), pch=c(1,3,10,19), bty="n", title = expression(paste(f[i],",",~f[o])), y.intersp = 1.1)

abline(h=0.01, lty=2)
text(0,1, label="b")
mtext(text=expression(V[SM]), side=1, line=-1, outer=TRUE)
dev.off()


# Conduct Sensitivity Analysis ----------------------------------------------------

# Tight Model
model_low = fullmodel_R(pars=params3_fit)

model_low2 = log(model_low[117,2:15,]+1, base=10)

sens_matrix = matrix(NA, nrow=20, ncol=14)

for(i in 1:length(params3_fit)){
  
  params_sens = params3_fit
  params_sens[i] = params3_fit[i]*1.01
  
  bttmparam = abs(log(params_sens[i]+1, base=10)-log(params3_fit[i]+1, base=10))
  
  model_sens = fullmodel_R(pars=params_sens)
  
  model_sens2 = log(model_sens[117,2:15,]+1, base=10)
  
  sens_matrix[i,]=apply(abs(model_sens2 - model_low2),1, mean)/bttmparam
  
}

# Diffuse Model
model_low = fullmodel_R(pars=params2_fit)

model_low2 = log(model_low[117,2:15,]+1, base=10)

sens_matrix_diffuse = matrix(NA, nrow=20, ncol=14)

for(i in 1:length(params2_fit)){
  
  params_sens = params2_fit
  params_sens[i] = params2_fit[i]*1.01
  
  bttmparam = abs(log(params_sens[i]+1, base=10)-log(params2_fit[i]+1, base=10))
  
  model_sens = fullmodel_R(pars=params_sens)
  
  model_sens2 = log(model_sens[117,2:15,]+1, base=10)
  
  sens_matrix_diffuse[i,]=apply(abs(model_sens2 - model_low2),1, mean)/bttmparam
  
}


sens_matrix[which(sens_matrix==Inf)] = NA
sens_matrix2 = sens_matrix[,c(-7,-14)]

row.names(sens_matrix2) = names(params3_fit)
colnames(sens_matrix2) = c("L", "M", "W", "N", "S", "P", "L15", "M15","W15", "N15", "S15", "P15")

sens_matrix_diffuse[which(sens_matrix_diffuse==Inf)] = NA
sens_matrix_diffuse2 = sens_matrix_diffuse[,c(-7,-14)]
row.names(sens_matrix_diffuse2) = names(params2_fit)
colnames(sens_matrix_diffuse2) = c("L", "M", "W", "N", "S", "P", "L15", "M15","W15", "N15", "S15", "P15")

View(sens_matrix2)
View(sens_matrix_diffuse2)

# Run model over 100 years to test long-term dynamics -----------------------------------------

fullmodel_R_longtime <- function(pars, io=inoutfull, time=5000){
  
  fullmodel <-function(t, y,pars){
    
    with(as.list(c(pars,y)),{
      
      mfr = (Vlm*L*M/(Klm + M))/((Vsm*S*M/(Ksm + M)) + Vlm*L*M/(Klm + M))
      
      addL = cv2*(tp*P+ tw*W) + io$L0[rep]/365
      addL15 = cv2*(tp*P15 + tw*W15) + io$L150[rep]/365
      
      dL = -Vlm*L*M/(Klm + M) - Vlw*L*W/(Klw + L) + addL
      dM = SUE*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - Vlw*mfr*M*W/(Klw + L) - tm*M
      dW = SUEW*(Vlw*L*W/(Klw + L) + Vlw*mfr*M*W/(Klw + L)) - tw*W
      dN = IN - q*N - fi*N + fo*S + (1-SUE)*(Vlm*L*M/(Klm + M) + Vsm*S*M/(Ksm + M)) - Vpf*N*P/(Kpf + N) 
      
      dS = tm*M - Vsm*S*M/(Ksm + M) + fi*N - fo*S + cv1*(tp*P+ tw*W)
      dP = Vpf*N*P/(Kpf + N) - tp*P
      dD = 0
      
      dL15 = -Vlm*L15*M/(Klm + M) - Vlw*L15*W/(Klw + L) + addL15
      dM15 = SUE*(Vlm*L15*M/(Klm + M) + Vsm*S15*M/(Ksm + M)) - Vlw*mfr*M15*W/(Klw + L) - tm*M15
      dW15 = SUEW*(Vlw*L15*W/(Klw + L) + Vlw*mfr*M15*W/(Klw + L)) - tw*W15
      dN15 = IN*Rin -q*N15 +  (1-SUE)*(Vlm*L15*M/(Klm + M) + Vsm*S15*M/(Ksm + M)) - Vpf*N15*P/(Kpf + N)
      
      dS15 = tm*M15 - Vsm*S15*M/(Ksm + M) + cv1*(tp*P15 + tw*W15)
      dP15 = Vpf*N15*P/(Kpf + N) - tp*P15
      dD15 = 0
      
      return(list(c(dL, dM, dW, dN, dS, dP, dD, 
                    dL15, dM15, dW15, dN15, dS15, dP15, dD15)))
      
    }
    )
  }
  
  reps=40
  
  output = array(NA, dim=c(time,15,reps))
  
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
    
    output[,,rep] = ode(y=y, times = 1:time, func = fullmodel, parms = pars)
    
  }
  
  output
}

params4 = params2_fit
out_Linput = fullmodel_R_longtime(pars=params4, time=365*100)

# params4["cv1"] = 1
# params4["cv2"] = 0
# out = fullmodel_R_longtime(pars=params4, time=365*10)
# 
# out2 = out
# #dim(out)

out = out_Linput


# Fit the long model

lpars_long = log(c(Vlm = 0.025, Vlw = 0.015, Vsm = 0.1, Vpf = 0.15,tm = 0.1,tp = 0.03))

fullmodelcost_long <- function(lpars_long){
  
  pars_long = c(exp(lpars_long),
                Klm = 0.005,
                Klw = 30,
                Ksm = 2500,
                Kpf = 10,
                SUE = 0.50,
                SUEW = 0.8,
                q = 0.2,
                IN=3,
                Rin=0.00367738555411437,
                tw = 0.01,
                fi=0.6,
                fo=0.002,
                cv1=1,
                cv2=0)
  
  
  # Get model output minus time, detritus, and detritus 15N
  modelout = fullmodel_R_longtime(pars=pars_long, time=365*10)[365*10,2:7,]
  
  if(min(modelout) < 0) return(1e6)
  
  
  #Objective function, each variable rescaled by mean value
  sum(t(matrix(rep(c(rep(1,6), rep(1, 6)),40),byrow=T,nrow=40,ncol=6))*((t(fielddata[,1:6]) - modelout)/t(matrix(rep(apply(fielddata[,1:6],2,sd),40),byrow=T,nrow=40,ncol=6)))^2)
  
}

fullmodelcost_long(lpars_long)

optim_long = optim(lpars_long, fullmodelcost_long) # took 4-5 hours to optimize

write.csv(optim_long$par, file="optim_long_par.csv")

out_long_fit = fullmodel_R_longtime(pars=c(exp(optim_long$par),Klm = 0.005,
                                           Klw = 30,
                                           Ksm = 2500,
                                           Kpf = 10,
                                           SUE = 0.50,
                                           SUEW = 0.8,
                                           q = 0.2,
                                           IN=3,
                                           Rin=0.00367738555411437,
                                           tw = 0.01,
                                           fi=0.6,
                                           fo=0.002,
                                           cv1=1,
                                           cv2=0), time=365*100)


# par(mfrow=c(3,2))
# for(i in 2:7){
#   plot(out[,i,1]~out[,1,1], type="l",xlab="Days",
#        lwd=2, col=1, ylab=c("Time","L", "M", "W", "N", "S", "P", "D", "L15", "M15")[i],
#        ylim=c(min(out[,i,1:30]),
#               max(out[,i,1:30])))
#   for(u in 2:30){
#     points(out[,i,u]~out[,1,u], type="l",
#            lwd=2, col=u)
#   }
# }
# 
# out[365*10,2:7,]

out = out_long_fit

pdf("model_longrun2_27Mar2018.pdf", width=7, height=8)

cvoutput = array(NA, dim=c(100,6,9))

sample=seq(1,100)*365
treatments = inoutfull[,c("Label", "Herbivory", "Isopod")]

for(i in 1:100){
  for(j in 2:7){
    
    treatments["data"] = out[sample[i],j,]
    
    meanlong = aggregate(data~., data=treatments, mean)
    sdlong = aggregate(data~., data=treatments, sd)
    
    cvoutput[i,(j-1),1:8] = sdlong$data/meanlong$data
  }
  
  cvoutput[i,,9] = apply(out[sample[i],2:7,],1, sd)/apply(out[sample[i],2:7,],1, mean)
}

titvec = c("Litter", "Microbes", "Isopods", "Inorganic Nitrogen", "Soil Nitrogen", "Plants")

par(mfrow=c(3,2), oma=c(1,1,0,0), mar=c(4,4,3,1))
for(i in 1:6){
  plot(cvoutput[,i,1], type="l", ylim = c(0, max(cvoutput[,i,])), xlab="", 
       ylab="", main=titvec[i])
  points(cvoutput[,i,2], type="l", col="red")
  points(cvoutput[,i,3], type="l", lty=2)
  points(cvoutput[,i,4], type="l", lty=2, col="red")
  points(cvoutput[,i,5], type="b")
  points(cvoutput[,i,6], type="b", col="red")
  points(cvoutput[,i,7], type="b", lty=2)
  points(cvoutput[,i,8], type="b", col="red", lty=2)
  points(cvoutput[,i,9], lwd=2, type="b", col="blue", lty=3)
  
  if(i==1){legend(x=6, y=0.15, legend = c("Nothing", "Fertilized", "Herbivory", "Isopods", "Overall"), 
                  lty=c(1,1,2,1,3), col=c("black", "red", "black","black", "blue"), pch=c(NA,NA,NA,1,1), bty="n",
                  lwd=c(1,1,1,1,2), title="Treatment Group", y.intersp = 0.75)}
}

mtext(text="Coefficient of Variation", side=2, outer = T, line =-1)
mtext(text="Year", side=1, outer = T, line =-1)


R2output = matrix(NA, nrow=100, ncol=6)

sample=seq(1,100)*365
treatments = inoutfull[,c("Label", "Herbivory", "Isopod")]

for(i in 1:100){
  for(j in 2:7){
    
    treatments["data"] = out[sample[i],j,]
    
    R2output[i, (j-1)] =summary(lm(data~Label+Herbivory+Isopod, data=treatments))$r.squared
    
  }
}

par(mfrow=c(3,2), oma=c(1,1,0,0), mar=c(4,4,3,1))
for(i in 1:6){
  plot(R2output[,i], xlab="",ylim=c(0.15, 1), 
       ylab="", main=titvec[i], type="b", col="blue", lwd=2, lty=3)
}

mtext(text=expression(R^2), side=2, outer = T, line =-1)
mtext(text="Year", side=1, outer = T, line =-1)

dev.off()

pdf("model_longrun2_R2_27Mar2018.pdf", width=7, height=7)
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(5, 4, 4, 2) + 0.1)
plot(R2output[,2], xlab="Year",ylim=c(0, 1), 
     ylab=expression(R^2), type="b", col="blue", lwd=2, lty=3, pch=0)
points(R2output[,3], xlab="Year",ylim=c(0.15, 1), 
       ylab=expression(R^2), type="b", col="grey", lwd=2, lty=3, pch=1)
points(R2output[,4], xlab="Year",ylim=c(0.15, 1), 
       ylab=expression(R^2), type="b", col="red", lwd=2, lty=3, pch=2)
points(R2output[,5], xlab="Year",ylim=c(0.15, 1), 
       ylab=expression(R^2), type="b", col="brown", lwd=2, lty=3, pch=7)
points(R2output[,6], xlab="Year",ylim=c(0.15, 1), 
       ylab=expression(R^2), type="b", col="green", lwd=2, lty=3, pch=10)
abline(v = 10, lty = 2, col="black")
arrows(10, 1, 20, 1, length=0.1)
text(25, 0.95, "Projection beyond fit")
legend("bottomright", legend=c("Microbes", "Isopods", "Inorganic N", "Soil N", "Plants"), lwd=2, lty=3,
       pch=c(0,1,2,7,10), col=c("blue", "grey", "red", "brown", "green"), bty="n", title="State Variable")
dev.off()



# Check end for equilibrium
max(out[36500,2:15,]-out[(36500-1),2:15,])

treatments = inoutfull[,c("Label", "Herbivory", "Isopod")]

yearofinterest = 100

pdf("model_longrun_data_27Mar2018_yr100.pdf", width=7, height=7)
par(mfrow=c(3,2), mar=c(3,3,3,1))

for(j in 2:5){
  
  treatments["data"] = out[yearofinterest*365,j,]
  
  boxplot(data~Label*Herbivory*Isopod, data=treatments, main=titvec[(j-1)],
          names=NA, border=c("black", "red","black", "red","black", "red","black", "red"),
          col=(c("white", "white", "lightgray","lightgray","white","white","lightgray","lightgray")),
          boxlty=c(1,1,1,1,3,3,3,3))
  
}
for(j in 6:7){
  
  treatments["data"] = out[yearofinterest*365,j,]
  
  boxplot(data~Label*Herbivory*Isopod, data=treatments, main=titvec[(j-1)],
          names=c("Ctrl", "F", "H", "F + H", "I", "F + I", "H+I", "F+H+I"), 
          border=c("black", "red","black", "red","black", "red","black", "red"),
          col=(c("white", "white", "lightgray","lightgray","white","white","lightgray","lightgray")),
          boxlty=c(1,1,1,1,3,3,3,3))
  
}
dev.off()

