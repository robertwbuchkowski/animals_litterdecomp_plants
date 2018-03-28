# Required plotting functions to run the dynamical systems model
plotcompare <- function(model, data, number =7, reps = c(1,8), gt = "both"){
  
  model = model[,c(-8,-15),]
  
  #combine the SOM with M and N because that was the measured pool
  
  data$S15 = data$S15 - data$M15 - data$N15
  data$S = data$S - data$M - data$N
  
  if(gt=="both"|gt=="trend"){
  
  if(number > 9) par(mfrow=c(4,4), mar = c(4,3,3,1))
  if(number <=9) par(mfrow=c(3,3))
  if(number <=6) par(mfrow=c(3,2))
  if(number ==1) par(mfrow=c(1,1))
  
  for(j in 1:number){
    plot(model[,(j+1),reps[1]]~model[,1,reps[1]], type="l",
         lwd=2, col=1, main=names(data)[j], xlab="Time", ylab="mg/cage",
         ylim=c(min(min(model[,(j+1),1:8]),min(data[1:8,j])),
                max(max(model[,(j+1),1:8]),max(data[1:8,j]))))
    points(117, data[reps[1],j], col=1, pch=19)
    
    for(u in (reps[1]+1):reps[2]){
      points(model[,(j+1),u]~model[,1,u], type="l",
             lwd=2, col=u)
      points(117, data[u,j], col=u, pch=19)
    }
  }
  }
  
  if(gt=="both"|gt=="compare"){
    
  if(number > 9) par(mfrow=c(4,4), mar = c(4,3,3,1))
  if(number <=9) par(mfrow=c(3,3))
  if(number <=6) par(mfrow=c(3,2))
  if(number ==1) par(mfrow=c(1,1))
  
  model2 = t(model[117,-1,])
  data = as.matrix(data)
  for(j in 1:number){
    maxtot = max(max(data[,j],max(model2[,j])))
    ylimit = c(0, maxtot)
    xlimit = c(0, maxtot)
    plot(model2[,j]~data[,j], type="p", col="black", pch=19,
         main=colnames(data)[j], xlab="Data", ylab="Model", ylim = ylimit,
         xlim = xlimit)
    abline(0,1, lwd=2, col="red")
    tt = paste("R2 =", round(summary(lm(model2[,j]~data[,j]))$r.squared,2))
    tt2 = paste("m =", round(lm(model2[,j]~data[,j])$coefficients[2],2))
    text(x = maxtot*0.2, y=maxtot*0.8, label=tt)
    text(x = maxtot*0.2, y=maxtot*0.7, label=tt2)
    
  }

  }
  
  par(mfrow=c(1,1))
  
}

plotcompareF15 <- function(model, data, number =6, reps = c(1,8)){
  
  data = data[,7:12]/data[,1:6]
  
  data[apply(data,2, is.nan)] = 0
  
  model = model[,c(-8,-15),]
  
  model = model[,8:13,]/model[,2:7,]
  
  model[is.nan(model)] =0
  
  ts = 1:117
  
  if(number > 9) par(mfrow=c(4,4), mar = c(4,3,3,1))
  if(number <=9) par(mfrow=c(3,3))
  if(number <=6) par(mfrow=c(3,2))
  
  for(j in 1:number){
    plot(model[,j,reps[1]]~ts, type="l",
         lwd=2, col=1, main=names(data)[j], xlab="Time", ylab="mg/cage",
         ylim=c(min(min(model[,j,reps[1]:reps[2]]),min(data[reps[1]:reps[2],j])),
                max(max(model[,j,reps[1]:reps[2]]),max(data[reps[1]:reps[2],j]))))
    points(117, data[reps[1],j], col=1, pch=19)
    
    for(u in (reps[1]+1):reps[2]){
    points(model[,j,u]~ts, type="l",
           lwd=2, col=u)
    points(117, data[u,j], col=u, pch=19)
    }
  }
  
  if(number > 9) par(mfrow=c(4,4), mar = c(4,3,3,1))
  if(number <=9) par(mfrow=c(3,3))
  if(number <=6) par(mfrow=c(3,2))
  
  model2 = t(model[117,,])
  data = as.matrix(data)
  for(j in 1:number){
    maxtot = max(max(data[,j],max(model2[,j])))
    ylimit = c(0, maxtot)
    xlimit = c(0, maxtot)
    plot(model2[,j]~data[,j], type="p", col="black", pch=19,
         main=colnames(data)[j], xlab="Data", ylab="Model", ylim = ylimit,
         xlim = xlimit)
    abline(0,1, lwd=2, col="red")
    tt = paste("R2 =", round(summary(lm(model2[,j]~data[,j]))$r.squared,2))
    tt2 = paste("m =", round(lm(model2[,j]~data[,j])$coefficients[2],2))
    text(x = maxtot*0.2, y=maxtot*0.8, label=tt)
    text(x = maxtot*0.2, y=maxtot*0.7, label=tt2)
   
  }
  
  par(mfrow=c(1,1))
  
}

plot_figureLP2_scaled <- function(data, model1, model2, fname= paste0("model_scaled_",as.character(Sys.Date()),"LP.pdf")){
  
  #Model 2 should be realistic cycling      
  isopod1 = model1[117,4,] + model1[117,11,]
  litter1 = (model1[1,2,] + model1[1,9,]- model1[117,2,] - model1[117,9,])/117
  #nmin1 = (model1[1,5,] + model1[1,12,] - model1[117,5,]- model1[117,12,])/117
  nmin1 = apply(0.1*model1[91:117,5,]/(10+ model1[91:117,5,]),2,sum)/26
  plant1 = (model1[117,7,] + model1[117,14,])/colSums(model1[117,,])
  
  isopod2 = model2[117,4,] + model2[117,11,]
  litter2 = (model2[1,2,] + model2[1,9,]- model2[117,2,] - model2[117,9,])/117
  #nmin2 = (model2[1,5,] + model2[1,12,] - model2[117,5,]- model2[117,12,])/117
  nmin2 = apply(0.1*model2[91:117,5,]/(10+ model2[91:117,5,]),2,sum)/26
  plant2 = (model2[117,7,] + model2[117,14,])/colSums(model2[117,,])
  
  isopoddata = data$W
  litterdata = (data$L0-data$L)/117
  #nmindata = (data$N-data$N0)/117
  nmindata = data$NS_Oct*1000
  plantdata = data$P/rowSums(data[,1:12])     
  Label = data$Label
  Herbivory = data$Herbivory
  
  toplot = as.data.frame(cbind(isopod1,litter1,nmin1,plant1,isopod2,litter2,nmin2,plant2,
                               isopoddata,litterdata,nmindata,plantdata,Label,Herbivory))
  
  toplot$Label = as.factor(toplot$Label)
  toplot$Herbivory = as.factor(toplot$Herbivory)
  
  note1a = bquote(R^2==.(round(summary(lm(litter1~isopod1 + Label + Herbivory))$r.squared,2)))
  note1b = bquote(R^2==.(round(summary(lm(litter2~isopod2 + Label + Herbivory))$r.squared,2)))
  note1c = bquote(R^2==.(round(summary(lm(litterdata~isopoddata + Label + Herbivory))$r.squared,2)))
  note3a = bquote(R^2==.(round(summary(lm(plant1~litter1))$r.squared,2)))
  note3b = bquote(R^2==.(round(summary(lm(plant2~litter2))$r.squared,2)))
  note3c = bquote(R^2==.(round(summary(lm(plantdata~litterdata))$r.squared,2)))
  
  toplot$Labelp = ifelse(toplot$Label==1,"red", "black")
  toplot$Herbivoryp = ifelse(toplot$Herbivory==1,10,19)
  
  isopodmax = round(max(rbind(isopod1, isopod2, isopoddata)),0) + 2
  decompmax = round(max(rbind(litter1, litter2, litterdata)),0)
  plantmax = max(rbind(plant1, plant2, plantdata)) + 0.02
  
  isopodmin = 0
  decompmin = round(min(rbind(litter1, litter2, litterdata)),0)
  plantmin = min(rbind(plant1, plant2, plantdata))
  
  pdf(fname, width=11, height=6.5)  
  par(mfrow=c(2,3), mar = c(5,5,4,2))
  plot(litter1~isopod1, col=Labelp, pch=Herbivoryp,data=toplot, ylab=expression(Litter~Decomp~(mg[N]~day^-1)),
       xlab=expression(Isopods~(mg[N])), ylim=c(decompmin, decompmax), xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.85, y=decompmax*0.95, note1a, cex=1.5)
  text(x= isopodmax*0.85, y=decompmax*0.9, paste("Isopod",formatC(lm(litter1~isopod1 + Label + Herbivory)$coefficients[2], format = "e", digits = 2)), cex=1)
  text(x= isopodmax*0.85, y=decompmax*0.85, paste("Fertilized",round(lm(litter1~isopod1 + Label + Herbivory)$coefficients[3],2)), cex=1)
  text(x= isopodmax*0.85, y=decompmax*0.8, paste("Herbivory",round(lm(litter1~isopod1 + Label + Herbivory)$coefficients[4],2)), cex=1)
  text(x= isopodmin, y=decompmax*0.98, "d", cex=1.5)
  
  plot(litter2~isopod2, col=Labelp, pch=Herbivoryp,data=toplot, ylab="",
       xlab=expression(Isopods~(mg[N])), ylim=c(decompmin, decompmax), xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.85, y=decompmax*0.95, note1b, cex=1.5)
  text(x= isopodmax*0.85, y=decompmax*0.9, paste("Isopod",formatC(lm(litter2~isopod2 + Label + Herbivory)$coefficients[2], format = "e", digits = 2)), cex=1)
  text(x= isopodmax*0.85, y=decompmax*0.85, paste("Fertilized",round(lm(litter2~isopod2 + Label + Herbivory)$coefficients[3],2)), cex=1)
  text(x= isopodmax*0.85, y=decompmax*0.8, paste("Herbivory",round(lm(litter2~isopod2 + Label + Herbivory)$coefficients[4],2)), cex=1)
  text(x= isopodmin, y=decompmax*0.98, "e", cex=1.5)
  
  plot(litterdata~isopoddata, col=Labelp, pch=Herbivoryp,data=toplot, ylab="",
       xlab=expression(Isopods~(mg[N])), ylim=c(decompmin, decompmax), xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.85, y=decompmax*0.95, note1c, cex=1.5)
  text(x= isopodmax*0.85, y=decompmax*0.9, paste("Isopod",formatC(lm(litterdata~isopoddata + Label + Herbivory)$coefficients[2], format = "e", digits = 2)), cex=1)
  text(x= isopodmax*0.85, y=decompmax*0.85, paste("Fertilized",round(lm(litterdata~isopoddata + Label + Herbivory)$coefficients[3],2)), cex=1)
  text(x= isopodmax*0.85, y=decompmax*0.8, paste("Herbivory",round(lm(litterdata~isopoddata + Label + Herbivory)$coefficients[4],2)), cex=1)
  text(x= isopodmin, y=decompmax*0.98, "f", cex=1.5)
  
  
  plot(plant1~litter1, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Litter~Decomp~(mg[N]~day^-1)), 
       ylab= expression(Proportion~Plants~(mg[N][Plants]~mg[N][Total]^-1)), xlim=c(decompmin, decompmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(y= plantmax*0.95, x=decompmax*0.9, note3a, cex=1.5)
  text(y= plantmax*0.88, x=decompmax*0.9, paste("Decomp.",formatC(lm(plant1~litter1)$coefficients[2], format = "e", digits = 2)), cex=1)
  text(y= plantmax*0.95, x=decompmin*1, "g", cex=1.5)
  
  legend("left", legend = c("Herbivory", "No Herbivory", "Fertilized", "Unfertilized"),
         pch=c(10,19,19,19), col=c("black", "black", "red", "black"), bty="n", cex=1.5)
  
  plot(plant2~litter2, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Litter~Decomp~(mg[N]~day^-1)), 
       ylab= "", xlim=c(decompmin, decompmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(y= plantmax*0.95, x=decompmax*0.9, note3b, cex=1.5)
  text(y= plantmax*0.88, x=decompmax*0.9, paste("Decomp.",formatC(lm(plant2~litter2)$coefficients[2], format = "e", digits = 2)), cex=1)
  text(y= plantmax*0.95, x=decompmin*1, "h", cex=1.5)
  
  
  plot(plantdata~litterdata, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Litter~Decomp~(mg[N]~day^-1)), 
       ylab= "", xlim=c(decompmin, decompmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(y= plantmax*0.95, x=decompmax*0.9, note3c, cex=1.5)
  text(y= plantmax*0.88, x=decompmax*0.9, paste("Decomp.",formatC(lm(plantdata~litterdata)$coefficients[2], format = "e", digits = 2)), cex=1)
  text(y= plantmax*0.95, x=decompmin*1, "i", cex=1.5)
  
  dev.off()
  
}

plot_figureLP2_plantunscaled <- function(data, model1, model2, fname= paste0("model_unscaled_",as.character(Sys.Date()),"LP.pdf")){
  
  #Model 2 should be realistic cycling      
  isopod1 = model1[117,4,] + model1[117,11,]
  litter1 = (model1[1,2,] + model1[1,9,]- model1[117,2,] - model1[117,9,])/117
  #nmin1 = (model1[1,5,] + model1[1,12,] - model1[117,5,]- model1[117,12,])/117
  nmin1 = apply(0.1*model1[91:117,5,]/(10+ model1[91:117,5,]),2,sum)/26
  plant1 = (model1[117,7,] + model1[117,14,])#/colSums(model1[117,,])
  
  isopod2 = model2[117,4,] + model2[117,11,]
  litter2 = (model2[1,2,] + model2[1,9,]- model2[117,2,] - model2[117,9,])/117
  #nmin2 = (model2[1,5,] + model2[1,12,] - model2[117,5,]- model2[117,12,])/117
  nmin2 = apply(0.1*model2[91:117,5,]/(10+ model2[91:117,5,]),2,sum)/26
  plant2 = (model2[117,7,] + model2[117,14,])#/colSums(model2[117,,])
  
  isopoddata = data$W
  litterdata = (data$L0-data$L)/117
  #nmindata = (data$N-data$N0)/117
  nmindata = data$NS_Oct*1000
  plantdata = data$P#/rowSums(data[,1:12])     
  Label = data$Label
  Herbivory = data$Herbivory
  
  toplot = as.data.frame(cbind(isopod1,litter1,nmin1,plant1,isopod2,litter2,nmin2,plant2,
                               isopoddata,litterdata,nmindata,plantdata,Label,Herbivory))
  
  toplot$Label = as.factor(toplot$Label)
  toplot$Herbivory = as.factor(toplot$Herbivory)
  
  note1a = bquote(R^2==.(round(summary(lm(litter1~isopod1 + Label + Herbivory))$r.squared,2)))
  note1b = bquote(R^2==.(round(summary(lm(litter2~isopod2 + Label + Herbivory))$r.squared,2)))
  note1c = bquote(R^2==.(round(summary(lm(litterdata~isopoddata + Label + Herbivory))$r.squared,2)))
  note3a = bquote(R^2==.(round(summary(lm(plant1~litter1))$r.squared,2)))
  note3b = bquote(R^2==.(round(summary(lm(plant2~litter2))$r.squared,2)))
  note3c = bquote(R^2==.(round(summary(lm(plantdata~litterdata))$r.squared,2)))
  
  toplot$Labelp = ifelse(toplot$Label==1,"red", "black")
  toplot$Herbivoryp = ifelse(toplot$Herbivory==1,10,19)
  
  isopodmax = round(max(rbind(isopod1, isopod2, isopoddata)),0) + 2
  decompmax = round(max(rbind(litter1, litter2, litterdata)),0)
  plantmax = max(rbind(plant1, plant2, plantdata)) + 0.02
  
  isopodmin = 0
  decompmin = round(min(rbind(litter1, litter2, litterdata)),0)
  plantmin = min(rbind(plant1, plant2, plantdata))
  
  pdf(fname, width=11, height=4.5)  
  par(mfrow=c(1,3), mar = c(5,5,4,2))
  
  #expression(Proportion~Plants~(mg[N][Plants]~mg[N][Total]^-1))
  plot(plant1~litter1, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Litter~Decomp~(mg[N]~day^-1)),main="Tight Cycling", 
       ylab= expression(Plants~(mg[N])), xlim=c(decompmin, decompmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(y= plantmax*0.95, x=decompmax*0.9, note3a, cex=1.5)
  text(y= plantmax*0.88, x=decompmax*0.9, paste("Decomp.",formatC(lm(plant1~litter1)$coefficients[2], format = "e", digits = 2)), cex=1)
  text(y= plantmax*0.95, x=decompmin*1, "a", cex=1.5)
  
  legend("left", legend = c("Herbivory", "No Herbivory", "Fertilized", "Unfertilized"),
         pch=c(10,19,19,19), col=c("black", "black", "red", "black"), bty="n", cex=1.5)
  
  plot(plant2~litter2, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Litter~Decomp~(mg[N]~day^-1)), main="Diffuse Cycling",
       ylab= "", xlim=c(decompmin, decompmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(y= plantmax*0.95, x=decompmax*0.9, note3b, cex=1.5)
  text(y= plantmax*0.88, x=decompmax*0.9, paste("Decomp.",formatC(lm(plant2~litter2)$coefficients[2], format = "e", digits = 2)), cex=1)
  text(y= plantmax*0.95, x=decompmin*1, "b", cex=1.5)
  
  
  plot(plantdata~litterdata, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Litter~Decomp~(mg[N]~day^-1)), main="Field Data",
       ylab= "", xlim=c(decompmin, decompmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(y= plantmax*0.95, x=decompmax*0.9, note3c, cex=1.5)
  text(y= plantmax*0.88, x=decompmax*0.9, paste("Decomp.",formatC(lm(plantdata~litterdata)$coefficients[2], format = "e", digits = 2)), cex=1)
  text(y= plantmax*0.95, x=decompmin*1, "c", cex=1.5)
  
  dev.off()
  
}

plot_figureISO2 <- function(data, model1, model2, fname= paste0("model_",as.character(Sys.Date()),"ISO.pdf")){
  
  #Model 2 should be realistic cycling      
  isopod1 = model1[117,4,] + model1[117,11,]
  litter1 = (model1[1,2,] + model1[1,9,]- model1[117,2,] - model1[117,9,])/117
  #nmin1 = (model1[1,5,] + model1[1,12,] - model1[117,5,]- model1[117,12,])/117
  nmin1 = apply(0.1*model1[91:117,5,]/(10+ model1[91:117,5,]),2,sum)/26
  plant1 = model1[117,7,] + model1[117,14,]
  
  isopod2 = model2[117,4,] + model2[117,11,]
  litter2 = (model2[1,2,] + model2[1,9,]- model2[117,2,] - model2[117,9,])/117
  #nmin2 = (model2[1,5,] + model2[1,12,] - model2[117,5,]- model2[117,12,])/117
  nmin2 = apply(0.1*model2[91:117,5,]/(10+ model2[91:117,5,]),2,sum)/26
  plant2 = model2[117,7,] + model2[117,14,]
  
  isopoddata = data$W
  litterdata = (data$L0-data$L)/117
  #nmindata = (data$N-data$N0)/117
  nmindata = data$NS_Oct*1000
  plantdata = data$P        
  Label = data$Label
  Herbivory = data$Herbivory
  
  toplot = as.data.frame(cbind(isopod1,litter1,nmin1,plant1,isopod2,litter2,nmin2,plant2,
                               isopoddata,litterdata,nmindata,plantdata,Label,Herbivory))
  
  toplot$Label = as.factor(toplot$Label)
  toplot$Herbivory = as.factor(toplot$Herbivory)
  
  note1a = bquote(R^2==.(round(summary(lm(litter1~isopod1*Label*Herbivory))$r.squared,2)))
  note1b = bquote(R^2==.(round(summary(lm(litter2~isopod2*Label*Herbivory))$r.squared,2)))
  note1c = bquote(R^2==.(round(summary(lm(litterdata~isopoddata*Label*Herbivory))$r.squared,2)))
  note2a = bquote(R^2==.(round(summary(lm(nmin1~isopod1*Label*Herbivory))$r.squared,2)))
  note2b = bquote(R^2==.(round(summary(lm(nmin2~isopoddata*Label*Herbivory))$r.squared,2)))
  note2c = bquote(R^2==.(round(summary(lm(nmindata~litterdata*Label*Herbivory))$r.squared,2)))
  note3a = bquote(R^2==.(round(summary(lm(plant1~isopod1*Label*Herbivory))$r.squared,2)))
  note3b = bquote(R^2==.(round(summary(lm(plant2~isopod2*Label*Herbivory))$r.squared,2)))
  note3c = bquote(R^2==.(round(summary(lm(plantdata~isopoddata*Label*Herbivory))$r.squared,2)))
  
  toplot$Labelp = ifelse(toplot$Label==1,"red", "black")
  toplot$Herbivoryp = ifelse(toplot$Herbivory==1,10,19)
  
  isopodmax = round(max(rbind(isopod1, isopod2, isopoddata)),0)
  decompmax = round(max(rbind(litter1, litter2, litterdata)),0)
  plantmax = round(max(rbind(plant1, plant2, plantdata)),0) + 1
  nminmax = round(max(rbind(nmin1, nmin2)),2)
  
  isopodmin = 0
  decompmin = round(min(rbind(litter1, litter2, litterdata)),0)
  plantmin = round(min(rbind(plant1, plant2, plantdata)),0) + 1
  nminmin = round(min(rbind(nmin1, nmin2)),2)
  
  pdf(fname, width=11, height=7.5)  
  par(mfrow=c(3,3), mar = c(5,5,4,2))
  
  #Litter
  plot(litter1~isopod1, col=Labelp, pch=Herbivoryp,data=toplot, ylab=expression(Litter~Decomp~(mg[N]~day^-1)), 
       main="a. Tight Cycling",
       xlab="", ylim=c(decompmin, decompmax), xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.9, y=decompmax*0.95, note1a, cex=1.5)
  
  plot(litter2~isopod2, col=Labelp, pch=Herbivoryp,data=toplot, ylab=expression(Litter~Decomp~(mg[N]~day^-1)), 
       main="b. Realistic Cycling",
       xlab="", ylim=c(decompmin, decompmax), xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.9, y=decompmax*0.95, note1b, cex=1.5)
  
  plot(litterdata~isopoddata, col=Labelp, pch=Herbivoryp,data=toplot, ylab=expression(Litter~Decomp~(mg[N]~day^-1)), 
       main="c. Field Data",
       xlab="", ylim=c(decompmin, decompmax), xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.9, y=decompmax*0.95, note1c, cex=1.5)
  
  #Nmin
  plot(nmin1~isopod1, col=Labelp, pch=Herbivoryp,data=toplot, ylab=expression(Soil~Min~(mg[N]~day^-1)), 
       xlab="", xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.9, y=0.0014, note2a, cex=1.5)
  
  plot(nmin2~isopod2, col=Labelp, pch=Herbivoryp,data=toplot, ylab=expression(Soil~Min~(mg[N]~day^-1)), 
       xlab="", xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.9, y=0.0272, note2b, cex=1.5)
  
  plot(nmindata~isopoddata, col=Labelp, pch=Herbivoryp,data=toplot, ylab=expression(Soil~Min~(mg[N]~day^-1~cm^-2)), 
       xlab="", xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.9, y=0.2, note2c, cex=1.5)
  
  #Plant same scales
  # plot(plant1~isopod1, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Isopods~(mg[N])), 
  #      ylab= expression(Plants~(mg[N])), xlim=c(isopodmin, isopodmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  # text(y= plantmax*0.85, x=isopodmax*0.9, note3a, cex=1.5)
  # legend("topleft", legend = c("Herbivory", "No Herbivory", "Fertilized", "Unfertilized"),
  #        pch=c(10,19,19,19), col=c("black", "black", "red", "black"), bty="n", cex=1.5)
  # 
  # plot(plant2~isopod2, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Isopods~(mg[N])), 
  #      ylab= expression(Plants~(mg[N])), xlim=c(isopodmin, isopodmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  # text(y= plantmax*0.85, x=isopodmax*0.9, note3b, cex=1.5)
  # 
  # plot(plantdata~isopoddata, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Isopods~(mg[N])), 
  #      ylab= expression(Plants~(mg[N])), xlim=c(isopodmin, isopodmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  # text(y= plantmax*0.85, x=isopodmax*0.9, note3c, cex=1.5)
  
  #Plant different scales
  plot(plant1~isopod1, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Isopods~(mg[N])), 
       ylab= expression(Plants~(mg[N])), xlim=c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(y= 100, x=isopodmax*0.9, note3a, cex=1.5)
  
  plot(plant2~isopod2, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Isopods~(mg[N])), 
       ylab= expression(Plants~(mg[N])), xlim=c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(y= 240, x=isopodmax*0.9, note3b, cex=1.5)
  
  plot(plantdata~isopoddata, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Isopods~(mg[N])), 
       ylab= expression(Plants~(mg[N])), xlim=c(isopodmin, isopodmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(y= plantmax*0.85, x=isopodmax*0.9, note3c, cex=1.5)
  
  dev.off()
  
}

plot_figureISO_realcompare <- function(data, model1, model2, fname= paste0("model_",as.character(Sys.Date()),"ISO.pdf")){
  
  #Model 1 should be actual inputs
  #Model 2 should be homogenized inputs
  isopod1 = model1[117,4,] + model1[117,11,]
  litter1 = (model1[1,2,] + model1[1,9,]- model1[117,2,] - model1[117,9,])/117
  #nmin1 = (model1[1,5,] + model1[1,12,] - model1[117,5,]- model1[117,12,])/117
  nmin1 = apply(0.1*model1[91:117,5,]/(10+ model1[91:117,5,]),2,sum)/26
  soiln1 = model1[117,5,] + model1[117,12,]
  plant1 = model1[117,7,] + model1[117,14,]
  
  isopod2 = model2[117,4,] + model2[117,11,]
  litter2 = (model2[1,2,] + model2[1,9,]- model2[117,2,] - model2[117,9,])/117
  #nmin2 = (model2[1,5,] + model2[1,12,] - model2[117,5,]- model2[117,12,])/117
  nmin2 = apply(0.1*model2[91:117,5,]/(10+ model2[91:117,5,]),2,sum)/26
  soiln2 = model2[117,5,] + model2[117,12,]
  plant2 = model2[117,7,] + model2[117,14,]
  
  
  toplot = as.data.frame(cbind(isopod1,litter1,nmin1,soiln1,plant1,isopod2,litter2,nmin2,soiln2,plant2,Label,Herbivory))
  
  toplot$Label = as.factor(toplot$Label)
  toplot$Herbivory = as.factor(toplot$Herbivory)
  
  note1a = bquote(R^2==.(round(summary(lm(litter1~isopod1*Label*Herbivory))$r.squared,2)))
  note1b = bquote(R^2==.(round(summary(lm(litter2~isopod2*Label*Herbivory))$r.squared,2)))
  note2a = bquote(R^2==.(round(summary(lm(nmin1~isopod1*Label*Herbivory))$r.squared,2)))
  note2b = bquote(R^2==.(round(summary(lm(nmin2~isopod2*Label*Herbivory))$r.squared,2)))
  note3a = bquote(R^2==.(round(summary(lm(plant1~isopod1*Label*Herbivory))$r.squared,2)))
  note3b = bquote(R^2==.(round(summary(lm(plant2~isopod2*Label*Herbivory))$r.squared,2)))
  
  note4a = bquote(R^2==.(round(summary(lm(soiln1~isopod1*Label*Herbivory))$r.squared,2)))
  note4b = bquote(R^2==.(round(summary(lm(soiln2~isopod2*Label*Herbivory))$r.squared,2)))
  
  toplot$Labelp = ifelse(toplot$Label==1,"red", "black")
  toplot$Herbivoryp = ifelse(toplot$Herbivory==1,10,19)
  
  isopodmax = round(max(rbind(isopod1, isopod2)),0)
  decompmax = round(max(rbind(litter1, litter2)),0)
  plantmax = round(max(rbind(plant1, plant2)),0) + 1
  nminmax = round(max(rbind(nmin1, nmin2)),2)
  soilnmax = round(max(rbind(soiln1, soiln2)),2)
  
  isopodmin = 0
  decompmin = round(min(rbind(litter1, litter2)),0)
  plantmin = round(min(rbind(plant1, plant2)),0) + 1
  nminmin = round(min(rbind(nmin1, nmin2)),2)
  soilnmin = round(min(rbind(soiln1, soiln2)),2)
  
  pdf(fname, width=6, height=7.5)  
  par(mfrow=c(4,2), mar = c(4,5,3,1), oma=c(3,0,0,0))
  
  #Litter
  plot(litter1~isopod1, col=Labelp, pch=Herbivoryp,data=toplot, ylab=expression(Litter~Decomp~(mg[N]~day^-1)), 
       main="a. Heterogenous Start",
       xlab="", ylim=c(decompmin, decompmax), xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.9, y=decompmax*0.95, note1a, cex=1)
  
  plot(litter2~isopod2, col=Labelp, pch=Herbivoryp,data=toplot, ylab="", 
       main="b. Homogenous Start",
       xlab="", ylim=c(decompmin, decompmax), xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.9, y=decompmax*0.95, note1b, cex=1)
  
  # plot(litterdata~isopoddata, col=Labelp, pch=Herbivoryp,data=toplot, ylab=expression(Litter~Decomp~(mg[N]~day^-1)), 
  #      main="c. Field Data",
  #      xlab="", ylim=c(decompmin, decompmax), xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  # text(x= isopodmax*0.9, y=decompmax*0.95, note1c, cex=1.5)
  
  #Nmin
  plot(nmin1~isopod1, col=Labelp, pch=Herbivoryp,data=toplot, ylab=expression(Soil~Min~(mg[N]~day^-1)), 
       xlab="", ylim=c(nminmin, nminmax), xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.9, y=0.025, note2a, cex=1)
  
  plot(nmin2~isopod2, col=Labelp, pch=Herbivoryp,data=toplot, ylab="", 
       xlab="", ylim=c(nminmin, nminmax), xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.9, y=0.045, note2b, cex=1)
  
  #Soil N
  plot(soiln1~isopod1, col=Labelp, pch=Herbivoryp,data=toplot, ylab=expression(Soil~N~(mg[N])), 
       xlab="", ylim=c(soilnmin, soilnmax), xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.9, y=3, note4a, cex=1)
  
  plot(soiln2~isopod2, col=Labelp, pch=Herbivoryp,data=toplot, ylab="", 
       xlab="", ylim=c(soilnmin, soilnmax), xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(x= isopodmax*0.9, y=7, note4b, cex=1)
  
  # plot(nmindata~isopoddata, col=Labelp, pch=Herbivoryp,data=toplot, ylab=expression(Soil~Min~(mg[N]~day^-1~cm^-2)), 
  #      xlab="", xlim= c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  # text(x= isopodmax*0.9, y=0.2, note2c, cex=1.5)
  
  #Plant same scales
  # plot(plant1~isopod1, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Isopods~(mg[N])), 
  #      ylab= expression(Plants~(mg[N])), xlim=c(isopodmin, isopodmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  # text(y= plantmax*0.85, x=isopodmax*0.9, note3a, cex=1.5)
  # legend("topleft", legend = c("Herbivory", "No Herbivory", "Fertilized", "Unfertilized"),
  #        pch=c(10,19,19,19), col=c("black", "black", "red", "black"), bty="n", cex=1.5)
  # 
  # plot(plant2~isopod2, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Isopods~(mg[N])), 
  #      ylab= expression(Plants~(mg[N])), xlim=c(isopodmin, isopodmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  # text(y= plantmax*0.85, x=isopodmax*0.9, note3b, cex=1.5)
  # 
  # plot(plantdata~isopoddata, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Isopods~(mg[N])), 
  #      ylab= expression(Plants~(mg[N])), xlim=c(isopodmin, isopodmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  # text(y= plantmax*0.85, x=isopodmax*0.9, note3c, cex=1.5)
  
  #Plant different scales
  plot(plant1~isopod1, col=Labelp, pch=Herbivoryp,data=toplot, xlab="", 
       ylab= expression(Plants~(mg[N])), ylim=c(plantmin, plantmax), xlim=c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(y= 400, x=isopodmax*0.9, note3a, cex=1)
  
  plot(plant2~isopod2, col=Labelp, pch=Herbivoryp,data=toplot, xlab="", 
       ylab= "", ylim=c(plantmin, plantmax), xlim=c(isopodmin, isopodmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  text(y= 400, x=isopodmax*0.9, note3b, cex=1)
  
  # plot(plantdata~isopoddata, col=Labelp, pch=Herbivoryp,data=toplot, xlab=expression(Isopods~(mg[N])), 
  #      ylab= expression(Plants~(mg[N])), xlim=c(isopodmin, isopodmax), ylim= c(plantmin, plantmax), cex=2, cex.axis=1.2, cex.lab=1.2)
  # text(y= plantmax*0.85, x=isopodmax*0.9, note3c, cex=1.5)
  
  mtext(text=expression(Isopods~(mg[N])), side=1, line=0, outer=TRUE)
  
  dev.off()
  
}