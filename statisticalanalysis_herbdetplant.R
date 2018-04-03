# Statistical Analysis of Data from 

#load data
library(readr)
ch2 <- read_csv("Chapter2_complieddata.csv", 
                col_types = cols(Herbivory = col_factor(levels = c("0", 
                                                                   "1")), Label = col_factor(levels = c("0", 
                                                                                                        "1")), Treatment = col_factor(levels = c("1", 
                                                                                                                                                 "2", "3", "4"))))
perN <- read_csv("perN_balance.csv")
perC <- read_csv("perC_balance.csv")

# Calculate the decomposition of litter N
Ldecomptime = ifelse(ch2$Label==0|ch2$Plot <9, 289-170, 290-170)
ch2["LitterN2"] = ch2$Litter*perN$Litter/100
ch2["LNDecomp"] = 1000*(7*(ch2$perN/100) - ch2$LitterN2)/Ldecomptime # mg-N day-1

#Subset out the control plots

ch2exp = subset(ch2, Plot <41)

head(ch2exp)


# Run model selection -----------------------------------------------------

# create center data frame for calculating vif
library(car)
variablestocenter = c("Isopods","N_Sol_June", "SIR_June", "Worms","NE_June", "NS_June")

ch2exp_center = ch2exp

for(i in 1:length(variablestocenter)) {
  
  ch2exp_center[variablestocenter[i]] = ch2exp_center[variablestocenter[i]] - colMeans(ch2exp_center[variablestocenter[i]], na.rm=T)
}

#*********************************************
# Plants
#*********************************************

vif(lm(Plant~(Herbivory+Label+Isopods+N_Sol_June)^2 +SIR_June+Worms + NE_June, data=ch2exp_center))<10

P1 = lm(Plant~(Herbivory+Label+Isopods+N_Sol_June)^2 + SIR_June + Worms + NE_June, data=ch2exp); summary(P1)
P2 = update(P1, .~. -Label:N_Sol_June); summary(P2)
P3 = update(P2, .~. -Label:Isopods); summary(P3)
P4 = update(P3, .~. -Herbivory:Label); summary(P4)
P5 = update(P4, .~. -NE_June); summary(P5)
P6 = update(P5, .~. -Herbivory:Isopods); summary(P6)
P6 = update(P6, .~. -Label); summary(P6)

# check if removal of high isopods and worm mesocosms matters
summary(lm(formula(P6), data=subset(ch2exp, Isopods <20 & Worms <20)))

#*********************************************
# Litter Decomposition
#*********************************************

vif(lm(LNDecomp~(Herbivory+Label+Isopods)^2 +N_Sol_June+SIR_June+Worms + NE_June, data=ch2exp_center))<10

LND1 = lm(LNDecomp~(Herbivory+Label+Isopods)^2 +N_Sol_June+SIR_June+Worms + NE_June, data=ch2exp); summary(LND1)
LND2 = update(LND1, .~. -NE_June); summary(LND2)
LND3 = update(LND2, .~. -N_Sol_June); summary(LND3)
LND4 = update(LND3, .~. -Worms); summary(LND4)
LND5 = update(LND4, .~. -Herbivory:Isopods); summary(LND5)
LND6 = update(LND5, .~. -Herbivory:Label); summary(LND6)

# check if removal of high isopods and worm mesocosms matters
summary(lm(formula(LND6), data=subset(ch2exp, Isopods <20 & Worms <20)))

#*********************************************
# Extractable Nitrogen
#*********************************************
NE1 = lm(NE_Oct~(Herbivory+Label+Isopods+ NE_June)^2+N_Sol_June+SIR_June+Worms, data=ch2exp); summary(NE1)
NE2 = update(NE1, .~. -Herbivory:Isopods); summary(NE2)
NE3 = update(NE2, .~. -SIR_June); summary(NE3)
NE4 = update(NE3, .~. -Herbivory:Label); summary(NE4)
NE5 = update(NE4, .~. -Label:Isopods); summary(NE5)
NE6 = update(NE5, .~. -Herbivory:NE_June); summary(NE6)
NE7 = update(NE6, .~. -Isopods:NE_June); summary(NE7)
NE8 = update(NE7, .~. -Isopods); summary(NE8)
NE9 = update(NE8, .~. -Herbivory); summary(NE9)
NE10 = update(NE9, .~. -Label:NE_June); summary(NE10)
NE11 = update(NE10, .~. -Worms); summary(NE11)
NE12 = update(NE11, .~. -NE_June); summary(NE12)
NE13 = update(NE12, .~. -Herbivory:SIR_June); summary(NE13)

# check if removal of high isopods and worm mesocosms matters
summary(lm(formula(NE13), data=subset(ch2exp, Isopods <20 & Worms <20)))

#*********************************************
# Ion exchange membrane Nitrogen
#*********************************************

vif(lm(NS_Oct~(Herbivory+Label+Isopods+ NS_June)^2+N_Sol_June+SIR_June+Worms, data=ch2exp_center))<10

NS1 = lm(NS_Oct~(Herbivory+Label+Isopods+ NS_June)^2+N_Sol_June+SIR_June+Worms, data=ch2exp); summary(NS1)
NS2 = update(NS1, .~. -Herbivory:Isopods); summary(NS2)
NS3 = update(NS2, .~. -Label:Isopods); summary(NS3)
NS4 = update(NS3, .~. -Herbivory:Label); summary(NS4)
NS5 = update(NS4, .~. -N_Sol_June); summary(NS5)
NS6 = update(NS5, .~. -SIR_June); summary(NS6)
NS7 = update(NS6, .~. -Herbivory:NS_June); summary(NS7)
NS8 = update(NS7, .~. -Label:NS_June); summary(NS8)
NS9 = update(NS8, .~. -Label); summary(NS9)
NS10 = update(NS9, .~. -Isopods:NS_June ); summary(NS10)
NS11 = update(NS10, .~. -NS_June); summary(NS11)
NS12 = update(NS11, .~. -Isopods); summary(NS12)
NS13 = update(NS12, .~. -Worms); summary(NS13)
# check if removal of high isopods and worm mesocosms matters
summary(lm(formula(NS13), data=subset(ch2exp, Isopods <20 & Worms <20)))

#*********************************************
# Roots
#*********************************************

vif(lm(Roots~(Herbivory+Label+Isopods+N_Sol_June)^2 +SIR_June+Worms + NE_June, data=ch2exp_center))<10

R1 = lm(Roots~(Herbivory+Label+Isopods+N_Sol_June)^2 + SIR_June + Worms + NE_June, data=ch2exp); summary(R1)
R2 = update(R1, .~. -Herbivory:N_Sol_June); summary(R2)
R3 = update(R2, .~. -Isopods:N_Sol_June); summary(R3)
R4 = update(R3, .~. -Label:Isopods); summary(R4)
R5 = update(R4, .~. -Worms); summary(R5)
R6 = update(R5, .~. -SIR_June); summary(R6)
R7 = update(R6, .~. -NE_June); summary(R7)
R8 = update(R7, .~. -Herbivory:Label); summary(R8)
R9 = update(R8, .~. -Label:N_Sol_June); summary(R9)
R10 = update(R9, .~. -N_Sol_June); summary(R10)
R11 = update(R10, .~. -Isopods:Herbivory); summary(R11)
R12 = update(R11, .~. -Herbivory); summary(R12)
R13 = update(R12, .~. -Isopods); summary(R13)

# check if removal of high isopods and worm mesocosms matters
summary(lm(formula(R13), data=subset(ch2exp, Isopods <20 & Worms <20)))

#*********************************************
# Total Plant Biomass
#*********************************************

vif(lm(Plant + Roots~(Herbivory+Label+Isopods+N_Sol_June)^2 +SIR_June+Worms + NE_June, data=ch2exp_center))<10

TPR1 = lm(Plant + Roots~(Herbivory+Label+Isopods+N_Sol_June)^2 + SIR_June + Worms + NE_June, data=ch2exp); summary(TPR1)
TPR2 = update(TPR1, .~. -Herbivory:N_Sol_June); summary(TPR2)
TPR3 = update(TPR2, .~. -Isopods:N_Sol_June); summary(TPR3)
TPR4 = update(TPR3, .~. -Label:Isopods); summary(TPR4)
TPR5 = update(TPR4, .~. -Worms); summary(TPR5)
TPR6 = update(TPR5, .~. -SIR_June); summary(TPR6)
TPR7 = update(TPR6, .~. -NE_June); summary(TPR7)
TPR8 = update(TPR7, .~. -Herbivory:Label); summary(TPR8)
TPR9 = update(TPR8, .~. -Label:N_Sol_June); summary(TPR9)
TPR10 = update(TPR9, .~. -N_Sol_June); summary(TPR10)
TPR11 = update(TPR10, .~. -Isopods:Herbivory); summary(TPR11)
TPR12 = update(TPR11, .~. -Herbivory); summary(TPR12)
TPR13 = update(TPR12, .~. -Isopods); summary(TPR13)

# check if removal of high isopods and worm mesocosms matters
summary(lm(formula(TPR13), data=subset(ch2exp, Isopods <20 & Worms <20)))

#*********************************************
# Forbs
#*********************************************

vif(lm(Forbs~(Herbivory+Label+Isopods+N_Sol_June)^2 +SIR_June+Worms + NE_June, data=ch2exp_center))<10

F1 = lm(Forbs~(Herbivory+Label+Isopods+N_Sol_June)^2 + SIR_June + Worms + NE_June, data=ch2exp); summary(F1)
F2 = update(F1, .~. -Label:N_Sol_June); summary(F2)
F3 = update(F2, .~. -Herbivory:Label); summary(F3)
F4 = update(F3, .~. -NE_June); summary(F4)
F5 = update(F4, .~. -SIR_June); summary(F5)
F6 = update(F5, .~. -Label:Isopods); summary(F6)
F7 = update(F6, .~. -Herbivory:Isopods); summary(F7)
F8 = update(F7, .~. -Worms); summary(F8)
F9 = update(F8, .~. -Isopods:N_Sol_June); summary(F9)
F10 = update(F9, .~. -Isopods); summary(F10)

# check if removal of high isopods and worm mesocosms matters
summary(lm(formula(F10), data=subset(ch2exp, Isopods <20 & Worms <20)))

#*********************************************
# Grass
#*********************************************

vif(lm(Grass~(Herbivory+Label+Isopods+N_Sol_June)^2 +SIR_June+Worms + NE_June, data=ch2exp_center))<10

G1 = lm(Grass~(Herbivory+Label+Isopods+N_Sol_June)^2 + SIR_June + Worms + NE_June, data=ch2exp); summary(G1)
G2 = update(G1, .~. -Herbivory:Isopods); summary(G2)
G3 = update(G2, .~. -Herbivory:N_Sol_June); summary(G3)
G4 = update(G3, .~. -Label:N_Sol_June); summary(G4)
G5 = update(G4, .~. -Label:Isopods); summary(G5)
G6 = update(G5, .~. -Herbivory:Label); summary(G6)
G7 = update(G6, .~. -Isopods:N_Sol_June); summary(G7)
G8 = update(G7, .~. -Isopods); summary(G8)
G9 = update(G8, .~. -N_Sol_June); summary(G9)
G10 = update(G9, .~. -NE_June); summary(G10)
G11 = update(G10, .~. -Herbivory); summary(G11)
G12 = update(G11, .~. -Worms); summary(G12)

# check if removal of high isopods and worm mesocosms matters
summary(lm(formula(G12), data=subset(ch2exp, Isopods <20 & Worms <20)))


# Plot Model Results ------------------------------------------------------

#*********************************************
# Create predictions
#*********************************************

for(i in 0:1){
  for(j in 0:1){
    isorange = range(subset(ch2exp, Herbivory==i & Label==j)$Isopods)
    isoseq = seq(isorange[1], isorange[2], by = 1)
    assign(paste0("Treat",i,j),cbind(isoseq, rep(i, length(isoseq)), rep(j, length(isoseq))))
  }
}

newdata = as.data.frame(rbind(Treat00, Treat01, Treat10, Treat11))
names(newdata) = c("Isopods", "Herbivory", "Label")
newdata["Worms"] = rep(mean(ch2exp$Worms), times=dim(newdata)[1])
newdata["SIR_June"] = rep(mean(ch2exp$SIR_June, na.rm=T), times=dim(newdata)[1])
newdata["N_Sol_June"] = rep(2, times=dim(newdata)[1])
newdata["Herbivory"] = as.factor(newdata$Herbivory)
newdata["Label"] = as.factor(newdata$Label)
newdata["NE_June"] = rep(mean(ch2exp$NE_June), times=dim(newdata)[1])
newdata["NS_June"] = rep(mean(ch2exp$NS_June), times=dim(newdata)[1])

#New Model
newdata["Npredict"] = predict(P6, newdata)
newdata["Npredictse"] = predict(P6, newdata, se.fit=T)$se.fit

newdata["NpredictNE"] = predict(NE13, newdata)
newdata["NpredictseNE"] = predict(NE13, newdata, se.fit=T)$se.fit

newdata["NpredictL"] = predict(LND6, newdata)
newdata["NpredictseL"] = predict(LND6, newdata, se.fit=T)$se.fit

#*********************************************
# Plot Results
#*********************************************

pdf("Figure2.pdf", width=8, height=6)

par(oma=c(3,0,0,0),mar=c(3,5,2,2),mfrow=c(2,2), cex.lab=1.2)

#Litter
plot(LNDecomp~Isopods, data=subset(ch2exp, Herbivory==0), pch=19, col=Label,
     xlab="Isopods (#)", ylab=expression(Litter~Decomp~(mg[N]~day^-1)), main="No Herbivory", type="n",
     ylim=c(1,2), xlim=c(0,20))
text(x=0, y=1.9, label="a", cex=2)
subofdata = subset(newdata, Label==0 & Herbivory==0)
lines(subofdata$Isopods,
      subofdata$NpredictL, cex=2, lty=1)
polygon(c(subofdata$Isopods,rev(subofdata$Isopods)),
        c(subofdata$NpredictL+subofdata$NpredictseL,rev(subofdata$NpredictL-subofdata$NpredictseL)),
        col=scales::alpha("grey",.5))


subofdata = subset(newdata, Label==1 & Herbivory==0)
lines(subofdata$Isopods,
      subofdata$NpredictL, cex=2, lty=1, col="red")
polygon(c(subofdata$Isopods,rev(subofdata$Isopods)),
        c(subofdata$NpredictL+subofdata$NpredictseL,rev(subofdata$NpredictL-subofdata$NpredictseL)),
        col=scales::alpha("red",.5))
points(subset(ch2exp, Herbivory==0)$Isopods, subset(ch2exp, Herbivory==0)$LNDecomp, pch=19,col=subset(ch2exp, Herbivory==0)$Label, cex=1.5)

plot(LNDecomp~Isopods, data=subset(ch2exp, Herbivory==1), pch=19, col=Label,
     xlab="Isopods (#)", ylab="", main="Herbivory", type="n",
     ylim=c(1,2), xlim=c(0,20))
text(x=0, y=1.9, label="b", cex=2)
subofdata = subset(newdata, Label==0 & Herbivory==1)
lines(subofdata$Isopods,
      subofdata$NpredictL, cex=2, lty=2)
polygon(c(subofdata$Isopods,rev(subofdata$Isopods)),
        c(subofdata$NpredictL+subofdata$NpredictseL,rev(subofdata$NpredictL-subofdata$NpredictseL)),
        col=scales::alpha("grey",.5))

subofdata = subset(newdata, Label==1 & Herbivory==1)
lines(subofdata$Isopods,
      subofdata$NpredictL, cex=2, lty=2, col="red")
polygon(c(subofdata$Isopods,rev(subofdata$Isopods)),
        c(subofdata$NpredictL+subofdata$NpredictseL,rev(subofdata$NpredictL-subofdata$NpredictseL)),
        col=scales::alpha("red",.5))
points(subset(ch2exp, Herbivory==1)$Isopods, subset(ch2exp, Herbivory==1)$LNDecomp, pch=10,col=subset(ch2exp, Herbivory==1)$Label, cex=1.5)

#NE Oct
plot(NE_Oct~Isopods, data=subset(ch2exp, Herbivory==0), pch=19, col=Label,
     xlab="Isopods (#)", ylab=expression(Soil~N~(mu*g[N]~g[DMES]^-1)), type="n",
     ylim=c(0,5), xlim=c(0,20))
text(x=0, y=4.5, label="c", cex=2)
subofdata = subset(newdata, Label==0 & Herbivory==0)
lines(subofdata$Isopods,
      subofdata$NpredictNE, cex=2, lty=1)
polygon(c(subofdata$Isopods,rev(subofdata$Isopods)),
        c(subofdata$NpredictNE+subofdata$NpredictseNE,rev(subofdata$NpredictNE-subofdata$NpredictseNE)),
        col=scales::alpha("grey",.5))

subofdata = subset(newdata, Label==1 & Herbivory==0)
lines(subofdata$Isopods,
      subofdata$NpredictNE, cex=2, lty=1, col="red")
polygon(c(subofdata$Isopods,rev(subofdata$Isopods)),
        c(subofdata$NpredictNE+subofdata$NpredictseNE,rev(subofdata$NpredictNE-subofdata$NpredictseNE)),
        col=scales::alpha("red",.5))
points(subset(ch2exp, Herbivory==0)$Isopods, subset(ch2exp, Herbivory==0)$NE_Oct, pch=19, col=subset(ch2exp, Herbivory==0)$Label, cex=1.5)

plot(NE_Oct~Isopods, data=subset(ch2exp, Herbivory==1), pch=19, col=Label,
     xlab="Isopods (#)", ylab="", type="n",
     ylim=c(0,5), xlim=c(0,20))
text(x=0, y=4.5, label="d", cex=2)
subofdata = subset(newdata, Label==0 & Herbivory==1)
lines(subofdata$Isopods,
      subofdata$NpredictNE, cex=2, lty=2)
polygon(c(subofdata$Isopods,rev(subofdata$Isopods)),
        c(subofdata$NpredictNE+subofdata$NpredictseNE,rev(subofdata$NpredictNE-subofdata$NpredictseNE)),
        col=scales::alpha("grey",.5))

subofdata = subset(newdata, Label==1 & Herbivory==1)
lines(subofdata$Isopods,
      subofdata$NpredictNE, cex=2, lty=2, col="red")
polygon(c(subofdata$Isopods,rev(subofdata$Isopods)),
        c(subofdata$NpredictNE+subofdata$NpredictseNE,rev(subofdata$NpredictNE-subofdata$NpredictseNE)),
        col=scales::alpha("red",.5))
points(subset(ch2exp, Herbivory==1)$Isopods, subset(ch2exp, Herbivory==1)$NE_Oct, pch=10, col=subset(ch2exp, Herbivory==1)$Label, cex=1.5)

mtext(text="Isopods (#)",side=1,line=0,outer=TRUE)
dev.off()

pdf("Figure3.pdf", width=8, height=6)

par(oma=c(1,1,0,0),mar=c(3,5,2,2),mfrow=c(2,2), cex.lab=1.2)

#Plants
plot(Plant~Isopods, data=subset(ch2exp, Herbivory==0), pch=19, col=Label,
     xlab="Isopods (#)", ylab="", type="n",main="No Herbivory",
     ylim=c(0,14), xlim=c(0,20))
text(x=0, y=13.5, label="a", cex=2)
subofdata = subset(newdata, Label==0 & Herbivory==0)
lines(subofdata$Isopods,
      subofdata$Npredict, cex=2, lty=1)
polygon(c(subofdata$Isopods,rev(subofdata$Isopods)),
        c(subofdata$Npredict+subofdata$Npredictse,rev(subofdata$Npredict-subofdata$Npredictse)),
        col=scales::alpha("grey",.5))

subofdata = subset(newdata, Label==1 & Herbivory==0)
lines(subofdata$Isopods,
      subofdata$Npredict, cex=2, lty=1, col="red")
polygon(c(subofdata$Isopods,rev(subofdata$Isopods)),
        c(subofdata$Npredict+subofdata$Npredictse,rev(subofdata$Npredict-subofdata$Npredictse)),
        col=scales::alpha("red",.5))
points(subset(ch2exp, Herbivory==0)$Isopods, subset(ch2exp, Herbivory==0)$Plant, pch=19, col=subset(ch2exp, Herbivory==0)$Label, cex=1.5)

plot(Plant~Isopods, data=subset(ch2exp, Herbivory==1), pch=19, col=Label,
     xlab="Isopods (#)", ylab="", type="n",main="Herbivory",
     ylim=c(0,14), xlim=c(0,20))
text(x=0, y=13.5, label="b", cex=2)
subofdata = subset(newdata, Label==0 & Herbivory==1)
lines(subofdata$Isopods,
      subofdata$Npredict, cex=2, lty=2)
polygon(c(subofdata$Isopods,rev(subofdata$Isopods)),
        c(subofdata$Npredict+subofdata$Npredictse,rev(subofdata$Npredict-subofdata$Npredictse)),
        col=scales::alpha("grey",.5))

subofdata = subset(newdata, Label==1 & Herbivory==1)
lines(subofdata$Isopods,
      subofdata$Npredict, cex=2, lty=2, col="red")
polygon(c(subofdata$Isopods,rev(subofdata$Isopods)),
        c(subofdata$Npredict+subofdata$Npredictse,rev(subofdata$Npredict-subofdata$Npredictse)),
        col=scales::alpha("red",.5))
points(subset(ch2exp, Herbivory==1)$Isopods, subset(ch2exp, Herbivory==1)$Plant, pch=10, col=subset(ch2exp, Herbivory==1)$Label, cex=1.5)

Isopodrange = seq(0,12, by=1)

Wormmean = rep(mean(ch2exp$Worms), times=length(Isopodrange))
SIRmean = rep(mean(ch2exp$SIR_June, na.rm=T), times=length(Isopodrange))
N_Sol_June1 = rep(1, times=length(Isopodrange))
Label1 = rep(1, times=length(Isopodrange))

#Herbivory 0
newdatplant1 = as.data.frame(cbind(Isopodrange, Wormmean, SIRmean, rep(0, length(Isopodrange)), N_Sol_June1))
names(newdatplant1) = c("Isopods", "Worms", "SIR_June", "Herbivory", "N_Sol_June")
newdatplant1["Herbivory"] = as.factor(newdatplant1$Herbivory)

newplant1 = predict(P6, newdata=newdatplant1, se.fit=T)$fit

newdatplant2 = newdatplant1
newdatplant2["N_Sol_June"] = newdatplant1$N_Sol_June + 1
newplant2 = predict(P6, newdata=newdatplant2, se.fit=T)$fit

newdatplant3 = newdatplant1
newdatplant3["N_Sol_June"] = newdatplant1$N_Sol_June + 2
newplant3 = predict(P6, newdata=newdatplant3, se.fit=T)$fit

newdatplant4 = newdatplant1
newdatplant4["N_Sol_June"] = newdatplant1$N_Sol_June + 3
newplant4 = predict(P6, newdata=newdatplant4, se.fit=T)$fit

newdatplant5 = newdatplant1
newdatplant5["N_Sol_June"] = 0
newplant5 = predict(P6, newdata=newdatplant5, se.fit=T)$fit

plot(newplant1~Isopodrange, pch=19, type="b", ylim=c(0,14), col="brown", 
     ylab="", xlab="", xlim=c(0,20))
text(x=0, y=13.5, label="c", cex=2)
points(newplant5~Isopodrange, pch=19, type="b", col="black")
points(newplant2~Isopodrange, pch=19, type="b", col="orange")
points(newplant3~Isopodrange, pch=19, type="b", col="red")
points(newplant4~Isopodrange, pch=19, type="b", col="green")

legend("bottomright", legend=c(0,1,2,3,4),
       col=c("black", "brown", "orange", "red", "green", "darkgreen"),
       pch=19, title=expression(italic(Solidago)), bty="n",
       x.intersp = 0.5, y.intersp = 1)


# Herbivory 1
newdatplant1 = as.data.frame(cbind(Isopodrange, Wormmean, SIRmean, Label1, N_Sol_June1))
names(newdatplant1) = c("Isopods", "Worms", "SIR_June", "Herbivory", "N_Sol_June")
newdatplant1["Herbivory"] = as.factor(newdatplant1$Herbivory)

newplant1 = predict(P6, newdata=newdatplant1, se.fit=T)$fit

newdatplant2 = newdatplant1
newdatplant2["N_Sol_June"] = newdatplant1$N_Sol_June + 1
newplant2 = predict(P6, newdata=newdatplant2, se.fit=T)$fit

newdatplant3 = newdatplant1
newdatplant3["N_Sol_June"] = newdatplant1$N_Sol_June + 2
newplant3 = predict(P6, newdata=newdatplant3, se.fit=T)$fit

newdatplant4 = newdatplant1
newdatplant4["N_Sol_June"] = newdatplant1$N_Sol_June + 3
newplant4 = predict(P6, newdata=newdatplant4, se.fit=T)$fit

newdatplant5 = newdatplant1
newdatplant5["N_Sol_June"] = 0
newplant5 = predict(P6, newdata=newdatplant5, se.fit=T)$fit

plot(newplant1~Isopodrange, pch=19, type="b", ylim=c(0,14), col="brown",
     ylab="", xlab="", xlim=c(0,20))
text(x=0, y=13.5, label="d", cex=2)
points(newplant5~Isopodrange, pch=19, type="b", col="black")
points(newplant2~Isopodrange, pch=19, type="b", col="orange")
points(newplant3~Isopodrange, pch=19, type="b", col="red")
points(newplant4~Isopodrange, pch=19, type="b", col="green")
mtext(text="Aboveground Plant Biomass (g)",side=2,line=-1,outer=TRUE)
mtext(text="Isopods (#)",side=1,line=0,outer=TRUE)
dev.off()

