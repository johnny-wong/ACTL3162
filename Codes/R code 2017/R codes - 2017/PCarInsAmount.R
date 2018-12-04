# Fitting a G.L.M. (severity) model to the privately owned, comprehensive insurance claims data in 1975.
#
# read the data file
PCarIns.tmp <- read.csv("PrivateCarIns1975-Data.csv")
attach(PCarIns.tmp)
PCarIns <- PCarIns.tmp[PCarIns.tmp$Numb.Claims>0,] # remove the 3 categories with no claims
detach(PCarIns.tmp)
attach(PCarIns)

# understand and summarize the data.
print(summary(PCarIns))
# convert to categorical / treat Bonus as continuous
Pol.Age <- factor(Pol.Age)
Car.Group <- factor(Car.Group)
Veh.Age <- factor(Veh.Age)

#prepare for plots
par(mfrow=c(3,3))

# alternative glm models
#pcarins.glm <- glm(Avg.Claims~Pol.Age+Car.Group+Veh.Age+Pol.Age*Veh.Age+Pol.Age*Car.Group+Car.Group*Veh.Age, weights=Numb.Claims, family=Gamma)

#chosen model - comment out model above to run it
pcarins.glm <- glm(Avg.Claims~Pol.Age+Car.Group+Veh.Age, weights=Numb.Claims, family=Gamma)

#loot at glm summary with correlations
print(summary(pcarins.glm,corr=T))

#analysis of the deviance table
print(anova(pcarins.glm, test="Chi"))

#another diagnostic tool: try to drop each factor separately and see how it goes
print(drop1(pcarins.glm,text="Chisq"))

# analysis of residuals

 
plot(Avg.Claims,fitted(pcarins.glm),xlim=c(0,900),ylim=c(0,900),xlab="observed amounts",ylab="fitted values",main="Observed vs Predicted",pch=20)
abline(0,1)
plot(fitted(pcarins.glm),resid(pcarins.glm,type="deviance"),xlab="fitted values",ylab="deviance residuals",main="Fitted vs Residuals",pch=20)
abline(0,0)
qqnorm(resid(pcarins.glm,type="pearson"),xlim=c(-4,4),ylim=c(-4,4),xlab="quantiles of Std Normal",ylab="Pearson residuals",pch=20)
qqline(resid(pcarins.glm))

# diagnostic plots
library(boot)
# pcarins.glm.diag <- glm.diag(pcarins.glm)
# print(names(pcarins.glm.diag))
# glm.diag.plots(pcarins.glm,pcarins.glm.diag)
# plot(pcarins.glm)

ts.plot(pcarins.glm.diag$cook,type="h",gpars=list(xlab="observation",ylab="Cook statistic",main="Influential Observations",pch=20))

# plot predicted values with each regressor variable
# plot(Pol.Age,resid(pcarins.glm))
# plot(Car.Group,resid(pcarins.glm))
# plot(Veh.Age,resid(pcarins.glm))


# dataset for prediction (we predict those with missing claims)
PCarIns.pred <- PCarIns.tmp[PCarIns.tmp$Numb.Claims==0,]
Pol.Age.pred <- factor(PCarIns.pred$Pol.Age)
Cpol.Age.pred <- factor(PCarIns.pred$Cpol.Age)
Car.Group.pred <- factor(PCarIns.pred$Car.Group)
Veh.Age.pred <- factor(PCarIns.pred$Veh.Age)

predict.na <- predict.glm(pcarins.glm, data.frame(Pol.Age=Pol.Age.pred,Cpol.Age=Cpol.Age.pred,Car.Group=Car.Group.pred,Veh.Age=Veh.Age.pred),se.fit=TRUE,type="response")
print(predict.na)

par()

detach(PCarIns)

