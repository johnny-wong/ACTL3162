# Fitting a G.L.M. (frequency) model to the privately owned, comprehensive insurance claims data in 1975.
#
# read the data file

PCarIns.tmp <- read.csv("PrivateCarIns1975-Data.csv")
attach(PCarIns.tmp)
PCarIns <- PCarIns.tmp[PCarIns.tmp$Numb.Claims>0,]
attach(PCarIns)

# understand and summarize the data.
print(summary(PCarIns))
# convert to categorical
Pol.Age <- factor(Pol.Age)
Car.Group <- factor(Car.Group)
Veh.Age <- factor(Veh.Age)

# hist(Avg.Claims,breaks=50,prob=T)
# alternative glm models
#pcarins.glm <- glm(Numb.Claims~Pol.Age+Car.Group+Veh.Age, family=poisson)
pcarins.glm <- glm(Numb.Claims~Pol.Age+Car.Group+Veh.Age+Pol.Age*Car.Group+Pol.Age*Veh.Age+Car.Group*Veh.Age, family=poisson)
print(summary(pcarins.glm))
print(anova(pcarins.glm, test="Chi"))

# plot(Numb.Claims,predict(pcarins.glm))
 
par(mfrow=c(2,2))
plot(Numb.Claims,fitted(pcarins.glm),xlab="observed amounts",ylab="fitted values",main="Observed vs Predicted",pch=20)
abline(0,1)
plot(fitted(pcarins.glm),resid(pcarins.glm,type="deviance"),xlab="fitted values",ylab="deviance residuals",main="Fitted vs Residuals",pch=20)
abline(0,0)
qqnorm(resid(pcarins.glm,type="pearson"), xlab="quantiles of Std Normal",ylab="Pearson residuals",pch=20)
qqline(resid(pcarins.glm))

# diagnostic plots
library(boot)
pcarins.glm.diag <- glm.diag(pcarins.glm)
# print(names(pcarins.glm.diag))
# glm.diag.plots(pcarins.glm,pcarins.glm.diag)
# plot(pcarins.glm)

ts.plot(pcarins.glm.diag$cook,type="h",gpars=list(xlab="observation",ylab="Cook statistic",main="Influential Observations",pch=20))

# plot predicted values with each regressor variable
# plot(Pol.Age,resid(pcarins.glm))
# plot(Car.Group,resid(pcarins.glm))
# plot(Veh.Age,resid(pcarins.glm))
