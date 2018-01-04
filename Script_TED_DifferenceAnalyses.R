# Elif Kardas - 2017
##########################################################
#
#
#        Change in TED analysis for all associations 
#
#
##########################################################
# I. STATISTICAL TEST - QUALITATIVE DESCRIPTION
##########################################################


library(readr)
require(lattice)
Assoc <- read_delim("C:/Users/Elif Ka/OneDrive/Mémoire/Excel tables/Associations' TED table - OUI.csv", 
                                         ";", escape_double = FALSE, col_types = cols(Association = col_factor(levels = c("1", 
                                                                                                                          "2", "3")), Concentration = col_factor(levels = c("0", 
                                                                                                                                                                            "1", "10", "100")), TED = col_number(), 
                                                                                      Time = col_factor(levels = c("1","2", 
                                                                                                                   "3","4"))), trim_ws = TRUE)

attach(Assoc)


xyplot(log(TED)~Time|Concentration, groups=Association, data=Assoc[-135,],
       layout=c(2,2),
       index.cond=list(c(3,4,1,2)),
       auto.key =list(space="top", columns=3, 
                               title="Association", cex.title=1), 
       par.settings = list(superpose.symbol = list(pch = c(1,2,3), 
                                                   cex = 0.6, 
                                                   col = c("#FF0066", 
                                                           "#66FF00","#0066CC"))),
       strip = strip.custom(style = 1, bg="darkolivegreen3", strip.names = c(TRUE, TRUE)))

Assoc$Time<- as.numeric(Assoc$Time) #to set "Time" as a numeric variable
Assoc$Concentration<-as.factor(Assoc$Concentration) # to set the other as factors
Assoc$Association <-as.factor(Assoc$Association)


##########################################################
# II. STATISTICAL TEST - QUANTITATIVE DESCRIPTION
##########################################################
# II.1.          APPLICATION CONDITIONS
##########################################################
# 1. tests of the application conditions for global TED (All associations)
#The linear model is: 
lm_assoc_NN <- lm(log(TED)~Time*Concentration*Association, data=Assoc)
### Residuals = normally distributed ?
### for: lm_assocw <- lm(log(TED)~Time*Concentration*Association, data=Assoc)
residus_NN <- residuals(lm_assoc_NN) # NN for non normally distributed
shapiro.test(residus_NN) # NN for non normally distributed - shapiro.test gives a p-value < 0.05
hist(residus_NN) # not a normal distribution
plot(lm_assoc_NN) # line 135 seems to be out of the normal Q-Q -> we decided to take it off, see explaination in master thesis
### for: lm_assoc <- lm(log(TED)~Time*Concentration*Association, data=Assoc[-135,])
lm_assoc <- lm(log(TED)~Time*Concentration*Association, data=Assoc[-135,]) # with a deleted line
residus <- residuals(lm_assoc)
shapiro.test(residus) # Residuals are normally distirbuted as p-value > 0.05
hist(residus) # seems to be a better distribution
plot(lm_assoc) 
##########################################################
# II.2.            BEST MODEL
##########################################################
Assocb <- Assoc[-135,]
lm_assoc <- lm(log(TED)~Time*Concentration*Association, data=Assocb)
nullmodel<- lm(log(TED)~1, data=Assoc)
summary(lm_assoc)
step(nullmodel, scope=lm_assoc, direction="both")
library(MuMIn)
options(na.action = "na.fail")
dredge(lm_assoc)

##########################################################
#
#
#        Change in TED analysis for each species 
#
#
##########################################################
# Load of all species in association data
library(readr)
spall <- read_delim("C:/Users/Elif Ka/OneDrive/Mémoire/Excel tables/Species TED in associations without analysis.csv", 
                                                           ";", escape_double = FALSE, col_types = cols(Association = col_factor(levels = c("1", 
                                                                                                                                            "2", "3")), Concentration = col_factor(levels = c("0", 
                                                                                                                                                                                              "1", "10", "100")), Species = col_factor(levels = c("1", 
                                                                                                                                                                                                                                                  "2", "3")), TED = col_number(), Time = col_number()), 
                                                           trim_ws = TRUE)
##################################
#
#           Species 1
#
##################################
# I. STATISTICAL TEST - QUALITATIVE DESCRIPTION
##########################################################
sp1 <- subset(spall, spall$Species==1)

attach(sp1)
xyplot(log(TED)~Time|Concentration, groups=Association, data=sp1, 
       auto.key =list(space="top", columns=3, 
                      title="Association", cex.title=1), 
       par.settings = list(superpose.symbol = list(pch = c(1,2,3), 
                                                   cex = 0.6, 
                                                   col = c("#FF0066", 
                                                           "#66FF00","#0066CC"))))

sp1$Time<- as.numeric(sp1$Time) #to set "Time" as a numeric variable
sp1$Concentration<-as.factor(sp1$Concentration)
sp1$Association <-as.factor(sp1$Association)

##########################################################
# II. STATISTICAL TEST - QUANTITATIVE DESCRIPTION
##########################################################
# II.1.          APPLICATION CONDITIONS
##########################################################
#The linear model is: 
lm_sp1 <- lm(log(TED)~Time*Concentration*Association, data=sp1)
plot(lm_sp1)
# 1. tests of the application conditions for each species
### Residuals = normally distributed ?
### for: lm_assocw <- lm(log(TED)~Time*Concentration*Association, data=Assoc)
residus <- residuals(lm_sp1) # NN for non normally distributed
shapiro.test(residus) # Residuals are normally distirbuted as p-value > 0.05 = 0.5883
hist(residus) 
plot(residus)
##########################################################
# II.2.            BEST MODEL
##########################################################
# 2. Try to find the best model
lm_sp1 <- lm(log(TED)~Time*Concentration*Association, data=sp1)
nullmodel<- lm(log(TED)~1, data=sp1)
summary(lm_sp1)
step(nullmodel, scope=lm_assoc, direction="both")

library(MuMIn)
options(na.action = "na.fail")
dredge(lm_sp1)

##################################
#
#           Species 2
#
##################################
# I. STATISTICAL TEST - QUALITATIVE DESCRIPTION
##########################################################
sp2 <- subset(spall, spall$Species==2)

attach(sp2)
xyplot(log(TED)~Time|Concentration, groups=Association, data=sp2, 
       auto.key =list(space="top", columns=3, 
                      title="Association", cex.title=1), 
       par.settings = list(superpose.symbol = list(pch = c(1,2,3), 
                                                   cex = 0.6, 
                                                   col = c("#FF0066", 
                                                           "#66FF00","#0066CC"))))

sp2$Time<- as.numeric(sp2$Time) #to set "Time" as a numeric variable
sp2$Concentration<-as.factor(sp2$Concentration)
sp2$Association <-as.factor(sp2$Association)

##########################################################
# II. STATISTICAL TEST - QUANTITATIVE DESCRIPTION
##########################################################
# II.1.          APPLICATION CONDITIONS
##########################################################
#The linear model is: 
lm_sp2 <- lm(log(TED)~Time*Concentration*Association, data=sp2)
plot(lm_sp2)
##########################################################
# II.2.            BEST MODEL
##########################################################
# 2. Try to find the best model
lm_sp2 <- lm(log(TED)~Time*Concentration*Association, data=sp2)
nullmodel<- lm(log(TED)~1, data=sp2)
summary(lm_sp2)
step(nullmodel, scope=lm_assoc, direction="both")
library(MuMIn)
options(na.action = "na.fail")
dredge(lm_sp2)
##################################
#
#           Species 3
#
##################################
# I. STATISTICAL TEST - QUALITATIVE DESCRIPTION
##########################################################
sp3 <- subset(spall, spall$Species==3)

attach(sp3)
xyplot(log(TED+1)~Time|Concentration, groups=Association, data=sp3, 
       auto.key =list(space="top", columns=3, 
                      title="Association", cex.title=1), 
       par.settings = list(superpose.symbol = list(pch = c(1,2,3), 
                                                   cex = 0.6, 
                                                   col = c("#FF0066", 
                                                           "#66FF00","#0066CC"))))

sp3$Time<- as.numeric(sp3$Time) #to set "Time" as a numeric variable
sp3$Concentration<-as.factor(sp3$Concentration)
sp3$Association <-as.factor(sp3$Association)

##########################################################
# II. STATISTICAL TEST - QUANTITATIVE DESCRIPTION
##########################################################
# II.1.          APPLICATION CONDITIONS
##########################################################
#The linear model is: 
sp3 <- na.omit(sp3)
lm_sp3 <- lm(log(TED+1)~Time*Concentration*Association, data=sp3)
plot(lm_sp3)
##########################################################
# II.2.            BEST MODEL
##########################################################
# 2. Try to find the best model
lm_sp3 <- lm(log(TED+1)~Time*Concentration*Association, data=sp3)
nullmodel<- lm(log(TED+1)~1, data=sp3)
summary(lm_sp3)
step(nullmodel, scope=lm_assoc, direction="both")
library(MuMIn)
options(na.action = "na.fail")
dredge(lm_sp3)