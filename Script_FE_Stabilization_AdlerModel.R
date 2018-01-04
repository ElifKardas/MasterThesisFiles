# Elif Kardas - 2017
####################################################################################
#                                                                                  #
#           Stabilization and fitness differences of algae species                 #
#                  exposed to the mixture of pharmaceuticals                       #
#                 (calculation + analysis of the differences)                      #
#                                                                                  #
####################################################################################

####################################################################################
#                                                                                  #
#             A.  Stabilization and fitness differences calculations               #
#                                                                                  #
####################################################################################
# All the following values in the vectors have been calculated with Python. 
# See: https://github.com/ElifKardas/MasterThesisLambdaAlphaCalculationScripts
####################################################################################
#             A.1. For association 1- c0,c1,c10,c100                               #
####################################################################################

## Load the data 
l1 <- c(99.9999999999999,2.419107081218254152e+01, 2.283605497553037722e+01,
        7.664215307022573143e+00, 6.533259264632527596e+00, 6.058567750413423880e+00,
        1.027546032592450409e+01, 1.167290479085182753e+01, 7.475477678988427854e+00,
        5.023362034817069599e+01, 8.611707922666278492e+00, 1.806148581052184099e+01
        )
l2 <- c(30.8272190616564,3.624100420886670548e+01, 2.231884784827114387e+01,
        2.232230930575889971e+01, 1.455629061720881623e+01, 1.706202023697022341e+01,
        3.250366639058563578e+01,2.174485679638635816e+01, 2.689159579558246804e+01,
        3.754180978379002909e+01, 2.498481670767827367e+01, 3.468765972679111798e+01
        )

a12 <- c(3.602218692341460481e-26, 2.335607062802659301e-03, 2.348650314651990270e-02,
         1.105338991369805799e-17, 3.927714152158269826e-04, 1.757972109940056888e-20,
         4.184551532820311297e-21, 9.398612433213970538e-04, 8.995816353454879616e-24,
         2.782201661297010629e-04, 1.931376028667157463e-03, 2.512487601415229448e-21
         )

a21 <- c(3.595189503137681195e-03, 3.595189503137681195e-03, 2.981028461232772503e-04,
         2.162949244101708846e-03, 2.162949244101708846e-03, 2.162949244101708846e-03,
         3.479856063604090232e-03, 2.911633482363844023e-03, 3.479856063604090232e-03,
         4.260738188156296125e-03, 3.880050630807366750e-03, 4.260738188156296125e-03
         )
a11 <- c(3.466641376545980674e-02, 1.052147081439682406e-02, 8.013411464279527816e-03,
         2.333018871415890472e-03, 2.350385749518380495e-03, 2.433631553058584508e-03,
         3.884039453296155107e-03, 5.331274780456979855e-03, 3.085387685759056510e-03,
         2.427706107901490828e-02, 3.109218569415480009e-03, 7.894072880167198047e-03
         )
a22 <- c(3.891014119919793644e-03, 4.931400376467440401e-03, 2.608765225991797193e-03,
         2.899652854242486602e-03, 1.810471086099515450e-03, 1.965205596434788961e-03,
         4.083883605982628617e-03, 3.214096882981266282e-03, 3.275978963804103070e-03,
         4.811869561239899909e-03, 3.387959134860759128e-03, 4.772233101526628107e-03
         )

## Stabilization for species 1
denom1 <- 1+(a12/a22)*(l2-1)
stab_a1_sp1 <- l2/denom1
stab_a1_sp1
## Stabilization for species 2 
denom2 <- 1+(a21/a11)*(l1-1)
stab_a1_sp2<- l1/denom2
stab_a1_sp2

## Fitness differences for species 1
FD_a1_sp1 <- l1/l2
FD_a1_sp1
## Fitness differences for species 2
FD_a1_sp2<- l2/l1
FD_a1_sp2

####################################################################################
#             A.2. For association 2- c0,c1,c10,c100                               #
####################################################################################

## Load the data 
l1 <- c(99.9999999999999,2.419107081218254152e+01, 2.283605497553037722e+01,
        7.664215307022573143e+00, 6.533259264632527596e+00, 6.058567750413423880e+00,
        1.027546032592450409e+01, 1.167290479085182753e+01, 7.475477678988427854e+00,
        5.023362034817069599e+01, 8.611707922666278492e+00, 1.806148581052184099e+01
)
l2 <- l3 <- c(2.211149962822818438e+01, 2.221679366422421609e+01, 2.123783831627480367e+01,
              3.850557943450959808e+00, 3.950287575930283346e+00, 4.467166878639175387e+00,
              9.692570078353361041e+00, 9.076268394985822496e+00, 7.101851368792775787e+00,
              1.163952255137472669e+01, 6.698594803756167337e+00, 6.452089075894939185e+00
              )

a12 <- a13 <- c(1.976070223451066588e-02, 3.463023134156176588e-03,1.239197565682879270e-03,
                2.361134906329751676e-03, 1.417238089148222359e-04, 2.361134906329751676e-03,
                2.310723513430507410e-29, 2.342998206637909799e-04, 1.601076249162856871e-32,
                9.333135633316901741e-04, 1.883130969591412087e-30, 2.626469931060703127e-04
                )

a21 <- a31 <- c(5.485602824941562990e-03,5.485602824941227321e-03,5.485602824941562990e-03,
                2.577490034168191269e-03,2.577490034168193004e-03,2.577490034168193004e-03,
                8.547875477512198103e-03,8.547875477512199838e-03,8.547875477512199838e-03,
                8.104097994650588666e-03,8.104097994650588666e-03,8.104097994650588666e-03
                )

a11 <- c(3.466641376545980674e-02, 1.052147081439682406e-02, 8.013411464279527816e-03,
         2.333018871415890472e-03, 2.350385749518380495e-03, 2.433631553058584508e-03,
         3.884039453296155107e-03, 5.331274780456979855e-03, 3.085387685759056510e-03,
         2.427706107901490828e-02, 3.109218569415480009e-03, 7.894072880167198047e-03
)
a22 <- a33 <- c(5.978498416761827798e-03,5.793968353154222047e-03,4.832258537747503510e-03,
                2.297425169712937355e-03,2.522943223541755133e-03,2.841078134617889238e-03,
                1.194380294127509248e-02,9.091579557870665268e-03,6.495388344200125376e-03,
                1.320187490844237400e-02,6.809225310120870303e-03,6.268070063650764065e-03
                )
## Stabilization for species 1
denom1 <- 1+(a12/a22)*(l2-1)
stab_a2_sp1 <- l2/denom1
stab_a2_sp1
## Stabilization for species 3
denom2 <- 1+(a21/a11)*(l1-1)
stab_a2_sp3<- l1/denom2
stab_a2_sp3

## Fitness differences for species 1
FD_a2_sp1 <- l1/l2
FD_a2_sp1
## Fitness differences for species 3
FD_a2_sp3<- l2/l1
FD_a2_sp3

####################################################################################
#             A.3. For association 3- c0,c1,c10,c100                               #
####################################################################################

## Load the data 
l1 <- l2 <- c(30.8272190616564,3.624100420886670548e+01, 2.231884784827114387e+01,
        2.232230930575889971e+01, 1.455629061720881623e+01, 1.706202023697022341e+01,
        3.250366639058563578e+01,2.174485679638635816e+01, 2.689159579558246804e+01,
        3.754180978379002909e+01, 2.498481670767827367e+01, 3.468765972679111798e+01
        )
l2 <- l3 <- c(2.211149962822818438e+01, 2.221679366422421609e+01, 2.123783831627480367e+01,
              3.850557943450959808e+00, 3.950287575930283346e+00, 4.467166878639175387e+00,
              9.692570078353361041e+00, 9.076268394985822496e+00, 7.101851368792775787e+00,
              1.163952255137472669e+01, 6.698594803756167337e+00, 6.452089075894939185e+00
              )

a12 <- a23 <- c(1.963164784057737627e-32,1.412552516523240308e-30,2.402274013839225688e-29,
                3.003026666090067604e-31,3.425443189263727696e-33,2.956191549544889623e-32,
                8.611408981622866896e-05,1.312454981313247370e-30,1.516127063864731840e-35,
                8.129154737687635293e-32, 3.836294471026215503e-04,4.771392665877218039e-29
                )

a21 <- a32 <- c(5.485602824941562990e-03,5.485602824941562122e-03,5.485602824941560388e-03,
                2.577490033541577925e-03,2.577490034168193004e-03,2.577490034168189968e-03,
                8.547875477512199838e-03,8.547875474001622592e-03,8.547875477512199838e-03,
                8.104097994650548767e-03,8.104097994650588666e-03,8.104097992845815321e-03
                )

a11 <- a22 <- c(3.891014119919793644e-03, 4.931400376467440401e-03, 2.608765225991797193e-03,
            2.899652854242486602e-03, 1.810471086099515450e-03, 1.965205596434788961e-03,
            4.083883605982628617e-03, 3.214096882981266282e-03, 3.275978963804103070e-03,
            4.811869561239899909e-03, 3.387959134860759128e-03, 4.772233101526628107e-03
            )
a22 <- a33 <- c(5.978498416761827798e-03,5.793968353154222047e-03,4.832258537747503510e-03,
                2.297425169712937355e-03,2.522943223541755133e-03,2.841078134617889238e-03,
                1.194380294127509248e-02, 9.091579557870665268e-03,6.495388344200125376e-03,
                1.320187490844237400e-02,6.809225310120870303e-03,6.268070063650764065e-03
                )
## Stabilization for species 2
denom1 <- 1+(a12/a22)*(l2-1)
stab_a3_sp2 <- l2/denom1
stab_a3_sp2
## Stabilization for species 3
denom2 <- 1+(a21/a11)*(l1-1)
stab_a3_sp3<- l1/denom2
stab_a3_sp3

## Fitness differences for species 2
FD_a2_sp2 <- l1/l2
FD_a2_sp2
## Fitness differences for species 3
FD_a2_sp3<- l2/l1
FD_a2_sp3



####################################################################################
#                                                                                  #
#     B. Analysis of differences of stabilization and fitness differences values   #
#      of algae communities exposed to the mixture of pharmaceuticals              #
#                                                                                  #
####################################################################################

# Load the data 
library(readr)
Adler_FD_STAB <- read_delim("C:/Users/Elif Ka/Desktop/Adler_FE_STAB.csv", 
                            ";", escape_double = FALSE, col_types = cols(Association = col_factor(levels = c("1", 
                                                                                                             "2", "3")), Concentration = col_factor(levels = c("0", 
                                                                                                                                                               "1", "10", "100")), Species = col_factor(levels = c("1", 
                                                                                                                                                                                                                   "2", "3")), 
                                                                         Replicate = col_factor(levels = c("1", "2", "3")), Stabilization = col_number()), 
                            trim_ws = TRUE)

#####################################################
#                                                   #
#             B.1. Qualitative description          #
#                                                   #
#####################################################

library(lattice)
library(latticeExtra)

newfactor<-as.factor(paste(Adler_FD_STAB[1:72,]$Concentration, Adler_FD_STAB[1:72,]$Species))
# To see plot of association 1:
p1 <- xyplot(log(FD)~Stabilization|Association, data=Adler_FD_STAB[1:72,], group=newfactor, 
             ylab="log(Fitness equivalence)", xlab="Stabilization",
             index.cond=list(c(1,1)),
             layout=c(1,1),
             par.settings = list(superpose.symbol = list(pch = c(rep(1, 3),rep(2, 3), rep(3,3), rep(8,3)), 
                                                         cex = 1),
                                 bg=rep(c("brown2", "chartreuse2", "dodgerblue2"), 4)), 
             scales = list(x = list(at = c(0,10,20,30), 
                                    limits = c(-1,36)), 
                           alternating=1, 
                           tck=c(0.5,0)), 
             col=rep(c("brown2", "chartreuse2", "dodgerblue2"), 4),
             strip = strip.custom(style = 1, bg="darkolivegreen3", strip.names = c(TRUE, TRUE))
             )


p1
x= seq(0,40, by = 0.025)
y= exp(1/x)
p2 <- xyplot(y~x, type="l", col="black")
p2 
p1 + as.layer(p2)

# To see plot of association 2:
p1 <- xyplot(log(FD)~Stabilization|Association, data=Adler_FD_STAB[1:72,], group=newfactor, 
       ylab="log(Fitness equivalence)", xlab="Stabilization",
       index.cond=list(c(1,2)),
       layout=c(1,1),
       par.settings = list(superpose.symbol = list(pch = c(rep(1, 3),rep(2, 3), rep(3,3), rep(8,3)), 
                                                   cex = 1),
                           bg=rep(c("brown2", "chartreuse2", "dodgerblue2"), 4)), 
       scales = list(x = list(at = c(0,5,10), 
                              limits = c(-1,10)), 
                     alternating=1, 
                     tck=c(0.5,0)), 
       col=rep(c("brown2", "chartreuse2", "dodgerblue2"), 4),
       strip = strip.custom(style = 1, bg="darkolivegreen3", strip.names = c(TRUE, TRUE))
)

p1
x= seq(0,40, by = 0.025)
y= exp(1/x)
p2 <- xyplot(y~x, type="l", col="black")
p2 
p1 + as.layer(p2)

# To see plot of association 3:
p1 <- xyplot(log(FD)~Stabilization|Association, data=Adler_FD_STAB[1:72,], group=newfactor, 
             ylab="log(Fitness equivalence)", xlab="Stabilization",
             index.cond=list(c(1,3)),
             layout=c(1,1),
             par.settings = list(superpose.symbol = list(pch = c(rep(1, 3),rep(2, 3), rep(3,3), rep(8,3)), 
                                                         cex = 1),
                                 bg=rep(c("brown2", "chartreuse2", "dodgerblue2"), 4)), 
             scales = list(x = list( 
               at = c(0,10,20), 
               limits = c(-1,23)), 
               alternating=1, 
               tck=c(0.5,0)), 
             col=rep(c("brown2", "chartreuse2", "dodgerblue2"), 4),
             strip = strip.custom(style = 1, bg="darkolivegreen3", strip.names = c(TRUE, TRUE))
             )

p1
x= seq(0,40, by = 0.025)
y= exp(1/x)
p2 <- xyplot(y~x, type="l", col="black")
p2 
p1 + as.layer(p2)

# To manually make the legend (don't do it like this if you know lattice very well :-))
plot(x=c(1,2,3,4, 0.5, 2, 3.5), y=c(1,1,1,1,-1, -1, -1), pch=c(1,2,3,8,15 , 15,15), xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0,5), col=c(rep("black",4),"brown2", "chartreuse2", "dodgerblue2"), ylim=c(-1.5, 1.5))
text(x=c(1,2,3,4, 0.5, 2, 3.5), y=c(1,1,1,1,-1, -1, -1), labels=c("0", "1", "10", "100", "Species 1", "Species 2", "Species 3"), pos=4)
text(x=c(3,2.7), y=c(1.3,-0.7), labels = c("Concentration", "Species"), pos=2)

#####################################################
#                                                   #
#             B.2. Quantitative description         #
#                                                   #
#####################################################

####################################################
#   B.2.1. ANOVAs for association 1 and species 1  #
####################################################

Adler_FD_STAB_ass1_sp1 <- Adler_FD_STAB[1:12,]

## Stabilization 
ND_ANOVA <- aov(log(Stabilization)~Concentration, data=Adler_FD_STAB_ass1_sp1)
### test of normality
res <- residuals(ND_ANOVA)
shapiro.test(res) # p-value >0.05, NRH0, Normally distributed
### test of homoscedasticity
bartlett.test(log(Stabilization)~Concentration, data=Adler_FD_STAB_ass1_sp1) # p-value > 0.05, NRH0 > Variances are homogenous
### test itself
summary(ND_ANOVA)
anova(ND_ANOVA)
plot(Stabilization~Concentration, data=Adler_FD_STAB_ass1_sp1)

## FD
FD_ANOVA <- aov(log(FD)~Concentration, data=Adler_FD_STAB_ass1_sp1)
### test of normality
res <- residuals(FD_ANOVA)
shapiro.test(res) # p-value >0.05, NRH0, Normally distributed
### test of homoscedasticity
bartlett.test(log(FD)~Concentration, data=Adler_FD_STAB_ass1_sp1) # p-value > 0.05, NRH0 > Variances are homogenous
### test itself
summary(FD_ANOVA)
anova(FD_ANOVA)
plot(FD~Concentration, data=Adler_FD_STAB_ass1_sp1, ylab="Fitness equivalences")


####################################################
#   B.2.2. ANOVAs for association 1 and species 2  #
####################################################

Adler_FD_STAB_ass1_sp2 <- Adler_FD_STAB[13:24,]

## Stabilization 
ND_ANOVA <- aov(log(Stabilization)~Concentration, data=Adler_FD_STAB_ass1_sp2)
### test of normality
res <- residuals(ND_ANOVA)
shapiro.test(res) # p_value > 0.05 normally distributed 
### test of homoscedasticity
bartlett.test(log(Stabilization)~Concentration, data=Adler_FD_STAB_ass1_sp2) # variances are homogenous at the 0.1% threshold.
### test itself
summary(ND_ANOVA)
anova(ND_ANOVA)
plot(Stabilization~Concentration, data=Adler_FD_STAB_ass1_sp2)
title(main = "Strength of stabilization of association 1, depending on pharmaceutical concentrations")

## FD
FD_ANOVA <- aov(FD~Concentration, data=Adler_FD_STAB_ass1_sp2)
### test of normality
res <- residuals(FD_ANOVA)
shapiro.test(res) # p_value > 0.05 normally distributed
### test of homoscedasticity
bartlett.test(FD~Concentration, data=Adler_FD_STAB_ass1_sp2) # pvalue>0.05, variances are homogenous
### test itself
summary(FD_ANOVA)
anova(FD_ANOVA)
plot(FD~Concentration, data=Adler_FD_STAB_ass1_sp2, ylab="Fitness equivalences")

####################################################
#   B.2.3. ANOVAs for association 2 and species 1  #
####################################################

Adler_FD_STAB_ass2_sp1 <- Adler_FD_STAB[25:36,]

## Stabilization
ND_ANOVA <- aov(Stabilization~Concentration, data=Adler_FD_STAB_ass2_sp1)
### test of normality
res <- residuals(ND_ANOVA)
shapiro.test(res) # p-value >0.05, normally distributed
### test of homoscedasticity
bartlett.test(Stabilization~Concentration, data=Adler_FD_STAB_ass2_sp1) #variances are homogenous, p-value>0.05
### test itself
summary(ND_ANOVA)
anova(ND_ANOVA)
plot(Stabilization~Concentration, data=Adler_FD_STAB_ass2_sp1)
title(main = "Strength of stabilization of association 2, depending on pharmaceutical concentrations")

## FD
FD_ANOVA <- aov(log(FD)~Concentration, data=Adler_FD_STAB_ass2_sp1)
### test of normality
res <- residuals(FD_ANOVA)
shapiro.test(res) # normally distributed p-value >0.05
### test of homoscedasticity
bartlett.test(log(FD)~Concentration, data=Adler_FD_STAB_ass2_sp1) # variances are homogenous p-value>0.05
### test itself
summary(FD_ANOVA)
anova(FD_ANOVA)
plot(FD~Concentration, data=Adler_FD_STAB_ass2_sp1, ylab="Fitness equivalences")
title(main="Fitness equivalences of association 2, depending on pharmaceutical concentrations", ylab="Fitness equivalences")

####################################################
#   B.2.4. ANOVAs for association 2 and species 3  #
####################################################

Adler_FD_STAB_ass2_sp3 <- Adler_FD_STAB[37:48,]

## Stabilization
ND_ANOVA <- aov(log(Stabilization)~Concentration, data=Adler_FD_STAB_ass2_sp3)
### test of normality
res <- residuals(ND_ANOVA)
shapiro.test(res) # normally distributed p-value>0.05
### test of homoscedasticity
bartlett.test(log(Stabilization)~Concentration, data=Adler_FD_STAB_ass2_sp3) # variances are homogenous at the 0.1% threshold
### test itself
summary(ND_ANOVA)
anova(ND_ANOVA)
plot(Stabilization~Concentration, data=Adler_FD_STAB_ass2_sp3)

## FD
FD_ANOVA <- aov(FD~Concentration, data=Adler_FD_STAB_ass2_sp3)
### test of normality
res <- residuals(FD_ANOVA)
shapiro.test(res) # normally distributed p-value >0.05
### test of homoscedasticity
bartlett.test(FD~Concentration, data=Adler_FD_STAB_ass2_sp3) # variances are homogenous p-value >0.05
### test itself
summary(FD_ANOVA)
anova(FD_ANOVA)
plot(FD~Concentration, data=Adler_FD_STAB_ass2_sp3, ylab="Fitness equivalences")


####################################################
#   B.2.5. ANOVAs for association 3 and species 2  #
####################################################

Adler_FD_STAB_ass3_sp2 <- Adler_FD_STAB[49:60,]

## For Stab
ND_ANOVA <- aov(Stabilization~Concentration, data=Adler_FD_STAB_ass3_sp2)
### test of normality
res <- residuals(ND_ANOVA)
shapiro.test(res) # normally distributed p-value>0.05
### test of homoscedasticity
bartlett.test(Stabilization~Concentration, data=Adler_FD_STAB_ass3_sp2) # variances are homogenous at the 1% threshold 
### test itself
summary(ND_ANOVA)
anova(ND_ANOVA)
plot(Stabilization~Concentration, data=Adler_FD_STAB_ass3_sp2)
title(main = "Strength of stabilization of association 3, depending on pharmaceutical concentrations")

## For FD
FD_ANOVA <- aov(FD~Concentration, data=Adler_FD_STAB_ass3_sp2)
### test of normality
res <- residuals(FD_ANOVA)
shapiro.test(res) # normally distributed p-value>0.05
### test of homoscedasticity
bartlett.test(FD~Concentration, data=Adler_FD_STAB_ass3_sp2) # variances are homogenous p-value>0.05
### test itself
summary(FD_ANOVA)
anova(FD_ANOVA)
plot(FD~Concentration, data=Adler_FD_STAB_ass3_sp2, ylab="Fitness equivalences")
title(main="Fitness equivalences of association 3, depending on pharmaceutical concentrations", ylab="Fitness equivalences")

####################################################
#   B.2.6. ANOVAs for association 3 and species 3  #
####################################################

Adler_FD_STAB_ass3_sp3 <- Adler_FD_STAB[61:72,]

## Stabilization 
ND_ANOVA <- aov(Stabilization~Concentration, data=Adler_FD_STAB_ass3_sp3)
### test of normality
res<- residuals(ND_ANOVA)
shapiro.test(res) # normally distributed p-value>0.05
### test of homoscedasticity
bartlett.test(Stabilization~Concentration, data=Adler_FD_STAB_ass3_sp3) # variances are homogenous p-value >0.05
### test itself
summary(ND_ANOVA)
anova(ND_ANOVA)
plot(Stabilization~Concentration, data=Adler_FD_STAB_ass3_sp3)

## FD
FD_ANOVA <- aov(FD~Concentration, data=Adler_FD_STAB_ass3_sp3)
### test of normality
res <- residuals(FD_ANOVA)
shapiro.test(res) # normally distributed p-value>0.05
### test of homoscedasticity
bartlett.test(FD~Concentration, data=Adler_FD_STAB_ass3_sp3) # variances are homogenous p-value>0.05
### test itself
summary(FD_ANOVA)
anova(FD_ANOVA)
plot(FD~Concentration, data=Adler_FD_STAB_ass3_sp3, ylab="Fitness equivalences")
