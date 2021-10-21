###############################################################################
############################EXACT BETA-POISSON#################################
###############################################################################

#Clean the R environment
rm(list=ls())


#Require necessary packages
if("rstudioapi" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("rstudioapi"); require(rstudioapi)} else{require(rstudioapi)}
if("ggplot2" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("ggplot2"); require(ggplot2)} else{require(ggplot2)}
if("scales" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("scales"); require(scales)} else{require(scales)}

#Establish the Working Directory
this.dir<-dirname(rstudioapi::getSourceEditorContext()$path)   
setwd(this.dir)


#Source the dataset script 
source('dataset.R')

#Source the script developed by Xie et al.
source('bp_code.R')



#############
#TOTAL CYSTS#
#############

#Set seed to make the simulations reproducible
set.seed(1)

#Proglottids

##Assign the arguments
attach(db_exp_total)
y = db_exp_total[exp_route=="proglottids",]$nInf
Ns = db_exp_total[exp_route=="proglottids",]$nExp
Ds = db_exp_total[exp_route=="proglottids",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_total)

##Prior parameters estimation by maximum likelihood
total.prog.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = total.prog.param[[1]], b = total.prog.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.total.prog <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                               UL = out$p975c, exp_route = "proglottids")

total.prog.sum <- list(alpha.mean = mean(out$alpha),
                       alpha.median = median(out$alpha),
                       alpha.p25c = quantile(out$alpha, probs = 0.025),
                       alpha.p975c = quantile(out$alpha, probs = 0.975),
                       beta.mean = mean(out$beta), 
                       beta.median = median(out$beta),
                       beta.p25c = quantile(out$beta, probs = 0.025),
                       beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.total.prog <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.total.prog, 'parameters_total_prog.csv')

#Eggs

##Assign the arguments
attach(db_exp_total)
y = db_exp_total[exp_route=="eggs",]$nInf
Ns = db_exp_total[exp_route=="eggs",]$nExp
Ds = db_exp_total[exp_route=="eggs",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_total)

##Prior parameters estimation by maximum likelihood
total.egg.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = total.egg.param[[1]], b = total.egg.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.total.egg <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                               UL = out$p975c, exp_route = "eggs")

total.egg.sum <- list(alpha.mean = mean(out$alpha),
                       alpha.median = median(out$alpha),
                       alpha.p25c = quantile(out$alpha, probs = 0.025),
                       alpha.p975c = quantile(out$alpha, probs = 0.975),
                       beta.mean = mean(out$beta), 
                       beta.median = median(out$beta),
                       beta.p25c = quantile(out$beta, probs = 0.025),
                       beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.total.egg <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.total.egg, 'parameters_total_egg.csv')

#Beetles

##Assign the arguments
attach(db_exp_total)
y = db_exp_total[exp_route=="beetles",]$nInf
Ns = db_exp_total[exp_route=="beetles",]$nExp
Ds = db_exp_total[exp_route=="beetles",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_total)

##Prior parameters estimation by maximum likelihood
total.bee.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = total.bee.param[[1]], b = total.bee.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.total.bee <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                              UL = out$p975c, exp_route = "beetles")

total.bee.sum <- list(alpha.mean = mean(out$alpha),
                      alpha.median = median(out$alpha),
                      alpha.p25c = quantile(out$alpha, probs = 0.025),
                      alpha.p975c = quantile(out$alpha, probs = 0.975),
                      beta.mean = mean(out$beta), 
                      beta.median = median(out$beta),
                      beta.p25c = quantile(out$beta, probs = 0.025),
                      beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.total.bee <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.total.bee, 'parameters_total_bee.csv')

#Carotid

##Assign the arguments
attach(db_oral_total)
y = db_oral_total[exp_route_b=="carotid",]$nInf
Ns = db_oral_total[exp_route_b=="carotid",]$nExp
Ds = db_oral_total[exp_route_b=="carotid",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_oral_total)

##Prior parameters estimation by maximum likelihood
total.car.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = total.car.param[[1]], b = total.car.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.total.car <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                              UL = out$p975c, exp_route = "carotid")

total.car.sum <- list(alpha.mean = mean(out$alpha),
                      alpha.median = median(out$alpha),
                      alpha.p25c = quantile(out$alpha, probs = 0.025),
                      alpha.p975c = quantile(out$alpha, probs = 0.975),
                      beta.mean = mean(out$beta), 
                      beta.median = median(out$beta),
                      beta.p25c = quantile(out$beta, probs = 0.025),
                      beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.total.car <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.total.car, 'parameters_total_car.csv')

#Oral

##Assign the arguments
attach(db_oral_total)
y = db_oral_total[exp_route_b=="oral",]$nInf
Ns = db_oral_total[exp_route_b=="oral",]$nExp
Ds = db_oral_total[exp_route_b=="oral",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_oral_total)

##Prior parameters estimation by maximum likelihood
total.oral.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = total.oral.param[[1]], b = total.oral.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.total.oral <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                               UL = out$p975c, exp_route = "oral")

total.oral.sum <- list(alpha.mean = mean(out$alpha),
                       alpha.median = median(out$alpha),
                       alpha.p25c = quantile(out$alpha, probs = 0.025),
                       alpha.p975c = quantile(out$alpha, probs = 0.975),
                       beta.mean = mean(out$beta), 
                       beta.median = median(out$beta),
                       beta.p25c = quantile(out$beta, probs = 0.025),
                       beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.total.oral <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.total.oral, 'parameters_total_oral.csv')

##PLOT DESIGN

###Merge the predicted values from each route of exposure
curve.total.1 <- rbind(curve.total.prog, curve.total.egg, curve.total.bee, 
                       curve.total.car)
curve.total.1$exp_route <- factor(curve.total.1$exp_route,
                                levels = c("proglottids", "eggs", "beetles", "carotid"))

###Storage the mean infective dose (ID50) of each route of exposure
ID50.total.1 <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                           ID50 = c(0.83, 7.52, 1.08, 0.88))
ID50.total.1$exp_route <- factor(ID50.total.1$exp_route, 
                                 levels = c("proglottids", "eggs", "beetles", "carotid"))

###PLOT 1
exp_route.labs <- c("PROGLOTTIDS", "EGGS", "BEETLES","CAROTID")
names(exp_route.labs) <- c("proglottids","eggs","beetles","carotid")

total_labels <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                           label = c("A","B","C","D"))
total_labels$exp_route <- factor(total_labels$exp_route,
                                 levels = c("proglottids","eggs","beetles","carotid"))

total.plot.1 <- ggplot(curve.total.1, aes(x = dose, y = p)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +
  geom_vline(data = ID50.total.1, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("P(inf)") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_total, aes(y = nInf/nExp), size = 3) +
  geom_errorbarh(data=geom_total, 
                 aes(y = nInf/nExp, xmax = dose + dose_sd, 
                     xmin = dose - dose_sd, height = .05)) +
  facet_wrap(~exp_route,
             labeller = labeller(exp_route = exp_route.labs),
             nrow = 4) +
  theme(strip.text.x = element_text(size = 18, face = "bold")) +
  geom_text(x = -2, y = 1, aes(label = label), data = total_labels)

total.plot.1



#################
#VESICULAR CYSTS#
#################

#Set seed to make the simulations reproducible
set.seed(1)

#Proglottids

##Assign the arguments
attach(db_exp_ves)
y = db_exp_ves[exp_route=="proglottids",]$nInf
Ns = db_exp_ves[exp_route=="proglottids",]$nExp
Ds = db_exp_ves[exp_route=="proglottids",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_ves)

##Prior parameters estimation by maximum likelihood
ves.prog.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = ves.prog.param[[1]], b = ves.prog.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.ves.prog <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                               UL = out$p975c, exp_route = "proglottids")

ves.prog.sum <- list(alpha.mean = mean(out$alpha),
                       alpha.median = median(out$alpha),
                       alpha.p25c = quantile(out$alpha, probs = 0.025),
                       alpha.p975c = quantile(out$alpha, probs = 0.975),
                       beta.mean = mean(out$beta), 
                       beta.median = median(out$beta),
                       beta.p25c = quantile(out$beta, probs = 0.025),
                       beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.ves.prog <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.ves.prog, 'parameters_ves_prog.csv')

#Eggs

##Assign the arguments
attach(db_exp_ves)
y = db_exp_ves[exp_route=="eggs",]$nInf
Ns = db_exp_ves[exp_route=="eggs",]$nExp
Ds = db_exp_ves[exp_route=="eggs",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_ves)

##Prior parameters estimation by maximum likelihood
ves.egg.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = ves.egg.param[[1]], b = ves.egg.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.ves.egg <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                              UL = out$p975c, exp_route = "eggs")

ves.egg.sum <- list(alpha.mean = mean(out$alpha),
                      alpha.median = median(out$alpha),
                      alpha.p25c = quantile(out$alpha, probs = 0.025),
                      alpha.p975c = quantile(out$alpha, probs = 0.975),
                      beta.mean = mean(out$beta), 
                      beta.median = median(out$beta),
                      beta.p25c = quantile(out$beta, probs = 0.025),
                      beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.ves.egg <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.ves.egg, 'parameters_ves_egg.csv')

#Beetles

##Assign the arguments
attach(db_exp_ves)
y = db_exp_ves[exp_route=="beetles",]$nInf
Ns = db_exp_ves[exp_route=="beetles",]$nExp
Ds = db_exp_ves[exp_route=="beetles",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_ves)

##Prior parameters estimation by maximum likelihood
ves.bee.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = ves.bee.param[[1]], b = ves.bee.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.ves.bee <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                              UL = out$p975c, exp_route = "beetles")

ves.bee.sum <- list(alpha.mean = mean(out$alpha),
                      alpha.median = median(out$alpha),
                      alpha.p25c = quantile(out$alpha, probs = 0.025),
                      alpha.p975c = quantile(out$alpha, probs = 0.975),
                      beta.mean = mean(out$beta), 
                      beta.median = median(out$beta),
                      beta.p25c = quantile(out$beta, probs = 0.025),
                      beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.ves.bee <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.ves.bee, 'parameters_ves_bee.csv')

#Carotid

##Assign the arguments
attach(db_oral_ves)
y = db_oral_ves[exp_route_b=="carotid",]$nInf
Ns = db_oral_ves[exp_route_b=="carotid",]$nExp
Ds = db_oral_ves[exp_route_b=="carotid",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_oral_ves)

##Prior parameters estimation by maximum likelihood
ves.car.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = ves.car.param[[1]], b = ves.car.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.ves.car <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                            UL = out$p975c, exp_route = "carotid")

ves.car.sum <- list(alpha.mean = mean(out$alpha),
                    alpha.median = median(out$alpha),
                    alpha.p25c = quantile(out$alpha, probs = 0.025),
                    alpha.p975c = quantile(out$alpha, probs = 0.975),
                    beta.mean = mean(out$beta), 
                    beta.median = median(out$beta),
                    beta.p25c = quantile(out$beta, probs = 0.025),
                    beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.ves.car <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.ves.car, 'parameters_ves_car.csv')

#Oral

##Assign the arguments
attach(db_oral_ves)
y = db_oral_ves[exp_route_b=="oral",]$nInf
Ns = db_oral_ves[exp_route_b=="oral",]$nExp
Ds = db_oral_ves[exp_route_b=="oral",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_oral_ves)

##Prior parameters estimation by maximum likelihood
ves.oral.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = ves.oral.param[[1]], b = ves.oral.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.ves.oral <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                             UL = out$p975c, exp_route = "oral")

ves.oral.sum <- list(alpha.mean = mean(out$alpha),
                     alpha.median = median(out$alpha),
                     alpha.p25c = quantile(out$alpha, probs = 0.025),
                     alpha.p975c = quantile(out$alpha, probs = 0.975),
                     beta.mean = mean(out$beta), 
                     beta.median = median(out$beta),
                     beta.p25c = quantile(out$beta, probs = 0.025),
                     beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.ves.oral <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.ves.oral, 'parameters_ves_oral.csv')

##PLOT DESIGN

###Merge the predicted values from each route of exposure
curve.ves.1 <- rbind(curve.ves.prog, curve.ves.egg, curve.ves.bee, 
                       curve.ves.car)
curve.ves.1$exp_route <- factor(curve.ves.1$exp_route,
                                  levels = c("proglottids", "eggs", "beetles", "carotid"))

###Storage the mean infective dose (ID50) of each route of exposure
ID50.ves.1 <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                           ID50 = c(1744.43, 829.83, 7.54, 0.88))
ID50.ves.1$exp_route <- factor(ID50.ves.1$exp_route, 
                                 levels = c("proglottids", "eggs", "beetles", "carotid"))

###PLOT 1
exp_route.labs <- c("PROGLOTTIDS", "EGGS", "BEETLES","CAROTID")
names(exp_route.labs) <- c("proglottids","eggs","beetles","carotid")

ves_labels <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                         label = c("A","B","C","D"))
ves_labels$exp_route <- factor(ves_labels$exp_route, 
                               levels = c("proglottids","eggs","beetles","carotid"))

ves.plot.1 <- ggplot(curve.ves.1, aes(x = dose, y = p)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +
  geom_vline(data = ID50.ves.1, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("P(inf)") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_ves, aes(y = nInf/nExp), size = 3) +
  geom_errorbarh(data=geom_ves, 
                 aes(y = nInf/nExp, xmax = dose + dose_sd, 
                     xmin = dose - dose_sd, height = .05)) +
  facet_wrap(~exp_route,
             labeller = labeller(exp_route = exp_route.labs),
             nrow = 4) +
  theme(strip.text.x = element_text(size = 18, face = "bold")) +
  geom_text(x = -2, y = 1, aes(label = label), data = ves_labels)

ves.plot.1



#############
#BRAIN CYSTS#
#############

#Set seed to make the simulations reproducible
set.seed(1)

#Proglottids

##Assign the arguments
attach(db_exp_brain)
y = db_exp_brain[exp_route=="proglottids",]$nInf
Ns = db_exp_brain[exp_route=="proglottids",]$nExp
Ds = db_exp_brain[exp_route=="proglottids",]$dose_b
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_brain)

##Prior parameters estimation by maximum likelihood
brain.prog.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = brain.prog.param[[1]], b = brain.prog.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.brain.prog <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                             UL = out$p975c, exp_route = "proglottids")

brain.prog.sum <- list(alpha.mean = mean(out$alpha),
                     alpha.median = median(out$alpha),
                     alpha.p25c = quantile(out$alpha, probs = 0.025),
                     alpha.p975c = quantile(out$alpha, probs = 0.975),
                     beta.mean = mean(out$beta), 
                     beta.median = median(out$beta),
                     beta.p25c = quantile(out$beta, probs = 0.025),
                     beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.brain.prog <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.brain.prog, 'parameters_brain_prog.csv')

#Eggs

##Assign the arguments
attach(db_exp_brain)
y = db_exp_brain[exp_route=="eggs",]$nInf
Ns = db_exp_brain[exp_route=="eggs",]$nExp
Ds = db_exp_brain[exp_route=="eggs",]$dose_b
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_brain)

##Prior parameters estimation by maximum likelihood
brain.egg.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = brain.egg.param[[1]], b = brain.egg.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.brain.egg <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                            UL = out$p975c, exp_route = "eggs")

brain.egg.sum <- list(alpha.mean = mean(out$alpha),
                    alpha.median = median(out$alpha),
                    alpha.p25c = quantile(out$alpha, probs = 0.025),
                    alpha.p975c = quantile(out$alpha, probs = 0.975),
                    beta.mean = mean(out$beta), 
                    beta.median = median(out$beta),
                    beta.p25c = quantile(out$beta, probs = 0.025),
                    beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.brain.egg <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.brain.egg, 'parameters_brain_egg.csv')

#Beetles

##Assign the arguments
attach(db_exp_brain)
y = db_exp_brain[exp_route=="beetles",]$nInf
Ns = db_exp_brain[exp_route=="beetles",]$nExp
Ds = db_exp_brain[exp_route=="beetles",]$dose_b
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_brain)

##Prior parameters estimation by maximum likelihood
brain.bee.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = brain.bee.param[[1]], b = brain.bee.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.brain.bee <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                            UL = out$p975c, exp_route = "beetles")

brain.bee.sum <- list(alpha.mean = mean(out$alpha),
                    alpha.median = median(out$alpha),
                    alpha.p25c = quantile(out$alpha, probs = 0.025),
                    alpha.p975c = quantile(out$alpha, probs = 0.975),
                    beta.mean = mean(out$beta), 
                    beta.median = median(out$beta),
                    beta.p25c = quantile(out$beta, probs = 0.025),
                    beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.brain.bee <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.brain.bee, 'parameters_brain_bee.csv')

#Carotid

##Assign the arguments
attach(db_oral_brain)
y = db_oral_brain[exp_route_b=="carotid",]$nInf
Ns = db_oral_brain[exp_route_b=="carotid",]$nExp
Ds = db_oral_brain[exp_route_b=="carotid",]$dose_b
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_oral_brain)

##Prior parameters estimation by maximum likelihood
brain.car.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = brain.car.param[[1]], b = brain.car.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.brain.car <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                              UL = out$p975c, exp_route = "carotid")

brain.car.sum <- list(alpha.mean = mean(out$alpha),
                      alpha.median = median(out$alpha),
                      alpha.p25c = quantile(out$alpha, probs = 0.025),
                      alpha.p975c = quantile(out$alpha, probs = 0.975),
                      beta.mean = mean(out$beta), 
                      beta.median = median(out$beta),
                      beta.p25c = quantile(out$beta, probs = 0.025),
                      beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.brain.car <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.brain.car, 'parameters_brain_car.csv')

#Oral

##Assign the arguments
attach(db_oral_brain)
y = db_oral_brain[exp_route_b=="oral",]$nInf
Ns = db_oral_brain[exp_route_b=="oral",]$nExp
Ds = db_oral_brain[exp_route_b=="oral",]$dose_b
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_oral_brain)

##Prior parameters estimation by maximum likelihood
brain.oral.param <- newBPoptim(Ds, Ns, y, approxi = FALSE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = brain.oral.param[[1]], b = brain.oral.param[[2]],
                 simuN=5000, Dsim,appxi = FALSE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list and in a csv file
curve.brain.oral <- data.frame(dose = Dsim, p = out$p50c, LL = out$p25c,
                               UL = out$p975c, exp_route = "oral")

brain.oral.sum <- list(alpha.mean = mean(out$alpha),
                       alpha.median = median(out$alpha),
                       alpha.p25c = quantile(out$alpha, probs = 0.025),
                       alpha.p975c = quantile(out$alpha, probs = 0.975),
                       beta.mean = mean(out$beta), 
                       beta.median = median(out$beta),
                       beta.p25c = quantile(out$beta, probs = 0.025),
                       beta.p975c = quantile(out$beta, probs = 0.975))

alphabeta.brain.oral <- data.frame(alpha = out$alpha, beta = out$beta)

write.csv(alphabeta.brain.oral, 'parameters_brain_oral.csv')


##PLOT DESIGN

###Merge the predicted values from each route of exposure
curve.brain.1 <- rbind(curve.brain.prog, curve.brain.egg, curve.brain.bee, 
                     curve.brain.car)
curve.brain.1$exp_route <- factor(curve.brain.1$exp_route,
                                levels = c("proglottids", "eggs", "beetles", "carotid"))

###Storage the mean infective dose (ID50) of each route of exposure
ID50.brain.1 <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                         ID50 = c(8.12e04, 5.72e05, 8.05e40, 34.88))
ID50.brain.1$exp_route <- factor(ID50.brain.1$exp_route, 
                               levels = c("proglottids", "eggs", "beetles", "carotid"))

###PLOT 1
exp_route.labs <- c("PROGLOTTIDS", "EGGS", "BEETLES","CAROTID")
names(exp_route.labs) <- c("proglottids","eggs","beetles","carotid")

brain_labels <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                           label = c("A","B","C","D"))
brain_labels$exp_route <- factor(brain_labels$exp_route, 
                                 levels = c("proglottids","eggs","beetles","carotid"))

brain.plot.1 <- ggplot(curve.brain.1, aes(x = dose, y = p)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +
  geom_vline(data = ID50.brain.1, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("P(inf)") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_brain, aes(y = nInf/nExp), size = 3) +
  geom_errorbarh(data=geom_brain, 
                 aes(y = nInf/nExp, xmax = dose + dose_sd, 
                     xmin = dose - dose_sd, height = .05)) +
  facet_wrap(~exp_route,
             labeller = labeller(exp_route = exp_route.labs),
             nrow = 4) +
  theme(strip.text.x = element_text(size = 18, face = "bold")) +
  geom_text(x = -2, y = 1, aes(label = label), data = brain_labels)

brain.plot.1



#ORAL PLOT
curve.total.oral$cyst_type <- "total"
curve.ves.oral$cyst_type <- "vesicular"
curve.brain.oral$cyst_type <- "brain"

curve.oral.1 <- rbind(curve.total.oral, curve.ves.oral, curve.brain.oral)

curve.oral.1$cyst_type <- factor(curve.oral.1$cyst_type, 
                                 levels = c("total", "vesicular","brain"))


ID50.oral.1 <- data.frame(cyst_type = c("total","vesicular","brain"),
                           ID50 = c(2.48, 59.31, 3.72e05))
ID50.oral.1$cyst_type <- factor(ID50.oral.1$cyst_type, 
                                levels = c("total","vesicular","brain"))

###PLOT
cyst_type.labs <- c("ANY CYST", "VIABLE CYSTS", "BRAIN CYSTS")
names(cyst_type.labs) <- c("total", "vesicular", "brain")

oral_labels <- data.frame(cyst_type = c("total","vesicular","brain"), 
                          label = c("A","B","C"))
oral_labels$cyst_type <- factor(oral_labels$cyst_type,
                                levels = c("total","vesicular","brain"))

oral.plot <- ggplot(curve.oral.1, aes(x = dose, y = p)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +
  geom_vline(data = ID50.oral.1, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("Probability of infection") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_oral, aes(y = nInf/nExp, shape = exp_route), size = 2) +
  geom_errorbarh(data=geom_oral, 
                 aes(y = nInf/nExp, xmax = dose + dose_sd, 
                     xmin = dose - dose_sd, height = .05)) +
  facet_wrap(~cyst_type,
             labeller = labeller(cyst_type = cyst_type.labs),
             nrow = 3) +
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  scale_shape_discrete(name = "Route of exposure",
                       labels = c("Proglottids", "Eggs", "Beetles")) +
  geom_text(x = -2, y = 1, aes(label = label), data = oral_labels)

oral.plot


#EGGS PLOT
curve.total.egg$cyst_type <- "total"
curve.ves.egg$cyst_type <- "vesicular"
curve.brain.egg$cyst_type <- "brain"

curve.egg <- rbind(curve.total.egg, curve.ves.egg, curve.brain.egg)

curve.egg$cyst_type <- factor(curve.egg$cyst_type, 
                              levels = c("total", "vesicular", "brain"))


ID50.egg <- data.frame(cyst_type = c("total","vesicular","brain"),
                          ID50 = c(8.49, 739, 6.4e04))
ID50.egg$cyst_type <- factor(ID50.egg$cyst_type,
                             levels = c("total","vesicular","brain"))

###PLOT
cyst_type.labs <- c("ANY CYST", "VIABLE CYSTS", "BRAIN CYSTS")
names(cyst_type.labs) <- c("total", "vesicular", "brain")

egg_labels <- data.frame(cyst_type = c("total","vesicular","brain"), 
                         label = c("A","B","C"))
egg_labels$cyst_type <- factor(egg_labels$cyst_type,
                               levels = c("total","vesicular","brain"))

egg.plot <- ggplot(curve.egg, aes(x = dose, y = p)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +
  geom_vline(data = ID50.egg, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("Probability of infection") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_egg, aes(y = nInf/nExp), size = 2) +
  facet_wrap(~cyst_type,
             labeller = labeller(cyst_type = cyst_type.labs),
             nrow = 3) +
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
  geom_text(x = -2, y = 1, aes(label = label), data = egg_labels)

egg.plot


#PROGLOTTIDS PLOT
curve.total.prog$cyst_type <- "total"
curve.ves.prog$cyst_type <- "vesicular"
curve.brain.prog$cyst_type <- "brain"

curve.prog <- rbind(curve.total.prog, curve.ves.prog, curve.brain.prog)

curve.prog$cyst_type <- factor(curve.prog$cyst_type, 
                              levels = c("total", "vesicular", "brain"))


ID50.prog <- data.frame(cyst_type = c("total","vesicular","brain"),
                       ID50 = c(0.83, 1707, 8.49e04))
ID50.prog$cyst_type <- factor(ID50.prog$cyst_type,
                             levels = c("total","vesicular","brain"))

###PLOT
cyst_type.labs <- c("ANY CYST", "VIABLE CYSTS", "BRAIN CYSTS")
names(cyst_type.labs) <- c("total", "vesicular", "brain")

prog_labels <- data.frame(cyst_type = c("total","vesicular","brain"), 
                         label = c("A","B","C"))
prog_labels$cyst_type <- factor(prog_labels$cyst_type,
                               levels = c("total","vesicular","brain"))

prog.plot <- ggplot(curve.prog, aes(x = dose, y = p)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +
  geom_vline(data = ID50.prog, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("Probability of infection") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_prog, aes(y = nInf/nExp), size = 2) +
  geom_errorbarh(data=geom_prog, 
                 aes(y = nInf/nExp, xmax = dose + dose_sd, 
                     xmin = dose - dose_sd, height = .05)) +
  facet_wrap(~cyst_type,
             labeller = labeller(cyst_type = cyst_type.labs),
             nrow = 3) +
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
  geom_text(x = -2, y = 1, aes(label = label), data = prog_labels)

prog.plot


#BEETLES PLOT
curve.total.bee$cyst_type <- "total"
curve.ves.bee$cyst_type <- "vesicular"
curve.brain.bee$cyst_type <- "brain"

curve.bee <- rbind(curve.total.bee, curve.ves.bee, curve.brain.bee)

curve.bee$cyst_type <- factor(curve.bee$cyst_type, 
                              levels = c("total", "vesicular", "brain"))


ID50.bee <- data.frame(cyst_type = c("total","vesicular","brain"),
                       ID50 = c(1.09, 8.49, 8.11e40))
ID50.bee$cyst_type <- factor(ID50.bee$cyst_type,
                             levels = c("total","vesicular","brain"))

###PLOT
cyst_type.labs <- c("ANY CYST", "VIABLE CYSTS", "BRAIN CYSTS")
names(cyst_type.labs) <- c("total", "vesicular", "brain")

bee_labels <- data.frame(cyst_type = c("total","vesicular","brain"), 
                         label = c("A","B","C"))
bee_labels$cyst_type <- factor(bee_labels$cyst_type,
                               levels = c("total","vesicular","brain"))

bee.plot <- ggplot(curve.bee, aes(x = dose, y = p)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +
  geom_vline(data = ID50.bee, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("Probability of infection") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_bee, aes(y = nInf/nExp), size = 2) +
  geom_errorbarh(data=geom_bee, 
                 aes(y = nInf/nExp, xmax = dose + dose_sd, 
                     xmin = dose - dose_sd, height = .05)) +
  facet_wrap(~cyst_type,
             labeller = labeller(cyst_type = cyst_type.labs),
             nrow = 3) +
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
  geom_text(x = -2, y = 1, aes(label = label), data = bee_labels)

bee.plot


#CAROTID PLOT
curve.total.car$cyst_type <- "total"
curve.ves.car$cyst_type <- "vesicular"
curve.brain.car$cyst_type <- "brain"

curve.car <- rbind(curve.total.car, curve.ves.car, curve.brain.car)

curve.car$cyst_type <- factor(curve.car$cyst_type, 
                              levels = c("total", "vesicular", "brain"))


ID50.car <- data.frame(cyst_type = c("total","vesicular","brain"),
                       ID50 = c(0.86, 0.91, 32.7))
ID50.car$cyst_type <- factor(ID50.car$cyst_type,
                             levels = c("total","vesicular","brain"))

###PLOT
cyst_type.labs <- c("ANY CYST", "VIABLE CYSTS", "BRAIN CYSTS")
names(cyst_type.labs) <- c("total", "vesicular", "brain")

car_labels <- data.frame(cyst_type = c("total","vesicular","brain"), 
                         label = c("A","B","C"))
car_labels$cyst_type <- factor(car_labels$cyst_type,
                               levels = c("total","vesicular","brain"))

car.plot <- ggplot(curve.car, aes(x = dose, y = p)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +
  geom_vline(data = ID50.car, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("Probability of infection") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_car, aes(y = nInf/nExp), size = 2) +
  facet_wrap(~cyst_type,
             labeller = labeller(cyst_type = cyst_type.labs),
             nrow = 3) +
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
  geom_text(x = -2, y = 1, aes(label = label), data = car_labels)

car.plot
