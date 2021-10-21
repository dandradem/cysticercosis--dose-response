###############################################################################
#########################APPROXIMATE BETA-POISSON##############################
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
total.prog.param <- newBPoptim(Ds, Ns, y, approxi = TRUE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = total.prog.param[[1]], b = total.prog.param[[2]],
                 simuN=5000, Dsim,appxi = TRUE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list
curve.total.prog <- data.frame(dose_a = Dsim, p = out$p50c, LL = out$p25c,
                               UL = out$p975c, exp_route = "proglottids")

total.prog.sum <- list(alpha.mean = mean(out$alpha),
                       alpha.median = median(out$alpha),
                       alpha.p25c = quantile(out$alpha, probs = 0.025),
                       alpha.p975c = quantile(out$alpha, probs = 0.975),
                       beta.mean = mean(out$beta), 
                       beta.median = median(out$beta),
                       beta.p25c = quantile(out$beta, probs = 0.025),
                       beta.p975c = quantile(out$beta, probs = 0.975))

#Eggs

##Assign the arguments
attach(db_exp_total)
y = db_exp_total[exp_route=="eggs",]$nInf
Ns = db_exp_total[exp_route=="eggs",]$nExp
Ds = db_exp_total[exp_route=="eggs",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_total)

##Prior parameters estimation by maximum likelihood
total.egg.param <- newBPoptim(Ds, Ns, y, approxi = TRUE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = total.egg.param[[1]], b = total.egg.param[[2]],
                 simuN=5000, Dsim,appxi = TRUE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list
curve.total.egg <- data.frame(dose_a = Dsim, p = out$p50c, LL = out$p25c,
                               UL = out$p975c, exp_route = "eggs")

total.egg.sum <- list(alpha.mean = mean(out$alpha),
                       alpha.median = median(out$alpha),
                       alpha.p25c = quantile(out$alpha, probs = 0.025),
                       alpha.p975c = quantile(out$alpha, probs = 0.975),
                       beta.mean = mean(out$beta), 
                       beta.median = median(out$beta),
                       beta.p25c = quantile(out$beta, probs = 0.025),
                       beta.p975c = quantile(out$beta, probs = 0.975))

#Beetles

##Assign the arguments
attach(db_exp_total)
y = db_exp_total[exp_route=="beetles",]$nInf
Ns = db_exp_total[exp_route=="beetles",]$nExp
Ds = db_exp_total[exp_route=="beetles",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_total)

##Prior parameters estimation by maximum likelihood
total.bee.param <- newBPoptim(Ds, Ns, y, approxi = TRUE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = total.bee.param[[1]], b = total.bee.param[[2]],
                 simuN=5000, Dsim,appxi = TRUE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list
curve.total.bee <- data.frame(dose_a = Dsim, p = out$p50c, LL = out$p25c,
                              UL = out$p975c, exp_route = "beetles")

total.bee.sum <- list(alpha.mean = mean(out$alpha),
                      alpha.median = median(out$alpha),
                      alpha.p25c = quantile(out$alpha, probs = 0.025),
                      alpha.p975c = quantile(out$alpha, probs = 0.975),
                      beta.mean = mean(out$beta), 
                      beta.median = median(out$beta),
                      beta.p25c = quantile(out$beta, probs = 0.025),
                      beta.p975c = quantile(out$beta, probs = 0.975))

#Carotid

##Assign the arguments
attach(db_exp_total)
y = db_exp_total[exp_route=="carotid",]$nInf
Ns = db_exp_total[exp_route=="carotid",]$nExp
Ds = db_exp_total[exp_route=="carotid",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_total)

##Prior parameters estimation by maximum likelihood
total.car.param <- newBPoptim(Ds, Ns, y, approxi = TRUE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = total.car.param[[1]], b = total.car.param[[2]],
                 simuN=5000, Dsim,appxi = TRUE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list
curve.total.car <- data.frame(dose_a = Dsim, p = out$p50c, LL = out$p25c,
                              UL = out$p975c, exp_route = "carotid")

total.car.sum <- list(alpha.mean = mean(out$alpha),
                      alpha.median = median(out$alpha),
                      alpha.p25c = quantile(out$alpha, probs = 0.025),
                      alpha.p975c = quantile(out$alpha, probs = 0.975),
                      beta.mean = mean(out$beta), 
                      beta.median = median(out$beta),
                      beta.p25c = quantile(out$beta, probs = 0.025),
                      beta.p975c = quantile(out$beta, probs = 0.975))

#PLOT DESIGN

#Merge the predicted values from each route of exposure
curve.total <- rbind(curve.total.prog,curve.total.car,curve.total.egg,
                     curve.total.bee)
curve.total$exp_route <- factor(curve.total$exp_route,
                                levels = c("proglottids", "eggs", "beetles", "carotid"))

#Storage the mean infective dose (ID50) of each route of exposure
ID50.total <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                         ID50 = c(NA, 5.85, NA, NA))
ID50.total$exp_route <- factor(ID50.total$exp_route, 
                               levels = c("proglottids", "eggs", "beetles", "carotid"))


##PLOT
exp_route.labs <- c("PROGLOTTIDS", "EGGS", "BEETLES","CAROTID")
names(exp_route.labs) <- c("proglottids", "eggs", "beetles","carotid")

letter_labels <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                            label = c("A","B","C","D"))
letter_labels$exp_route <- factor(letter_labels$exp_route, 
                                  levels = c("proglottids","eggs","beetles","carotid"))

total.plot <- ggplot(curve.total, aes(x = dose_a, y = p)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +
  geom_vline(data = ID50.total, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("Probability of infection") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_total, aes(x = dose, y = nInf/nExp), size = 2) +
  geom_errorbarh(data=geom_total, 
                 aes(x = dose, y = nInf/nExp, xmax = dose + dose_sd, 
                     xmin = dose - dose_sd, height = .05)) +
  facet_wrap(~exp_route,
             labeller = labeller(exp_route = exp_route.labs),
             nrow = 4) +
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
  geom_text(x = -2, y = 1, aes(label = label), data = letter_labels)

total.plot




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
ves.prog.param <- newBPoptim(Ds, Ns, y, approxi = TRUE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = ves.prog.param[[1]], b = ves.prog.param[[2]],
                 simuN=5000, Dsim,appxi = TRUE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list
curve.ves.prog <- data.frame(dose_a = Dsim, p = out$p50c, LL = out$p25c,
                               UL = out$p975c, exp_route = "proglottids")

ves.prog.sum <- list(alpha.mean = mean(out$alpha),
                       alpha.median = median(out$alpha),
                       alpha.p25c = quantile(out$alpha, probs = 0.025),
                       alpha.p975c = quantile(out$alpha, probs = 0.975),
                       beta.mean = mean(out$beta), 
                       beta.median = median(out$beta),
                       beta.p25c = quantile(out$beta, probs = 0.025),
                       beta.p975c = quantile(out$beta, probs = 0.975))

#Eggs

##Assign the arguments
attach(db_exp_ves)
y = db_exp_ves[exp_route=="eggs",]$nInf
Ns = db_exp_ves[exp_route=="eggs",]$nExp
Ds = db_exp_ves[exp_route=="eggs",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_ves)

##Prior parameters estimation by maximum likelihood
ves.egg.param <- newBPoptim(Ds, Ns, y, approxi = TRUE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = ves.egg.param[[1]], b = ves.egg.param[[2]],
                 simuN=5000, Dsim,appxi = TRUE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list
curve.ves.egg <- data.frame(dose_a = Dsim, p = out$p50c, LL = out$p25c,
                              UL = out$p975c, exp_route = "eggs")

ves.egg.sum <- list(alpha.mean = mean(out$alpha),
                      alpha.median = median(out$alpha),
                      alpha.p25c = quantile(out$alpha, probs = 0.025),
                      alpha.p975c = quantile(out$alpha, probs = 0.975),
                      beta.mean = mean(out$beta), 
                      beta.median = median(out$beta),
                      beta.p25c = quantile(out$beta, probs = 0.025),
                      beta.p975c = quantile(out$beta, probs = 0.975))

#Beetles

##Assign the arguments
attach(db_exp_ves)
y = db_exp_ves[exp_route=="beetles",]$nInf
Ns = db_exp_ves[exp_route=="beetles",]$nExp
Ds = db_exp_ves[exp_route=="beetles",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_ves)

##Prior parameters estimation by maximum likelihood
ves.bee.param <- newBPoptim(Ds, Ns, y, approxi = TRUE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = ves.bee.param[[1]], b = ves.bee.param[[2]],
                 simuN=5000, Dsim,appxi = TRUE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list
curve.ves.bee <- data.frame(dose_a = Dsim, p = out$p50c, LL = out$p25c,
                              UL = out$p975c, exp_route = "beetles")

ves.bee.sum <- list(alpha.mean = mean(out$alpha),
                      alpha.median = median(out$alpha),
                      alpha.p25c = quantile(out$alpha, probs = 0.025),
                      alpha.p975c = quantile(out$alpha, probs = 0.975),
                      beta.mean = mean(out$beta), 
                      beta.median = median(out$beta),
                      beta.p25c = quantile(out$beta, probs = 0.025),
                      beta.p975c = quantile(out$beta, probs = 0.975))

#Carotid

##Assign the arguments
attach(db_exp_ves)
y = db_exp_ves[exp_route=="carotid",]$nInf
Ns = db_exp_ves[exp_route=="carotid",]$nExp
Ds = db_exp_ves[exp_route=="carotid",]$dose_a
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_ves)

##Prior parameters estimation by maximum likelihood
ves.car.param <- newBPoptim(Ds, Ns, y, approxi = TRUE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = ves.car.param[[1]], b = ves.car.param[[2]],
                 simuN=5000, Dsim,appxi = TRUE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list
curve.ves.car <- data.frame(dose_a = Dsim, p = out$p50c, LL = out$p25c,
                              UL = out$p975c, exp_route = "carotid")

ves.car.sum <- list(alpha.mean = mean(out$alpha),
                      alpha.median = median(out$alpha),
                      alpha.p25c = quantile(out$alpha, probs = 0.025),
                      alpha.p975c = quantile(out$alpha, probs = 0.975),
                      beta.mean = mean(out$beta), 
                      beta.median = median(out$beta),
                      beta.p25c = quantile(out$beta, probs = 0.025),
                      beta.p975c = quantile(out$beta, probs = 0.975))

#PLOT DESIGN

#Merge the predicted values from each route of exposure
curve.ves <- rbind(curve.ves.prog,curve.ves.car,curve.ves.egg,
                     curve.ves.bee)
curve.ves$exp_route <- factor(curve.ves$exp_route,
                                levels = c("proglottids", "eggs", "beetles", "carotid"))

#Storage the mean infective dose (ID50) of each route of exposure
ID50.ves <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                         ID50 = c(1200, 739.07, 1.04, NA))
ID50.ves$exp_route <- factor(ID50.ves$exp_route, 
                               levels = c("proglottids", "eggs", "beetles", "carotid"))


##PLOT
exp_route.labs <- c("PROGLOTTIDS", "EGGS", "BEETLES","CAROTID")
names(exp_route.labs) <- c("proglottids", "eggs", "beetles","carotid")

ves.plot <- ggplot(curve.ves, aes(x = dose_a, y = p)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +
  geom_vline(data = ID50.ves, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("Probability of infection") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_ves, aes(x = dose, y = nInf/nExp), size = 2) +
  geom_errorbarh(data=geom_ves, 
                 aes(x = dose, y = nInf/nExp, xmax = dose + dose_sd, 
                     xmin = dose - dose_sd, height = .05)) +
  facet_wrap(~exp_route,
             labeller = labeller(exp_route = exp_route.labs),
             nrow = 4) +
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
  geom_text(x = -2, y = 1, aes(label = label), data = letter_labels)

ves.plot




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
brain.prog.param <- newBPoptim(Ds, Ns, y, approxi = TRUE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = brain.prog.param[[1]], b = brain.prog.param[[2]],
                 simuN=5000, Dsim,appxi = TRUE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list
curve.brain.prog <- data.frame(dose_b = Dsim, p = out$p50c, LL = out$p25c,
                               UL = out$p975c, exp_route = "proglottids")

brain.prog.sum <- list(alpha.mean = mean(out$alpha),
                       alpha.median = median(out$alpha),
                       alpha.p25c = quantile(out$alpha, probs = 0.025),
                       alpha.p975c = quantile(out$alpha, probs = 0.975),
                       beta.mean = mean(out$beta), 
                       beta.median = median(out$beta),
                       beta.p25c = quantile(out$beta, probs = 0.025),
                       beta.p975c = quantile(out$beta, probs = 0.975))

#Eggs

##Assign the arguments
attach(db_exp_brain)
y = db_exp_brain[exp_route=="eggs",]$nInf
Ns = db_exp_brain[exp_route=="eggs",]$nExp
Ds = db_exp_brain[exp_route=="eggs",]$dose_b
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_brain)

##Prior parameters estimation by maximum likelihood
brain.egg.param <- newBPoptim(Ds, Ns, y, approxi = TRUE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = brain.egg.param[[1]], b = brain.egg.param[[2]],
                 simuN=5000, Dsim,appxi = TRUE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list
curve.brain.egg <- data.frame(dose_b = Dsim, p = out$p50c, LL = out$p25c,
                              UL = out$p975c, exp_route = "eggs")

brain.egg.sum <- list(alpha.mean = mean(out$alpha),
                      alpha.median = median(out$alpha),
                      alpha.p25c = quantile(out$alpha, probs = 0.025),
                      alpha.p975c = quantile(out$alpha, probs = 0.975),
                      beta.mean = mean(out$beta), 
                      beta.median = median(out$beta),
                      beta.p25c = quantile(out$beta, probs = 0.025),
                      beta.p975c = quantile(out$beta, probs = 0.975))

#Beetles

##Assign the arguments
attach(db_exp_brain)
y = db_exp_brain[exp_route=="beetles",]$nInf
Ns = db_exp_brain[exp_route=="beetles",]$nExp
Ds = db_exp_brain[exp_route=="beetles",]$dose_b
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_brain)

##Prior parameters estimation by maximum likelihood
brain.bee.param <- newBPoptim(Ds, Ns, y, approxi = TRUE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = brain.bee.param[[1]], b = brain.bee.param[[2]],
                 simuN=5000, Dsim,appxi = TRUE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list
curve.brain.bee <- data.frame(dose_b = Dsim, p = out$p50c, LL = out$p25c,
                              UL = out$p975c, exp_route = "beetles")

brain.bee.sum <- list(alpha.mean = mean(out$alpha),
                      alpha.median = median(out$alpha),
                      alpha.p25c = quantile(out$alpha, probs = 0.025),
                      alpha.p975c = quantile(out$alpha, probs = 0.975),
                      beta.mean = mean(out$beta), 
                      beta.median = median(out$beta),
                      beta.p25c = quantile(out$beta, probs = 0.025),
                      beta.p975c = quantile(out$beta, probs = 0.975))

#Carotid

##Assign the arguments
attach(db_exp_brain)
y = db_exp_brain[exp_route=="carotid",]$nInf
Ns = db_exp_brain[exp_route=="carotid",]$nExp
Ds = db_exp_brain[exp_route=="carotid",]$dose_b
Dsim <- exp(seq(log(1.00e200), log(1.00e-02), length=10000))
detach(db_exp_brain)

##Prior parameters estimation by maximum likelihood
brain.car.param <- newBPoptim(Ds, Ns, y, approxi = TRUE)

##Perform the parameters simulation
out = ApprxiBPci(Ns, Ds, a = brain.car.param[[1]], b = brain.car.param[[2]],
                 simuN=5000, Dsim,appxi = TRUE, pic=FALSE )

##Storage the predicted values in a data.frame and the parameters summary in
##a list
curve.brain.car <- data.frame(dose_b = Dsim, p = out$p50c, LL = out$p25c,
                              UL = out$p975c, exp_route = "carotid")

brain.car.sum <- list(alpha.mean = mean(out$alpha),
                      alpha.median = median(out$alpha),
                      alpha.p25c = quantile(out$alpha, probs = 0.025),
                      alpha.p975c = quantile(out$alpha, probs = 0.975),
                      beta.mean = mean(out$beta), 
                      beta.median = median(out$beta),
                      beta.p25c = quantile(out$beta, probs = 0.025),
                      beta.p975c = quantile(out$beta, probs = 0.975))

#PLOT DESIGN

#Merge the predicted values from each route of exposure
curve.brain <- rbind(curve.brain.prog,curve.brain.car,curve.brain.egg,
                     curve.brain.bee)
curve.brain$exp_route <- factor(curve.brain$exp_route,
                                levels = c("proglottids", "eggs", "beetles", "carotid"))

#Storage the mean infective dose (ID50) of each route of exposure
ID50.brain <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                         ID50 = c(89000, 64280.73, 166.5, 0.2984))
ID50.brain$exp_route <- factor(ID50.brain$exp_route, 
                               levels = c("proglottids", "eggs", "beetles", "carotid"))


##PLOT
exp_route.labs <- c("PROGLOTTIDS", "EGGS", "BEETLES","CAROTID")
names(exp_route.labs) <- c("proglottids", "eggs", "beetles","carotid")

brain.plot <- ggplot(curve.brain, aes(x = dose_b, y = p)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +
  geom_vline(data = ID50.brain, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("Probability of infection") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_brain, aes(x = dose, y = nInf/nExp), size = 2) +
  geom_errorbarh(data=geom_brain, 
                 aes(x = dose, y = nInf/nExp, xmax = dose + dose_sd, 
                     xmin = dose - dose_sd, height = .05)) +
  facet_wrap(~exp_route,
             labeller = labeller(exp_route = exp_route.labs),
             nrow = 4) +
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
  geom_text(x = -2, y = 1, aes(label = label), data = letter_labels)

brain.plot
