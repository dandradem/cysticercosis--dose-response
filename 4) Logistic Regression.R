###############################################################################
#############################LOGISTIC REGRESSION###############################
###############################################################################

#Clean the R environment
rm(list=ls())


#Require necessary packages
if("rstudioapi" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("rstudioapi"); require(rstudioapi)} else{require(rstudioapi)}
if("car" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("car"); require(car)} else{require(car)}
if("nlstools" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("nlstools"); require(nlstools)} else{require(nlstools)}
if("ggplot2" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("ggplot2"); require(ggplot2)} else{require(ggplot2)}
if("scales" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("scales"); require(scales)} else{require(scales)}

#Establish the Working Directory
this.dir<-dirname(rstudioapi::getSourceEditorContext()$path)   
setwd(this.dir)

#Source the dataset script 
source('dataset.R')

#Create a function to calculate the mean infective dose (ID50)
lr.ID50 <- function(intercept,pardose) {
  ID50 <- (-intercept)/pardose
  return(ID50)
}



#############
#TOTAL CYSTS#
#############

#PROGLOTTIDS
total_prog <- glm(formula = total_cyst~dose_a, family = "binomial", 
                  data = db_bin[db_bin$exp_route=="proglottids",])

summary(total_prog)

##Obtain the confidence intervals
total.prog.ci <- confint.default(total_prog)


#EGGS
total_egg <- glm(formula = total_cyst~dose_a, family = "binomial", 
                  data = db_bin[db_bin$exp_route=="eggs",])

summary(total_egg)

##Obtain the confidence intervals
total.egg.ci <- confint.default(total_egg)


#BEETLES
total_bee <- glm(formula = total_cyst~dose_a, family = "binomial", 
                  data = db_bin[db_bin$exp_route=="beetles",])

summary(total_bee)

##Obtain the confidence intervals
total.bee.ci <- confint.default(total_bee)


#CAROTID
total_car <- glm(formula = total_cyst~dose_a, family = "binomial", 
                  data = db_bin[db_bin$exp_route=="carotid",])

summary(total_car)

##Obtain the confidence intervals
total.car.ci <- confint.default(total_car)



#Calculate the mean infective dose (ID50)
ID50.total.prog <- lr.ID50(total_prog$coefficients[1],total_prog$coefficients[2])
ID50.total.egg <- lr.ID50(total_egg$coefficients[1],total_egg$coefficients[2])
ID50.total.bee <- lr.ID50(total_bee$coefficients[1],total_bee$coefficients[2])
ID50.total.car <- lr.ID50(total_car$coefficients[1],total_car$coefficients[2])

##Create a data.frame to storage the calculated ID50
ID50.total <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"))
ID50.total$exp_route <- factor(ID50.total$exp_route, 
                                 levels = c("proglottids", "eggs", "beetles", "carotid"))
attach(ID50.total)
ID50.total$ID50 <- ifelse(exp_route=="proglottids",ID50.total.prog,
                         ifelse(exp_route=="carotid", ID50.total.car, 
                                ifelse(exp_route=="eggs", ID50.total.egg, ID50.total.bee)))
detach(ID50.total)


#PLOT DESIGN

##PROGLOTTIDS

###Extract the data for graphing the curve
curve.total.prog <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "proglottids"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.total.prog <- cbind(curve.total.prog, predict(total_prog, curve.total.prog, type = "link",
                                      se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.total.prog <- within(curve.total.prog, {
  p <- (1+exp(-fit))^-1
  LL <- (1+exp(-(fit - 1.96 * se.fit)))^-1
  UL <- (1+exp(-(fit + 1.96 * se.fit)))^-1
})


##EGGS

###Extract the data for graphing the curve
curve.total.egg <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "eggs"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.total.egg <- cbind(curve.total.egg, predict(total_egg, curve.total.egg, type = "link",
                                                    se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.total.egg <- within(curve.total.egg, {
  p <- (1+exp(-fit))^-1
  LL <- (1+exp(-(fit - 1.96 * se.fit)))^-1
  UL <- (1+exp(-(fit + 1.96 * se.fit)))^-1
})


##BEETLES

###Extract the data for graphing the curve
curve.total.bee <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "beetles"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.total.bee <- cbind(curve.total.bee, predict(total_bee, curve.total.bee, type = "link",
                                                    se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.total.bee <- within(curve.total.bee, {
  p <- (1+exp(-fit))^-1
  LL <- (1+exp(-(fit - 1.96 * se.fit)))^-1
  UL <- (1+exp(-(fit + 1.96 * se.fit)))^-1
})


##CAROTID

###Extract the data for graphing the curve
curve.total.car <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "carotid"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.total.car <- cbind(curve.total.car, predict(total_car, curve.total.car, type = "link",
                                                    se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.total.car <- within(curve.total.car, {
  p <- (1+exp(-fit))^-1
  LL <- (1+exp(-(fit - 1.96 * se.fit)))^-1
  UL <- (1+exp(-(fit + 1.96 * se.fit)))^-1
})


##MERGING THE DATA OF THE ROUTES OF EXPOSURE
curve.total <- rbind(curve.total.prog, curve.total.egg, curve.total.bee, curve.total.car)
curve.total$exp_route <- factor(curve.total$exp_route, 
                                levels = c("proglottids","eggs","beetles","carotid"))

#PLOT 1
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

#PROGLOTTIDS
ves_prog <- glm(formula = ves_cyst~dose_a, family = "binomial", 
                  data = db_bin[db_bin$exp_route=="proglottids",])

summary(ves_prog)

##Obtain the confidence intervals
ves.prog.ci <- confint.default(ves_prog)


#EGGS
ves_egg <- glm(formula = ves_cyst~dose_a, family = "binomial", 
                 data = db_bin[db_bin$exp_route=="eggs",])

summary(ves_egg)

##Obtain the confidence intervals
ves.egg.ci <- confint.default(ves_egg)


#BEETLES
ves_bee <- glm(formula = ves_cyst~dose_a, family = "binomial", 
                 data = db_bin[db_bin$exp_route=="beetles",])

summary(ves_bee)

##Obtain the confidence intervals
ves.bee.ci <- confint.default(ves_bee)


#CAROTID
ves_car <- glm(formula = ves_cyst~dose_a, family = "binomial", 
                 data = db_bin[db_bin$exp_route=="carotid",])

summary(ves_car)

##Obtain the confidence intervals
ves.car.ci <- confint.default(ves_car)



#Calculate the mean infective dose (ID50)
ID50.ves.prog <- lr.ID50(ves_prog$coefficients[1],ves_prog$coefficients[2])
ID50.ves.egg <- lr.ID50(ves_egg$coefficients[1],ves_egg$coefficients[2])
ID50.ves.bee <- lr.ID50(ves_bee$coefficients[1],ves_bee$coefficients[2])
ID50.ves.car <- lr.ID50(ves_car$coefficients[1],ves_car$coefficients[2])

##Create a data.frame to storage the calculated ID50
ID50.ves <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"))
ID50.ves$exp_route <- factor(ID50.ves$exp_route, 
                               levels = c("proglottids", "eggs", "beetles", "carotid"))
attach(ID50.ves)
ID50.ves$ID50 <- ifelse(exp_route=="proglottids",ID50.ves.prog,
                          ifelse(exp_route=="carotid", ID50.ves.car, 
                                 ifelse(exp_route=="eggs", ID50.ves.egg, ID50.ves.bee)))
detach(ID50.ves)


#PLOT DESIGN

##PROGLOTTIDS

###Extract the data for graphing the curve
curve.ves.prog <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "proglottids"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.ves.prog <- cbind(curve.ves.prog, predict(ves_prog, curve.ves.prog, type = "link",
                                                    se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.ves.prog <- within(curve.ves.prog, {
  p <- (1+exp(-fit))^-1
  LL <- (1+exp(-(fit - 1.96 * se.fit)))^-1
  UL <- (1+exp(-(fit + 1.96 * se.fit)))^-1
})


##EGGS

###Extract the data for graphing the curve
curve.ves.egg <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "eggs"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.ves.egg <- cbind(curve.ves.egg, predict(ves_egg, curve.ves.egg, type = "link",
                                                  se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.ves.egg <- within(curve.ves.egg, {
  p <- (1+exp(-fit))^-1
  LL <- (1+exp(-(fit - 1.96 * se.fit)))^-1
  UL <- (1+exp(-(fit + 1.96 * se.fit)))^-1
})


##BEETLES

###Extract the data for graphing the curve
curve.ves.bee <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "beetles"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.ves.bee <- cbind(curve.ves.bee, predict(ves_bee, curve.ves.bee, type = "link",
                                                  se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.ves.bee <- within(curve.ves.bee, {
  p <- (1+exp(-fit))^-1
  LL <- (1+exp(-(fit - 1.96 * se.fit)))^-1
  UL <- (1+exp(-(fit + 1.96 * se.fit)))^-1
})


##CAROTID

###Extract the data for graphing the curve
curve.ves.car <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "carotid"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.ves.car <- cbind(curve.ves.car, predict(ves_car, curve.ves.car, type = "link",
                                                  se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.ves.car <- within(curve.ves.car, {
  p <- (1+exp(-fit))^-1
  LL <- (1+exp(-(fit - 1.96 * se.fit)))^-1
  UL <- (1+exp(-(fit + 1.96 * se.fit)))^-1
})


##MERGING THE DATA OF THE ROUTES OF EXPOSURE
curve.ves <- rbind(curve.ves.prog, curve.ves.egg, curve.ves.bee, curve.ves.car)
curve.ves$exp_route <- factor(curve.ves$exp_route, 
                                levels = c("proglottids","eggs","beetles","carotid"))


#PLOT 1
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

#PROGLOTTIDS
brain_prog <- glm(formula = brain_cyst~dose_b, family = "binomial", 
                  data = db_bin[db_bin$exp_route=="proglottids",])

summary(brain_prog)

##Obtain the confidence intervals
brain.prog.ci <- confint.default(brain_prog)


#EGGS
brain_egg <- glm(formula = brain_cyst~dose_b, family = "binomial", 
                 data = db_bin[db_bin$exp_route=="eggs",])

summary(brain_egg)

##Obtain the confidence intervals
brain.egg.ci <- confint.default(brain_egg)


#BEETLES
brain_bee <- glm(formula = brain_cyst~dose_b, family = "binomial", 
                 data = db_bin[db_bin$exp_route=="beetles",])

summary(brain_bee)

##Obtain the confidence intervals
brain.bee.ci <- confint.default(brain_bee)


#CAROTID
brain_car <- glm(formula = brain_cyst~dose_b, family = "binomial", 
                 data = db_bin[db_bin$exp_route=="carotid",])

summary(brain_car)

##Obtain the confidence intervals
brain.car.ci <- confint.default(brain_car)



#Calculate the mean infective dose (ID50)
ID50.brain.prog <- lr.ID50(brain_prog$coefficients[1],brain_prog$coefficients[2])
ID50.brain.egg <- lr.ID50(brain_egg$coefficients[1],brain_egg$coefficients[2])
ID50.brain.bee <- lr.ID50(brain_bee$coefficients[1],brain_bee$coefficients[2])
ID50.brain.car <- lr.ID50(brain_car$coefficients[1],brain_car$coefficients[2])

##Create a data.frame to storage the calculated ID50
ID50.brain <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"))
ID50.brain$exp_route <- factor(ID50.brain$exp_route, 
                               levels = c("proglottids", "eggs", "beetles", "carotid"))
attach(ID50.brain)
ID50.brain$ID50 <- ifelse(exp_route=="proglottids",ID50.brain.prog,
                          ifelse(exp_route=="carotid", ID50.brain.car, 
                                 ifelse(exp_route=="eggs", ID50.brain.egg, ID50.brain.bee)))
detach(ID50.brain)


#PLOT DESIGN

##PROGLOTTIDS

###Extract the data for graphing the curve
curve.brain.prog <- data.frame(
  dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "proglottids"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.brain.prog <- cbind(curve.brain.prog, predict(brain_prog, curve.brain.prog, type = "link",
                                                    se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.brain.prog <- within(curve.brain.prog, {
  p <- (1+exp(-fit))^-1
  LL <- (1+exp(-(fit - 1.96 * se.fit)))^-1
  UL <- (1+exp(-(fit + 1.96 * se.fit)))^-1
})


##EGGS

###Extract the data for graphing the curve
curve.brain.egg <- data.frame(
  dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "eggs"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.brain.egg <- cbind(curve.brain.egg, predict(brain_egg, curve.brain.egg, type = "link",
                                                  se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.brain.egg <- within(curve.brain.egg, {
  p <- (1+exp(-fit))^-1
  LL <- (1+exp(-(fit - 1.96 * se.fit)))^-1
  UL <- (1+exp(-(fit + 1.96 * se.fit)))^-1
})


##BEETLES

###Extract the data for graphing the curve
curve.brain.bee <- data.frame(
  dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "beetles"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.brain.bee <- cbind(curve.brain.bee, predict(brain_bee, curve.brain.bee, type = "link",
                                                  se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.brain.bee <- within(curve.brain.bee, {
  p <- (1+exp(-fit))^-1
  LL <- (1+exp(-(fit - 1.96 * se.fit)))^-1
  UL <- (1+exp(-(fit + 1.96 * se.fit)))^-1
})


##CAROTID

###Extract the data for graphing the curve
curve.brain.car <- data.frame(
  dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "carotid"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.brain.car <- cbind(curve.brain.car, predict(brain_car, curve.brain.car, type = "link",
                                                  se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.brain.car <- within(curve.brain.car, {
  p <- (1+exp(-fit))^-1
  LL <- (1+exp(-(fit - 1.96 * se.fit)))^-1
  UL <- (1+exp(-(fit + 1.96 * se.fit)))^-1
})


##MERGING THE DATA OF THE ROUTES OF EXPOSURE
curve.brain <- rbind(curve.brain.prog, curve.brain.egg, curve.brain.bee, curve.brain.car)
curve.brain$exp_route <- factor(curve.brain$exp_route, 
                                levels = c("proglottids","eggs","beetles","carotid"))


#PLOT 1
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
