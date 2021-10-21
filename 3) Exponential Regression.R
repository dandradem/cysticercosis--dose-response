###############################################################################
##########################EXPONENTIAL REGRESSION###############################
###############################################################################

#Clean the R environment
rm(list=ls())


#Require necessary packages
if("rstudioapi" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("rstudioapi"); require(rstudioapi)} else{require(rstudioapi)}
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
exp.ID50 <- function(k) {
  ID50 <- -(log(0.5)/k)
  return(ID50)
}

#############
#TOTAL CYSTS#
#############

#PARAMETERS ESTIMATION

##Proglottids
total.prog <- nls(total_cyst~1-exp(-k*dose_a), 
                  data = db_bin[db_bin$exp_route=="proglottids",], 
                  start = list(k=0.0001))
sum.total.prog <- summary(total.prog)

##Eggs
total.egg <- nls(total_cyst~1-exp(-k*dose_a), 
                  data = db_bin[db_bin$exp_route=="eggs",], 
                  start = list(k=0.01109))
sum.total.egg <- summary(total.egg)

##Beetles
total.bee <- nls(total_cyst~1-exp(-k*dose_a), 
                 data = db_bin[db_bin$exp_route=="beetles",], 
                 start = list(k=0.01))
sum.total.bee <- summary(total.bee)

##Carotid
total.car <- nls(total_cyst~1-exp(-k*dose_a), 
                 data = db_bin[db_bin$exp_route=="carotid",], 
                 start = list(k=0.0001))
sum.total.car <- summary(total.car)

#Obtain the confidence intervals

total.prog.ci <- confint.default(total.prog)
total.egg.ci <- confint.default(total.egg)
total.bee.ci <- confint.default(total.bee)
total.car.ci <- confint.default(total.car)

#Calculate the mean infective dose (ID50)
ID50.total.prog <- exp.ID50(sum.total.prog$coefficients[1,1])
ID50.total.egg <- exp.ID50(sum.total.egg$coefficients[1,1])
ID50.total.bee <- exp.ID50(sum.total.bee$coefficients[1,1])
ID50.total.car <- exp.ID50(sum.total.car$coefficients[1,1])

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

#Extract the data for graphing the curve

##Proglottids
curve.total.prog <- expand.grid(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.total.prog$exp_route <- "proglottids"

curve.total.prog$p <- predict(total.prog,newdata = curve.total.prog)
curve.total.prog$LL <- 1-exp(-total.prog.ci[1]*curve.total.prog$dose_a)
curve.total.prog$UL <- 1-exp(-total.prog.ci[2]*curve.total.prog$dose_a)

##Eggs
curve.total.egg <- expand.grid(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.total.egg$exp_route <- "eggs"

curve.total.egg$p <- predict(total.egg,newdata = curve.total.egg)
curve.total.egg$LL <- 1-exp(-total.egg.ci[1]*curve.total.egg$dose)
curve.total.egg$UL <- 1-exp(-total.egg.ci[2]*curve.total.egg$dose)

##Beetles
curve.total.bee <- expand.grid(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.total.bee$exp_route <- "beetles"

curve.total.bee$p <- predict(total.bee,newdata = curve.total.bee)
curve.total.bee$LL <- 1-exp(-total.bee.ci[1]*curve.total.bee$dose_a)
curve.total.bee$UL <- 1-exp(-total.bee.ci[2]*curve.total.bee$dose_a)

###CAROTIDA
curve.total.car <- expand.grid(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.total.car$exp_route <- "carotid"

curve.total.car$p <- predict(total.car,newdata = curve.total.car,interval = "confidence")
curve.total.car$LL <- 1-exp(-total.car.ci[1]*curve.total.car$dose_a)
curve.total.car$UL <- 1-exp(-total.car.ci[2]*curve.total.car$dose_a)

#Merge the predicted values from each route of exposure
curve.total <- rbind(curve.total.prog,curve.total.car,curve.total.egg,
                     curve.total.bee)
curve.total$exp_route <- factor(curve.total$exp_route,
                                levels = c("proglottids", "eggs", "beetles", "carotid"))

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

#PARAMETERS ESTIMATION

##Proglottids
ves.prog <- nls(ves_cyst~1-exp(-k*dose_a), 
                  data = db_bin[db_bin$exp_route=="proglottids",], 
                  start = list(k=0.0001))
sum.ves.prog <- summary(ves.prog)

##Eggs
ves.egg <- nls(ves_cyst~1-exp(-k*dose_a), 
                 data = db_bin[db_bin$exp_route=="eggs",], 
                 start = list(k=0.01109))
sum.ves.egg <- summary(ves.egg)

##Beetles
ves.bee <- nls(ves_cyst~1-exp(-k*dose_a), 
                 data = db_bin[db_bin$exp_route=="beetles",], 
                 start = list(k=0.01))
sum.ves.bee <- summary(ves.bee)

##Carotid
ves.car <- nls(ves_cyst~1-exp(-k*dose_a), 
                 data = db_bin[db_bin$exp_route=="carotid",], 
                 start = list(k=0.0001))
sum.ves.car <- summary(ves.car)

#Obtain the confidence intervals

ves.prog.ci <- confint.default(ves.prog)
ves.egg.ci <- confint.default(ves.egg)
ves.bee.ci <- confint.default(ves.bee)
ves.car.ci <- confint.default(ves.car)

#Calculate the mean infective dose (ID50)
ID50.ves.prog <- exp.ID50(sum.ves.prog$coefficients[1,1])
ID50.ves.egg <- exp.ID50(sum.ves.egg$coefficients[1,1])
ID50.ves.bee <- exp.ID50(sum.ves.bee$coefficients[1,1])
ID50.ves.car <- exp.ID50(sum.ves.car$coefficients[1,1])

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

#Extract the data for graphing the curve

##Proglottids
curve.ves.prog <- expand.grid(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.ves.prog$exp_route <- "proglottids"

curve.ves.prog$p <- predict(ves.prog,newdata = curve.ves.prog)
curve.ves.prog$LL <- 1-exp(-ves.prog.ci[1]*curve.ves.prog$dose_a)
curve.ves.prog$UL <- 1-exp(-ves.prog.ci[2]*curve.ves.prog$dose_a)

##Eggs
curve.ves.egg <- expand.grid(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.ves.egg$exp_route <- "eggs"

curve.ves.egg$p <- predict(ves.egg,newdata = curve.ves.egg)
curve.ves.egg$LL <- 1-exp(-ves.egg.ci[1]*curve.ves.egg$dose)
curve.ves.egg$UL <- 1-exp(-ves.egg.ci[2]*curve.ves.egg$dose)

##Beetles
curve.ves.bee <- expand.grid(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.ves.bee$exp_route <- "beetles"

curve.ves.bee$p <- predict(ves.bee,newdata = curve.ves.bee)
curve.ves.bee$LL <- 1-exp(-ves.bee.ci[1]*curve.ves.bee$dose_a)
curve.ves.bee$UL <- 1-exp(-ves.bee.ci[2]*curve.ves.bee$dose_a)

###CAROTIDA
curve.ves.car <- expand.grid(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.ves.car$exp_route <- "carotid"

curve.ves.car$p <- predict(ves.car,newdata = curve.ves.car,interval = "confidence")
curve.ves.car$LL <- 1-exp(-ves.car.ci[1]*curve.ves.car$dose_a)
curve.ves.car$UL <- 1-exp(-ves.car.ci[2]*curve.ves.car$dose_a)

#Merge the predicted values from each route of exposure
curve.ves <- rbind(curve.ves.prog,curve.ves.car,curve.ves.egg,
                     curve.ves.bee)
curve.ves$exp_route <- factor(curve.ves$exp_route,
                                levels = c("proglottids", "eggs", "beetles", "carotid"))

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

#PARAMETERS ESTIMATION

##Proglottids
brain.prog <- nls(brain_cyst~1-exp(-k*dose_b), 
                  data = db_bin[db_bin$exp_route=="proglottids",], 
                  start = list(k=0.0001))
sum.brain.prog <- summary(brain.prog)

##Eggs
brain.egg <- nls(brain_cyst~1-exp(-k*dose_b), 
                 data = db_bin[db_bin$exp_route=="eggs",], 
                 start = list(k=0.0001))
sum.brain.egg <- summary(brain.egg)

##Beetles
brain.bee <- nls(brain_cyst~1-exp(-k*dose_b), 
                 data = db_bin[db_bin$exp_route=="beetles",], 
                 start = list(k=0.001))
sum.brain.bee <- summary(brain.bee)

##Carotid
brain.car <- nls(brain_cyst~1-exp(-k*dose_b), 
                 data = db_bin[db_bin$exp_route=="carotid",], 
                 start = list(k=0.0001))
sum.brain.car <- summary(brain.car)

#Obtain the confidence intervals

brain.prog.ci <- confint.default(brain.prog)
brain.egg.ci <- confint.default(brain.egg)
brain.bee.ci <- confint.default(brain.bee)
brain.car.ci <- confint.default(brain.car)

#Calculate the mean infective dose (ID50)
ID50.brain.prog <- exp.ID50(sum.brain.prog$coefficients[1,1])
ID50.brain.egg <- exp.ID50(sum.brain.egg$coefficients[1,1])
ID50.brain.bee <- exp.ID50(sum.brain.bee$coefficients[1,1])
ID50.brain.car <- exp.ID50(sum.brain.car$coefficients[1,1])

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

#Extract the data for graphing the curve

##Proglottids
curve.brain.prog <- expand.grid(dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.brain.prog$exp_route <- "proglottids"

curve.brain.prog$p <- predict(brain.prog,newdata = curve.brain.prog)
curve.brain.prog$LL <- 1-exp(-brain.prog.ci[1]*curve.brain.prog$dose_b)
curve.brain.prog$UL <- 1-exp(-brain.prog.ci[2]*curve.brain.prog$dose_b)

##Eggs
curve.brain.egg <- expand.grid(dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.brain.egg$exp_route <- "eggs"

curve.brain.egg$p <- predict(brain.egg,newdata = curve.brain.egg)
curve.brain.egg$LL <- 1-exp(-brain.egg.ci[1]*curve.brain.egg$dose)
curve.brain.egg$UL <- 1-exp(-brain.egg.ci[2]*curve.brain.egg$dose)

##Beetles
curve.brain.bee <- expand.grid(dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.brain.bee$exp_route <- "beetles"

curve.brain.bee$p <- predict(brain.bee,newdata = curve.brain.bee)
curve.brain.bee$LL <- 1-exp(-brain.bee.ci[1]*curve.brain.bee$dose_b)
curve.brain.bee$UL <- 1-exp(-brain.bee.ci[2]*curve.brain.bee$dose_b)

###CAROTIDA
curve.brain.car <- expand.grid(dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.brain.car$exp_route <- "carotid"

curve.brain.car$p <- predict(brain.car,newdata = curve.brain.car,interval = "confidence")
curve.brain.car$LL <- 1-exp(-brain.car.ci[1]*curve.brain.car$dose_b)
curve.brain.car$UL <- 1-exp(-brain.car.ci[2]*curve.brain.car$dose_b)

#Merge the predicted values from each route of exposure
curve.brain <- rbind(curve.brain.prog,curve.brain.car,curve.brain.egg,
                     curve.brain.bee)
curve.brain$exp_route <- factor(curve.brain$exp_route,
                                levels = c("proglottids", "eggs", "beetles", "carotid"))

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
