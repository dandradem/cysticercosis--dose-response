###############################################################################
########################CODE FOR CREATION OF DATASETS##########################
###############################################################################

#Require the necessary packages
if("tidyselect" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("tidyselect"); require(tidyselect)} else{require(tidyselect)}
if("dplyr" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("dplyr"); require(dplyr)} else{require(dplyr)}
if("tidyverse" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("tidyverse"); require(tidyverse)} else{require(tidyverse)}
if("extraDistr" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("extraDistr"); require(extraDistr)} else{require(extraDistr)}

#Import dataset

db <- read.csv2("dataset.csv", header=T)
db <- db[,-c(1,2,4)]
db$dose <- as.numeric(db$dose)

#Set the factor levels for the routes of exposure
db$exp_route <- factor(db$exp_route, 
                        levels = c("proglottids", "eggs", "beetles", "carotid"))

####################################
#ASSUMPTION OF 15 NEGATIVE CONTROLS#
####################################

#Eliminate the existing negative controls
attach(db)
db <- db[dose!=0,]
detach(db)

#Create the data of the 15 negative controls for each route of exposure
neg_control <- data.frame(exp_route=c(rep("eggs",15),rep("proglottids",15),
                                      rep("beetles",15), rep("carotid",15)),
                          dose=0, ves_cyst=0, deg_cyst=0, total_cyst=0, 
                          brain_cyst=0)

#Add the negative controls to the database
db <- rbind(db, neg_control)



###############################
#STANDARDIZATION OF DOSE UNITS#
###############################

##Set seeds and create the database for standardization
set.seed(2)

beetles <- c(28, 11, 17, 61, 58, 59, 121, 14, 17, 40, 235, 19, 21, 
               73, 21, 15, 36, 20, 44, 37, 68, 41, 7, 56, 8, 5, 5, 31,
               9, 60, 6, 4, 20, 3, 31)

proglottids <- rdunif(10000,30000,50000)


##Create the function for sampling the dose of T. solium eggs inoculated to each
##pig by its dose group
egg.sampling <- function(data, dose) {
  sam <- sample(data,dose,replace = T)
  return(sum(sam))
}

##Make a for loop to make the sampling for each pig involved in the studies

###BEETLES
d1 <- vector(length = 6)
for(n in 1:6){
  d1[n] <- egg.sampling(beetles, 1)
}


d3 <- vector(length = 6)
for(n in 1:6){
  d3[n] <- egg.sampling(beetles, 3)
}


d4 <- vector(length = 8)
for(n in 1:8){
  d4[n] <- egg.sampling(beetles, 4)
}


d6 <- vector(length = 30)
for(n in 1:30){
  d6[n] <- egg.sampling(beetles, 6)
}

####Calculate the mean and standard deviation for each group dose
beedose1 <- round(mean(d1))
beedose3 <- round(mean(d3))
beedose4 <- round(mean(d4))
beedose6 <- round(mean(d6))

beesd1 <- sd(d1)
beesd3 <- sd(d3)
beesd4 <- sd(d4)
beesd6 <- sd(d6)


###PROGLOTTIDS
d0.25 <- vector(length = 8)
for(n in 1:8){
  d0.25[n] <- egg.sampling(proglottids, 1)
}
d0.25 <- round(d0.25*0.25)

d0.5 <- vector(length = 8)
for(n in 1:8){
  d0.5[n] <- egg.sampling(proglottids, 1)
}
d0.5 <- round(d0.5*0.5)

####The following group dose is particular because there were pigs which were not
####involved in brain cysts counting, so create a specific group dose for
####the pigs which get brain cysts data

d1.0 <- vector(length = 30)
for(n in 1:30){
  d1.0[n] <- egg.sampling(proglottids, 1)
}

d1.0_b <- vector(length = 14)
for(n in 1:14){
  d1.0_b[n] <- egg.sampling(proglottids, 1)
}

####Calculate the mean and standard deviation for each group dose
prodose0.25 <- round(mean(d0.25))
prodose0.5 <- round(mean(d0.5))
prodose1.0 <- round(mean(d1.0))
prodose1.0_b <- round(mean(d1.0_b))

prosd0.25 <- sd(d0.25)
prosd0.5 <- sd(d0.5)
prosd1.0 <- sd(d1.0)
prosd1.0_b <- sd(d1.0_b)

###Add the standardized doses values to the data

db$dose_a <- ifelse(db$exp_route=="eggs", db$dose,
                      ifelse(db$exp_route=="carotid", db$dose, 
                             ifelse(db$exp_route=="proglottids", ifelse(db$dose==0.25,prodose0.25,
                                                                       ifelse(db$dose==0.5,prodose0.5,
                                                                              ifelse(db$dose==1,prodose1.0,0))),
                                    ifelse(db$exp_route=="beetles",ifelse(db$dose==1,beedose1,
                                                                               ifelse(db$dose==3,beedose3,
                                                                                      ifelse(db$dose==4,beedose4,
                                                                                             ifelse(db$dose==6,beedose6,0)))),0))))

db$dose_sd_a <- ifelse(db$exp_route=="proglottids", ifelse(db$dose==0.25, prosd0.25,
                                                            ifelse(db$dose==0.5, prosd0.5,
                                                                   ifelse(db$dose==1, prosd1.0,NA))),
                         ifelse(db$exp_route=="beetles", ifelse(db$dose==1, beesd1,
                                                                     ifelse(db$dose==3, beesd3,
                                                                            ifelse(db$dose==4, beesd4,
                                                                                   ifelse(db$dose==6, beesd6,NA)))),NA))


db$dose_b <- ifelse(db$exp_route=="eggs", db$dose,
                    ifelse(db$exp_route=="carotid", db$dose, 
                           ifelse(db$exp_route=="proglottids", ifelse(db$dose==0.25,prodose0.25,
                                                                      ifelse(db$dose==0.5,prodose0.5,
                                                                             ifelse(db$dose==1,prodose1.0_b,0))),
                                  ifelse(db$exp_route=="beetles",ifelse(db$dose==1,beedose1,
                                                                        ifelse(db$dose==3,beedose3,
                                                                               ifelse(db$dose==4,beedose4,
                                                                                      ifelse(db$dose==6,beedose6,0)))),0))))

db$dose_sd_b <- ifelse(db$exp_route=="proglottids", ifelse(db$dose==0.25, prosd0.25,
                                                           ifelse(db$dose==0.5, prosd0.5,
                                                                  ifelse(db$dose==1, prosd1.0_b,NA))),
                       ifelse(db$exp_route=="beetles", ifelse(db$dose==1, beesd1,
                                                              ifelse(db$dose==3, beesd3,
                                                                     ifelse(db$dose==4, beesd4,
                                                                            ifelse(db$dose==6, beesd6,NA)))),NA))


###########################
#DATABASE WITH CYSTS COUNT#
###########################

db <- db[,c(1,7,8,3,4,5,9,10,6)]



#################################
#DATABASE WITH BINOMIAL RESPONSE#
#################################

#BASE DE DATOS CON VARIABLE CUALITATIVA BINARIA
db_bin <- db[,-5]

###Volvemos a la variable "quistes totales" y "quistes cerebrales" una variable binaria
db_bin$total_cyst <- ifelse(db_bin$total_cyst>=1, 1, 0)
db_bin$ves_cyst <- ifelse(db_bin$ves_cyst>=1, 1, 0)
db_bin$brain_cyst <- ifelse(db_bin$brain_cyst>=1, 1, 0)



##############################################
#DATABASE WITH NUMBER OF EXPOSED AND INFECTED#
##############################################

#Create a dataframe with exposed and infected pigs for each type of developed
#cysts

db_exp_total <- db %>% 
  dplyr::select(exp_route,dose_a,dose_sd_a,total_cyst) %>%
  group_by(exp_route,dose_a) %>%
  mutate(nExp = length(total_cyst)) %>%
  filter(total_cyst>0) %>%
  mutate(nInf = length(total_cyst)) %>%
  dplyr::select(-total_cyst) %>%
  distinct() %>%
  as.data.frame() %>%
  add_row(exp_route = "proglottids", dose_a = 0, dose_sd_a = NA, nExp = 15, nInf = 0) %>%
  add_row(exp_route = "eggs", dose_a = 0, dose_sd_a = NA, nExp = 15, nInf = 0) %>%
  add_row(exp_route = "beetles", dose_a = 0, dose_sd_a = NA, nExp = 15, nInf = 0) %>%
  add_row(exp_route = "carotid", dose_a = 0, dose_sd_a = NA, nExp = 15, nInf = 0) %>%
  arrange(exp_route,dose_a)

db_exp_total$exp_route <- factor(db_exp_total$exp_route, 
                               levels = c("proglottids", "eggs", "beetles", "carotid"))

db_exp_ves <- db %>% 
  dplyr::select(exp_route,dose_a,dose_sd_a,ves_cyst) %>%
  group_by(exp_route,dose_a) %>%
  mutate(nExp = length(ves_cyst)) %>%
  filter(ves_cyst>0) %>%
  mutate(nInf = length(ves_cyst)) %>%
  dplyr::select(-ves_cyst) %>%
  distinct() %>%
  as.data.frame() %>%
  add_row(exp_route = "proglottids", dose_a = 0, dose_sd_a = NA, nExp = 15, nInf = 0) %>%
  add_row(exp_route = "eggs", dose_a = 0, dose_sd_a = NA, nExp = 15, nInf = 0) %>%
  add_row(exp_route = "beetles", dose_a = 0, dose_sd_a = NA, nExp = 15, nInf = 0) %>%
  add_row(exp_route = "carotid", dose_a = 0, dose_sd_a = NA, nExp = 15, nInf = 0) %>%
  arrange(exp_route,dose_a)

db_exp_ves$exp_route <- factor(db_exp_ves$exp_route, 
                               levels = c("proglottids", "eggs", "beetles", "carotid"))

db_exp_brain <- db[!is.na(db$brain_cyst),] %>% 
  dplyr::select(exp_route,dose_b,dose_sd_b,brain_cyst) %>%
  group_by(exp_route,dose_b) %>%
  mutate(nExp = length(brain_cyst)) %>%
  filter(brain_cyst>0) %>%
  mutate(nInf = length(brain_cyst)) %>%
  dplyr::select(-brain_cyst) %>%
  distinct() %>%
  as.data.frame() %>%
  add_row(exp_route = "proglottids", dose_b = 0, dose_sd_b = NA, nExp = 15, nInf = 0) %>%
  add_row(exp_route = "eggs", dose_b = 0, dose_sd_b = NA, nExp = 15, nInf = 0) %>%
  add_row(exp_route = "beetles", dose_b = 0, dose_sd_b = NA, nExp = 15, nInf = 0) %>%
  add_row(exp_route = "carotid", dose_b = 0, dose_sd_b = NA, nExp = 15, nInf = 0) %>%
  add_row(exp_route = "beetles", dose_b = beedose4, dose_sd_b = beesd4, nExp = 8, nInf = 0) %>%
  add_row(exp_route = "eggs", dose_b = 10, dose_sd_b = NA, nExp = 4, nInf = 0) %>%
  add_row(exp_route = "eggs", dose_b = 100, dose_sd_b = NA, nExp = 9, nInf = 0) %>%
  add_row(exp_route = "eggs", dose_b = 1000, dose_sd_b = NA, nExp = 11, nInf = 0) %>%
  arrange(exp_route,dose_b)

db_exp_brain$exp_route <- factor(db_exp_brain$exp_route, 
                               levels = c("proglottids", "eggs", "beetles", "carotid"))

#Create the data.frames to evaluate the oral inoculation

db_oral_total <- db_exp_total[db_exp_total$dose_a!=0,]
db_oral_total$exp_route_b <- ifelse(db_oral_total$exp_route=="carotid", "carotid",
                                    "oral")
db_oral_total <- db_oral_total[,c(6,2:5)] %>%
  add_row(exp_route_b = "oral", dose_a = 0, dose_sd_a = NA, nExp = 15, nInf = 0) %>%
  add_row(exp_route_b = "carotid", dose_a = 0, dose_sd_a = NA, nExp = 15, nInf = 0)

db_oral_total$exp_route_b <- factor(db_oral_total$exp_route_b, 
                                    levels = c("oral", "carotid"))


db_oral_ves <- db_exp_ves[db_exp_ves$dose_a!=0,]
db_oral_ves$exp_route_b <- ifelse(db_oral_ves$exp_route=="carotid", "carotid",
                                  "oral")
db_oral_ves <- db_oral_ves[,c(6,2:5)] %>%
  add_row(exp_route_b = "oral", dose_a = 0, dose_sd_a = NA, nExp = 15, nInf = 0) %>%
  add_row(exp_route_b = "carotid", dose_a = 0, dose_sd_a = NA, nExp = 15, nInf = 0)

db_oral_ves$exp_route_b <- factor(db_oral_ves$exp_route_b, 
                                    levels = c("oral", "carotid"))


db_oral_brain <- db_exp_brain[db_exp_brain$dose_b!=0,]
db_oral_brain$exp_route_b <- ifelse(db_oral_brain$exp_route=="carotid", "carotid",
                                    "oral")
db_oral_brain <- db_oral_brain[,c(6,2:5)] %>%
  add_row(exp_route_b = "oral", dose_b = 0, dose_sd_b = NA, nExp = 15, nInf = 0) %>%
  add_row(exp_route_b = "carotid", dose_b = 0, dose_sd_b = NA, nExp = 15, nInf = 0)

db_oral_brain$exp_route_b <- factor(db_oral_brain$exp_route_b, 
                                    levels = c("oral", "carotid"))

##################################
#DATABASE FOR GRAPHING GEOM-POINT#
##################################

geom_total <- db_exp_total[db_exp_total$dose_a!=0,]
geom_total$exp_route_b <- ifelse(geom_total$exp_route=="carotid","carotid",
                                 "oral")
geom_total$exp_route_b <- factor(geom_total$exp_route_b, 
                                 levels = c("oral", "carotid"))


geom_ves <- db_exp_ves[db_exp_ves$dose_a!=0,]
geom_ves$exp_route_b <- ifelse(geom_ves$exp_route=="carotid","carotid",
                                 "oral")
geom_ves$exp_route_b <- factor(geom_ves$exp_route_b, 
                                 levels = c("oral", "carotid"))

geom_brain <- db_exp_brain[db_exp_brain$dose_b!=0,]
geom_brain$exp_route_b <- ifelse(geom_brain$exp_route=="carotid","carotid",
                                 "oral")
geom_brain$exp_route_b <- factor(geom_brain$exp_route_b, 
                                 levels = c("oral", "carotid"))
