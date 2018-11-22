setwd(dirname(parent.frame(2)$ofile))

library(foreign)
library(rgenoud)
library(Matching)
library(rbounds)

multiple.genmatch = function(
  dataset, outcome.label, match.formula,
  random.runs=1, pop.size=100,
  max.generations=50, wait.generations=20,
  thread.count=4)
{
  outcome.vec = dataset[, outcome.label]
  var.labels = all.vars(match.formula)
  treatment.label = var.labels[1]
  treatment.vec = dataset[, treatment.label]
  covariates.labels = var.labels[2:length(var.labels)]
  
  covariates = c()
  for (i in c(1:length(covariates.labels))) {
    covariates = cbind(covariates, dataset[, covariates.labels[i]])
  }
  covariates = cbind(data.matrix(covariates))
  
  max.AMsmallest.p.value = -1
  current.best.result = list()
  for (i in c(1:random.runs)) {
    genmatch.result <- GenMatch(
      Tr=treatment.vec,
      X=covariates,
      BalanceMatrix=covariates,
      estimand="ATT",
      M=1, pop.size=pop.size,
      max.generations=max.generations,
      wait.generations=wait.generations,
      cluster=rep('localhost', thread.count)
    )
    match.result <- Match(
      Y=outcome.vec,
      Tr=treatment.vec,
      X=covariates,
      estimand="ATT",
      caliper=2,
      Weight.matrix=genmatch.result
    )
    balance <- MatchBalance(
      match.formula,
      data=dataset,
      match.out=match.result,
      nboots=500
    )
    if (balance$AMsmallest.p.value > max.AMsmallest.p.value) {
      current.best.result = list(genmatch.result, match.result, balance)
      max.AMsmallest.p.value = balance$AMsmallest.p.value
    }
  }
  return(current.best.result)
}



#################### import and clean dataset ##################

data.2014 <- read.dta("Primary_School_Data_2014.dta")
data.2014$JC_JE <- ifelse(data.2014$JC_JE == 'Si', 1, 0)
data.2014$JC <- ifelse(data.2014$JC == 'Si', 1, 0)
data.2014$JE <- ifelse(data.2014$JE == 'Si', 1, 0)

# using original Provincias for stratification
data.2014.nona.regional <- na.omit(data.2014[, c('Sector','Ambito','Provincia',  
                                                 'Frcn_QuintilIVSHogares',
                                                 'Frcn_QuintilNoAsist4a17', 
                                                 'Duracion', 'JC_JE', 'Mat_Total', 
                                                 'AxS', 'Sec_Total', 'TasaPromov11',
                                                 'TasaPromov12', 'TasaPromov13')])

####### coarsen regions, all data #######
data.2014$Provincia <- as.numeric(data.2014$Provincia)

#Pampas: Córdoba, Santa Fe, La Pampa, Buenos Aires, Ciudad de Buenos Aires == 1
data.2014$Provincia[data.2014$Provincia == 2] <- 1 #Buenos Aires
data.2014$Provincia[data.2014$Provincia == 4] <- 1 #Córdoba
data.2014$Provincia[data.2014$Provincia == 11] <- 1 #La Pampa
data.2014$Provincia[data.2014$Provincia == 21] <- 1 #Santa Fe

#Argentine Northwest: Jujuy, Salta, Tucumán, Catamarca == 3
data.2014$Provincia[data.2014$Provincia == 10] <- 3 #Jujuy
data.2014$Provincia[data.2014$Provincia == 17] <- 3 #Salta
data.2014$Provincia[data.2014$Provincia == 23] <- 3 #Tucuman

#Gran Chaco: Formosa, Chaco, Santiago del Estero == 6
data.2014$Provincia[data.2014$Provincia == 9] <- 6 #Formosa
data.2014$Provincia[data.2014$Provincia == 22] <- 6 #Santiago del Estero

#Mesopotamia (or Littoral): Misiones, Entre Ríos, Corrientes == 5
data.2014$Provincia[data.2014$Provincia == 8] <- 5 #Entre Rios
data.2014$Provincia[data.2014$Provincia == 14] <- 5 #Misiones

#Cuyo: San Juan, La Rioja, Mendoza, San Luis == 2
data.2014$Provincia[data.2014$Provincia == 13] <- 2 #Mendoza
data.2014$Provincia[data.2014$Provincia == 18] <- 2 #San Juan
data.2014$Provincia[data.2014$Provincia == 19] <- 2 #San Luis
data.2014$Provincia[data.2014$Provincia == 12] <- 2 #La Rioja

#Patagonia: Rio Negro, Neuquén, Chubut, Santa Cruz, Tierra del Fuego == 4
data.2014$Provincia[data.2014$Provincia == 16] <- 4 #Rio Negro
data.2014$Provincia[data.2014$Provincia == 15] <- 4 #Neuquen
data.2014$Provincia[data.2014$Provincia == 7] <- 4 #Chubut
data.2014$Provincia[data.2014$Provincia == 20] <- 4 #Santa Cruz
data.2014$Provincia[data.2014$Provincia == 24] <- 4 #Tierra del Fuego

data.2014$Provincia <- factor(data.2014$Provincia)
####### /coarsen regions, all data #######

data.2014.subset12 = subset(data.2014, data.2014$ExisteJEJC11==0)
data.2014.subset13 = subset(data.2014, data.2014$ExisteJEJC11==0 & data.2014$ExisteJEJC12==0)

data.2014.nona <- na.omit(data.2014[, c('Sector','Ambito','Provincia',
                                        'Frcn_QuintilIVSHogares',
                                        'Frcn_QuintilNoAsist4a17', 
                                        'Duracion', 'JC_JE', 'JC', 
                                        'JE','Mat_Total', 
                                        'AxS', 'Sec_Total', 'TasaPromov11',
                                        'TasaPromov12', 'TasaPromov13')])
data.2014.subset12.nona <- na.omit(data.2014.subset12[, c('Sector','Ambito','Provincia',
                                                          'Frcn_QuintilIVSHogares',
                                                          'Frcn_QuintilNoAsist4a17', 
                                                          'Duracion', 'JC_JE', 'JC', 
                                                          'JE','Mat_Total', 
                                                          'AxS', 'Sec_Total', 'TasaPromov11',
                                                          'TasaPromov12', 'TasaPromov13')])
data.2014.subset13.nona <- na.omit(data.2014.subset13[, c('Sector','Ambito','Provincia',
                                                          'Frcn_QuintilIVSHogares',
                                                          'Frcn_QuintilNoAsist4a17', 
                                                          'Duracion', 'JC_JE', 'JC', 
                                                          'JE','Mat_Total', 
                                                          'AxS', 'Sec_Total', 'TasaPromov11',
                                                          'TasaPromov12', 'TasaPromov13', 'TasaPromov14')])

#################### /import and clean dataset ##################

result.JC = multiple.genmatch(
  data.2014.nona,
  "TasaPromov12",
  JC~TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JE = multiple.genmatch(
  data.2014.nona,
  "TasaPromov12",
  JE~TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JCJE.11.12 = multiple.genmatch(
  data.2014.nona,
  "TasaPromov12",
  JC_JE~TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JC.12.13 = multiple.genmatch(
  data.2014.subset12.nona,
  "TasaPromov13",
  JC~TasaPromov12+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JE.12.13 = multiple.genmatch(
  data.2014.subset12.nona,
  "TasaPromov13",
  JE~TasaPromov12+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JCJE.12.13 = multiple.genmatch(
  data.2014.subset12.nona,
  "TasaPromov13",
  JC_JE~TasaPromov12+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JC.13.14 = multiple.genmatch(
  data.2014.subset13.nona,
  "TasaPromov14",
  JC~TasaPromov13+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JE.13.14 = multiple.genmatch(
  data.2014.subset13.nona,
  "TasaPromov14",
  JE~TasaPromov13+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JCJE.13.14 = multiple.genmatch(
  data.2014.subset13.nona,
  "TasaPromov14",
  JC_JE~TasaPromov13+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

######## Stratified Sector Effect (private vs. public) on treatment effect #########

data.2014.nona2 <- data.2014.nona
data.2014.nona2$Sector <- ifelse(data.2014.nona2$Sector == 'Privado', 1, 0)

Privado.data.2014.nona2 <- subset(data.2014.nona2, data.2014.nona2$Sector == 1) 
Estatal.data.2014.nona2 <- subset(data.2014.nona2, data.2014.nona2$Sector == 0)  

#Privado
result.JCJE.privado = multiple.genmatch(
  Privado.data.2014.nona2,
  "TasaPromov12",
  JC_JE ~ TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JC.privado = multiple.genmatch(
  Privado.data.2014.nona2,
  "TasaPromov12",
  JC ~ TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JE.privado = multiple.genmatch(
  Privado.data.2014.nona2,
  "TasaPromov12",
  JE ~ TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

#Estatal
result.JCJE.estatal = multiple.genmatch(
  Estatal.data.2014.nona2,
  "TasaPromov12",
  JC_JE ~ TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JC.estatal = multiple.genmatch(
  Estatal.data.2014.nona2,
  "TasaPromov12",
  JC ~ TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JE.estatal = multiple.genmatch(
  Estatal.data.2014.nona2,
  "TasaPromov12",
  JE ~ TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

######## /Stratified Sector Effect (private vs. public) on treatment effect #########

######## Stratified Ambito Effect (Urban vs. Rural) on treatment effect ########
data.2014.nona4 <- data.2014.nona
data.2014.nona4$Ambito <- ifelse(data.2014.nona4$Ambito == 'Urbano', 1, 0)

Urbano.data.2014.nona4 <- subset(data.2014.nona4, data.2014.nona4$Ambito == 1) 
Rural.data.2014.nona4 <- subset(data.2014.nona4, data.2014.nona4$Ambito == 0)

#Urbano
result.JCJE.urbano = multiple.genmatch(
  Urbano.data.2014.nona4,
  "TasaPromov12",
  JC_JE ~ TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JC.urbano = multiple.genmatch(
  Urbano.data.2014.nona4,
  "TasaPromov12",
  JC ~ TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JE.urbano = multiple.genmatch(
  Urbano.data.2014.nona4,
  "TasaPromov12",
  JE ~ TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

#Rural
result.JCJE.rural = multiple.genmatch(
  Rural.data.2014.nona4,
  "TasaPromov12",
  JC_JE ~ TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JC.rural = multiple.genmatch(
  Rural.data.2014.nona4,
  "TasaPromov12",
  JC ~ TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

result.JE.rural = multiple.genmatch(
  Rural.data.2014.nona4,
  "TasaPromov12",
  JE ~ TasaPromov11+Provincia+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)
######## /Stratified Ambito Effect (Urban vs. Rural) on treatment effect ########

######## Stratified Regional Effects ##########

data.2014.nona.regional$JC_JE <- ifelse(data.2014.nona.regional$JC_JE == 'Si', 1, 0)
data.2014.nona.regional$Provincia = as.numeric(data.2014.nona.regional$Provincia)

result.CBA = multiple.genmatch(
  CBA.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)

CBA.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 1) 
BA.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 2)
Catamarca.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 3)
Córdoba.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 4)
Corrientes.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 5)
Chaco.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 6)
Chubut.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 7)
Entre.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 8)
Formosa.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 9)
Jujuy.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 10)
LaPampa.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 11)
LaRioja.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 12)

Mendoza.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 13)
Misiones.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 14)
Neuquén.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 15)
Negro.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 16)
Salta.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 17)
SanJuan.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 18)
SanLuis.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 19)
SantaCruz.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 20)
SantaFe.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 21)
Santiago.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 22)
Tucumán.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 23)
Tierra.data.2014.nona.regional <- subset(data.2014.nona.regional, data.2014.nona.regional$Provincia == 24)

result.CBA = multiple.genmatch(
  CBA.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.BA = multiple.genmatch(
  BA.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Catamarca = multiple.genmatch(
  Catamarca.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Córdoba = multiple.genmatch(
  Córdoba.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Corrientes = multiple.genmatch(
  Corrientes.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Chaco = multiple.genmatch(
  Chaco.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Chubut = multiple.genmatch(
  Chubut.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Entre = multiple.genmatch(
  Entre.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Formosa = multiple.genmatch(
  Formosa.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Jujuy = multiple.genmatch(
  Jujuy.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.LaPampa = multiple.genmatch(
  LaPampa.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.LaRioja = multiple.genmatch(
  LaRioja.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Mendoza = multiple.genmatch(
  Mendoza.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Misiones = multiple.genmatch(
  Misiones.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Neuquén = multiple.genmatch(
  Neuquén.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Negro = multiple.genmatch(
  Negro.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Salta = multiple.genmatch(
  Salta.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.SanJuan = multiple.genmatch(
  SanJuan.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.SanLuis = multiple.genmatch(
  SanLuis.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.SantaCruz = multiple.genmatch(
  SantaCruz.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.SantaFe = multiple.genmatch(
  SantaFe.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Santiago = multiple.genmatch(
  Santiago.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Tucumán = multiple.genmatch(
  Tucumán.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


result.Tierra = multiple.genmatch(
  Tierra.data.2014.nona.regional,
  "TasaPromov12",
  JC_JE~TasaPromov11+Sector+Frcn_QuintilIVSHogares+Frcn_QuintilNoAsist4a17+Ambito+Duracion+Mat_Total+Sec_Total+AxS,
  thread.count = 8
)


######## /Stratified Regional Effects ##########

######## Result output ########



sink('~/r_results/result.JC.txt')
summary(result.JC[[2]])
print('AMsmallest.p.value')
print(result.JC[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JC[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JE.txt')
summary(result.JE[[2]])
print('AMsmallest.p.value')
print(result.JE[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JE[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.estatal.txt')
summary(result.estatal[[2]])
print('AMsmallest.p.value')
print(result.estatal[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.estatal[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.privado.txt')
summary(result.privado[[2]])
print('AMsmallest.p.value')
print(result.privado[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.privado[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.rural.txt')
summary(result.rural[[2]])
print('AMsmallest.p.value')
print(result.rural[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.rural[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.urbano.txt')
summary(result.urbano[[2]])
print('AMsmallest.p.value')
print(result.urbano[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.urbano[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.CBA.txt')
summary(result.CBA[[2]])
print('AMsmallest.p.value')
print(result.CBA[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.CBA[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.BA.txt')
summary(result.BA[[2]])
print('AMsmallest.p.value')
print(result.BA[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.BA[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Catamarca.txt')
summary(result.Catamarca[[2]])
print('AMsmallest.p.value')
print(result.Catamarca[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Catamarca[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Cordoba.txt')
summary(result.Cordoba[[2]])
print('AMsmallest.p.value')
print(result.Cordoba[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Cordoba[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Corrientes.txt')
summary(result.Corrientes[[2]])
print('AMsmallest.p.value')
print(result.Corrientes[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Corrientes[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Chaco.txt')
summary(result.Chaco[[2]])
print('AMsmallest.p.value')
print(result.Chaco[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Chaco[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Chubut.txt')
summary(result.Chubut[[2]])
print('AMsmallest.p.value')
print(result.Chubut[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Chubut[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Entre.txt')
summary(result.Entre[[2]])
print('AMsmallest.p.value')
print(result.Entre[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Entre[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Formosa.txt')
summary(result.Formosa[[2]])
print('AMsmallest.p.value')
print(result.Formosa[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Formosa[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Jujuy.txt')
summary(result.Jujuy[[2]])
print('AMsmallest.p.value')
print(result.Jujuy[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Jujuy[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.LaPampa.txt')
summary(result.LaPampa[[2]])
print('AMsmallest.p.value')
print(result.LaPampa[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.LaPampa[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.LaRioja.txt')
summary(result.LaRioja[[2]])
print('AMsmallest.p.value')
print(result.LaRioja[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.LaRioja[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Mendoza.txt')
summary(result.Mendoza[[2]])
print('AMsmallest.p.value')
print(result.Mendoza[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Mendoza[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Misiones.txt')
summary(result.Misiones[[2]])
print('AMsmallest.p.value')
print(result.Misiones[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Misiones[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Neuquen.txt')
summary(result.Neuquen[[2]])
print('AMsmallest.p.value')
print(result.Neuquen[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Neuquen[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Negro.txt')
summary(result.Negro[[2]])
print('AMsmallest.p.value')
print(result.Negro[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Negro[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Salta.txt')
summary(result.Salta[[2]])
print('AMsmallest.p.value')
print(result.Salta[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Salta[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.SanJuan.txt')
summary(result.SanJuan[[2]])
print('AMsmallest.p.value')
print(result.SanJuan[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.SanJuan[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.SanLuis.txt')
summary(result.SanLuis[[2]])
print('AMsmallest.p.value')
print(result.SanLuis[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.SanLuis[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.SantaCruz.txt')
summary(result.SantaCruz[[2]])
print('AMsmallest.p.value')
print(result.SantaCruz[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.SantaCruz[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.SantaFe.txt')
summary(result.SantaFe[[2]])
print('AMsmallest.p.value')
print(result.SantaFe[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.SantaFe[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Santiago.txt')
summary(result.Santiago[[2]])
print('AMsmallest.p.value')
print(result.Santiago[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Santiago[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Tucuman.txt')
summary(result.Tucuman[[2]])
print('AMsmallest.p.value')
print(result.Tucuman[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Tucuman[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.Tierra.txt')
summary(result.Tierra[[2]])
print('AMsmallest.p.value')
print(result.Tierra[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.Tierra[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JC.estatal.txt')
summary(result.JC.estatal[[2]])
print('AMsmallest.p.value')
print(result.JC.estatal[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JC.estatal[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JC.privado.txt')
summary(result.JC.privado[[2]])
print('AMsmallest.p.value')
print(result.JC.privado[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JC.privado[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JC.rural.txt')
summary(result.JC.rural[[2]])
print('AMsmallest.p.value')
print(result.JC.rural[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JC.rural[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JC.urbano.txt')
summary(result.JC.urbano[[2]])
print('AMsmallest.p.value')
print(result.JC.urbano[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JC.urbano[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JE.estatal.txt')
summary(result.JE.estatal[[2]])
print('AMsmallest.p.value')
print(result.JE.estatal[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JE.estatal[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JE.privado.txt')
summary(result.JE.privado[[2]])
print('AMsmallest.p.value')
print(result.JE.privado[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JE.privado[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JE.rural.txt')
summary(result.JE.rural[[2]])
print('AMsmallest.p.value')
print(result.JE.rural[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JE.rural[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JE.urbano.txt')
summary(result.JE.urbano[[2]])
print('AMsmallest.p.value')
print(result.JE.urbano[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JE.urbano[[3]]$AMsmallestVarName)
sink()



sink('~/r_results/result.JCJE.11.12.txt')
summary(result.JCJE.11.12[[2]])
print('AMsmallest.p.value')
print(result.JCJE.11.12[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JCJE.11.12[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JC.12.13.txt')
summary(result.JC.12.13[[2]])
print('AMsmallest.p.value')
print(result.JC.12.13[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JC.12.13[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JE.12.13.txt')
summary(result.JE.12.13[[2]])
print('AMsmallest.p.value')
print(result.JE.12.13[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JE.12.13[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JCJE.12.13.txt')
summary(result.JCJE.12.13[[2]])
print('AMsmallest.p.value')
print(result.JCJE.12.13[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JCJE.12.13[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JC.13.14.txt')
summary(result.JC.13.14[[2]])
print('AMsmallest.p.value')
print(result.JC.13.14[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JC.13.14[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JE.13.14.txt')
summary(result.JE.13.14[[2]])
print('AMsmallest.p.value')
print(result.JE.13.14[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JE.13.14[[3]]$AMsmallestVarName)
sink()


sink('~/r_results/result.JCJE.13.14.txt')
summary(result.JCJE.13.14[[2]])
print('AMsmallest.p.value')
print(result.JCJE.13.14[[3]]$AMsmallest.p.value)
print('AMsmallestVarName')
print(result.JCJE.13.14[[3]]$AMsmallestVarName)
sink()
