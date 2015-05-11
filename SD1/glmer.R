rm(list=ls())
require(lme4)
require(afex)
require(ggplot2)
setwd("~/Documents/SD1/Sample_info/")
dataglm=read.csv("Transmission_data.csv", header=T)
dataglm<-as.data.frame(dataglm)

#Diseased vs healthy (binomial)
mixed(Disease_state.1 ~ Dose + Site_dose + Site_collected + Site_dose:Site_collected + Dose:Site_dose + Tank + (1|Genotype),
      data = dataglm, family=binomial, type=3, method="LRT")

#time_survived.
mixed(Time_survived ~ Dose + Site_dose + Site_collected + Dose:Site_dose +Site_dose:Site_collected + (1|Genotype) + (1|Tank), 
      data = dataglm, type=3, method="LRT")


Diseased<-subset(dataglm, Dose=="Disease" & Disease_state=="Diseased")
mixed(Time_diseased ~ Site_dose + Site_collected + Site_dose:Site_collected+ (1|Genotype) + (1|Tank), 
      data = Diseased, type=3, method="LRT")

# group by Treatment and Family, calculate mean abundance and standard deviation
avgs <- ddply(dataglm, c("Site_dose", "Dose"), summarise, mean=mean(Time_survived))

# plot bar graph with standard deviation as error bars
means.barplot=qplot(y=mean, x=Dose, fill=Site_dose,
                    data=avgs, geom="bar", stat="identity",
                    position="dodge")

means.barplot

means.sem <- ddply(dataglm, c("Dose", "Site_dose"), summarise,
                   mean=mean(Time_survived), sem=sd(Time_survived)/sqrt(length(Time_survived)))
means.sem <- transform(means.sem, lower=mean-sem, upper=mean+sem)

means.barplot + geom_errorbar(aes(ymax=upper,
                                  ymin=lower),
                              position=position_dodge(0.9),
                              data=means.sem)





#try taking disease out
Diseased<-subset(dataglm, Dose=="Disease")
mixed(Disease_state.1 ~ Site_dose + Site_collected + (1|Genotype) + (1|Tank/Sample), 
      data = Diseased, family=binomial, type=3, method="LRT")


summary(mod.lme4)
aov4(mod.lme4)
rand(mod.lme4)