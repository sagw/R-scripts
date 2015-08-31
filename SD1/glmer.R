rm(list=ls())
require(lme4)
require(afex)
require(ggplot2)
require(plyr)
require(lsmeans)
require(multcompView)
setwd("~/Documents/SD1/Sample_info/")
dataglm=read.csv("Transmission_data.csv", header=T)
dataglm=as.data.frame(dataglm)

#Diseased vs healthy (binomial)
mixed=mixed(Disease_state.1 ~  Dose*Site_dose*Site_collected+
    (1|Tank)+(1|Genotype), data=dataglm, method="LRT", type=3)

lsmeans(mixeds,  cld~  Dose*Site_dose)
lsmip(mixeds,~  Dose*Site_dose)

hsd=HSD.test(mixed, trt= "Dose:Site_dose:Site_collected", group="TRUE")
mixed=mixed(Disease_state.1 ~  Dose+Site_dose+ Dose:Site_dose+
              (1|Tank)+(1|Genotype), data=dataglm, method="LRT", type=3)


#time_survived
mixeds=mixed(Time_survived ~ Dose*Site_dose*Site_collected + (1|Tank) + (1|Genotype), 
      data = dataglm, type=3, method="LRT")

Diseased<-subset(dataglm, Dose=="Disease" & Disease_state=="Diseased")
Diseased<-subset(dataglm, Dose=="Disease")
mixed(Time_diseased ~ Site_dose*Site_collected+ (1|Genotype) + (1|Tank), 
      data = Diseased)

# group by Treatment and Family, calculate mean abundance and standard deviation
avgs = ddply(Diseased, c("Site_dose", "Disease_state"), summarise, mean=mean(Time_survived))

# plot bar graph with standard deviation as error bars
means.barplot=qplot(y=mean, x=Disease_state, fill=Site_dose,
                    data=avgs, geom="bar", stat="identity",
                    position="dodge")

means.barplot

means.sem <- ddply(Diseased, c("Disease_state", "Site_dose"), summarise,
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


endoglm=merge(b1tendo, dataglm, by.x="row.names", by.y="Sample")

endoglm$count=as.numeric(endoglm$count)

endod=subset(endoglm, Dose=="Disease")
endolm=lm(Time_survived ~ count*Tank , endod)
dotplot(Time_survived~Genotype, Diseased, las=2, cex.axis=0.5)