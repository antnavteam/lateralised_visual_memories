library("beeswarm")

dati<-read.csv("table_mirror_means.csv",header=T,sep=",")
head(dati)
table(dati$cond)


# PWR analysis
#-------------PWR is between 0 (weak) and 1 (strong). 0.8 is excellent
library(pwr)




data1=(subset(dati,cond=='Trained_side_Realease_route')$before)
data2=(subset(dati,cond=='Trained_side_Realease_route')$mirror)

n1= length(data1)
n2 = length(data2)

var1 = var(data1)
var2 = var(data2)
pooled_std = sqrt( (var1+var2)/2) #var = square of std

mean1= mean(data1)
mean2= mean(data2)

effect_size= (mean1-mean2)/pooled_std
d = effect_size

pwr.t2n.test(n1,n2,d,alternative='two.sided')
pwr.t2n.test(n1,n2,d,alternative='greater')
pwr.t2n.test(n1,n2,d,alternative='less')


# n1 = 6
# n2 = 6
# d = 3.901028
# sig.level = 0.05
# power = 0.9999757
# alternative = two.sided



data1=(subset(dati,cond=='Trained_side_Realease_unfam')$before)
data2=(subset(dati,cond=='Trained_side_Realease_unfam')$mirror)

n1= length(data1)
n2 = length(data2)

var1 = var(data1)
var2 = var(data2)
pooled_std = sqrt( (var1+var2)/2) #var = square of std

mean1= mean(data1)
mean2= mean(data2)

effect_size= (mean1-mean2)/pooled_std
d = effect_size

pwr.t2n.test(n1,n2,d,alternative='two.sided')
pwr.t2n.test(n1,n2,d,alternative='greater')
pwr.t2n.test(n1,n2,d,alternative='less')

# n1 = 12
# n2 = 12
# d = 0.7817888
# sig.level = 0.05
# power = 0.4488735
# alternative = two.sided



data1=(subset(dati,cond=='Trained_straight_Realease_route')$before)
data2=(subset(dati,cond=='Trained_straight_Realease_route')$mirror)

n1= length(data1)
n2 = length(data2)

var1 = var(data1)
var2 = var(data2)
pooled_std = sqrt( (var1+var2)/2) #var = square of std

mean1= mean(data1)
mean2= mean(data2)

effect_size= (mean1-mean2)/pooled_std
d = effect_size

pwr.t2n.test(n1,n2,d,alternative='two.sided')
pwr.t2n.test(n1,n2,d,alternative='greater')
pwr.t2n.test(n1,n2,d,alternative='less')

# n1 = 12
# n2 = 12
# d = 1.216573
# power = 0.8126733
# sig.level = 0.05
# alternative = two.sided


#-------------------------------------One sample power analysis for natural sun above 0
library(pwr)

data1=(subset(dati,cond=='Trained_side_Realease_route')$before)
mu0 <- 0  # null hypothesis mean

mean1= mean(data1)
sd1 = sd(data1)
effect_size= (mean1-mu0)/sd1
d = effect_size

pwr.t.test(d = d, n = length(data1), sig.level = 0.05, type = "one.sample", alternative = "less")

# n = 6
# d = -2.20397
# sig.level = 0.05
# power = 0.998115
# alternative = less


data1=(subset(dati,cond=='Trained_side_Realease_unfam')$before)
mu0 <- 0  # null hypothesis mean

mean1= mean(data1)
sd1 = sd(data1)
effect_size= (mean1-mu0)/sd1
d = effect_size

pwr.t.test(d = d, n = length(data1), sig.level = 0.05, type = "one.sample", alternative = "less")

# n = 12
# d = 0.6164422
# sig.level = 0.05
# power = 0.0001262186
# alternative = less



data1=(subset(dati,cond=='Trained_straight_Realease_route')$before)
mu0 <- 0  # null hypothesis mean

mean1= mean(data1)
sd1 = sd(data1)
effect_size= (mean1-mu0)/sd1
d = effect_size

pwr.t.test(d = d, n = length(data1), sig.level = 0.05, type = "one.sample", alternative = "less")

# n = 12
# d = 0.7769941
# sig.level = 0.05
# power = 1.415139e-05
# alternative = less


##--------------------- Test if slope is more consistent for trainedA

dati$vartotest = dati$before - dati$mirror
beeswarm(vartotest ~ cond, data=dati)

library(nlme)
model<-lm(vartotest ~ cond, data=dati)
summary(model)


#-----------------Test for strength of the mirror effect independantly of slope direction (absolute slope) if stronger in Trained A 
dati$vartotest = abs(dati$before - dati$mirror)
beeswarm(vartotest ~ cond, data=dati)
boxplot(vartotest ~ cond, data=dati)

library(nlme)
model<-lm(vartotest ~ cond, data=dati)
summary(model)


plot(mean_mirror ~ mean_natural, data=subset(dati,group=='untrained'))
plot(mean_mirror ~ mean_natural, data=dati)

library(ggplot2)
ggplot(dati, aes(x=mean_mirror, y=mean_natural, color=group)) + geom_point()


