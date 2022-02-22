df = read.delim('/home/arnie/Documents/IITM/sem2/Dispersed Multiphaseflow(ME5310)/assignment 2/Particle size distribution.txt')
summary(df)
attach(df)
library(ggplot2)
library(stringr)
library(base)

# que 1 histogram of particle distribution
ggplot(df,aes(Particle.size.in.micron))+geom_histogram(bins=2,colour = 'black',fill=I("darkgrey"))+ggtitle("Number of bins 2")+labs(x="Particle size in micron",y="Frequency")+geom_freqpoly(stat = "bin",position = "stack")+theme_bw()
ggplot(df,aes(Particle.size.in.micron))+geom_histogram(bins=10,colour = 'black',fill=I("darkgrey"))+ggtitle("Number of bins 10")+labs(x="Particle size in micron",y="Frequency")+geom_freqpoly(stat = "bin",position = "stack")+theme_bw()
ggplot(df,aes(Particle.size.in.micron))+geom_histogram(bins=20,colour = 'black',fill=I("darkgrey"))+ggtitle("Number of bins 20")+labs(x="Particle size in micron",y="Frequency")+geom_freqpoly(stat = "bin",position = "stack")+theme_bw()
ggplot(df,aes(Particle.size.in.micron))+geom_histogram(bins=30,colour = 'black',fill=I("darkgrey"))+ggtitle("Number of bins 30")+labs(x="Particle size in micron",y="Frequency")+geom_freqpoly(stat = "bin",position = "stack")
  
# que 2 number probability
ggplot(df,aes(Particle.size.in.micron),after_stat(density))+geom_histogram(bins=30,colour="white")+ggtitle("Number of bins =30")+theme_bw()
b2=table(cut(x=df$Particle.size.in.micron,breaks = 2))
b10=table(cut(x=df$Particle.size.in.micron,breaks = 10))
b20=table(cut(x=df$Particle.size.in.micron,breaks = 20))
b30=table(cut(x=df$Particle.size.in.micron,breaks = 30))
total = length(Particle.size.in.micron)



probDense = function(table,total){
  prdist = c()
  df2=data.frame(table)
  for(i in 1:length(table))
    prdist[i]=table[i]/total
  
  # avg density
  df2$Var1 = gsub("\\(","",as.character(df2$Var1))
  df2$Var1 = gsub("\\]","",as.character(df2$Var1))
  df2[c("ldia",'udia')]=str_split_fixed(df2$Var1,",",2)
  df2$ldia = as.numeric(df2$ldia)
  df2$udia = as.numeric(df2$udia)
  df2$di = (df2$ldia+df2$udia)/2
  
  # number probability
  df2$numberprob=prdist
  df2$numberprob_cum = cumsum(df2$numberprob)
  
  # mass probability
  df2$massprob = df2$Freq*(df2$di)^3
  totalmass = sum(df2$massprob)
  print(paste("Total mass ",totalmass))
  df2$massprob = df2$massprob/totalmass
  df2$massprob_cum = cumsum(df2$massprob)
  df2$Freq_cum = cumsum(df2$Freq)
  
  # mean and mode number probability
  median =0
  mode = 0
  for(i in 1:length(df2$Freq)){
    if(df2$Freq_cum[i]>=sum(df2$Freq)/2){
      print(paste(" i == ",i))
      median = df2$ldi[i] + (sum(df2$Freq)/2-df2$Freq_cum[i-1])/df2$Freq[i]*(df2$udi[i]-df2$ldi[i])
      mode =  df2$ldi[i] +(df2$Freq[i]-df2$Freq[i-1])/(2*df2$Freq[i]-df2$Freq[i+1]-df2$Freq[i-1])*(df2$udi[i]-df2$ldi[i])
      break
    }
  }
  print(paste("Median :",median))
  print(paste("Mode :",mode))
  
  # Arithmatic mean Diameter, Surface mean dia, Volume mean dia
  sigma_ni = sum(df2$Freq)
  sigma_di2_ni = sum(df2$di^2*df2$Freq)
  sigma_di3_ni = sum(df2$di^3*df2$Freq)
  print(paste("AMD",sum(df2$Freq * df2$di)/sigma_ni))
  print(paste("SMD : ",(sigma_di2_ni/sigma_ni)**(0.5)))
  print(paste("VMD : ",(sigma_di3_ni/sigma_ni)**(1/3)))
  
  # variance
  AMD = sum(df2$Freq * df2$di)/sigma_ni
  sigma_xbar_xi = sum((df2$di-AMD)**2)
  print(paste("Variance: ",sigma_xbar_xi/sigma_ni))
  # write.csv(df2,"/home/arnie/Documents/IITM/sem2/Dispersed Multiphaseflow(ME5310)/assignment 2/bin10.csv")
}
probDense(b20,total)

df3 = read.csv("/home/arnie/Documents/IITM/sem2/Dispersed Multiphaseflow(ME5310)/assignment 2/bin20.csv")

df4 = subset(df3,select=c(di,numberprob_cum,massprob_cum))
df4$rosinrammler_num = log(-log(1-df4$numberprob_cum))
df4$rosinrammler_mass = log(-log(1-df4$massprob_cum))
df4$ln_di = log(df4$di)
df4 = df4[-(length(df4$di)),]
ggplot(df4,aes(x=ln_di,y=rosinrammler_num))+geom_point(col="black",size=2.5,pch=19,show.legend = TRUE)+geom_smooth(method = lm,se=FALSE,size=1,col='red',linetype='dashed',show.legend = TRUE,formula = y~x)+ggtitle("Rossin-Rammler Distribution for Cummulative Number probability")+labs(x="ln(D)",y="ln(-ln(1-F(D)))") 
regline = lm(df4$rosinrammler_num~df4$ln_di,df4)
# y = mx + c
m = coef(regline)[2]
c = coef(regline)[1]


# shape and Scale parameter
print(paste("Shape parameter = ",m))

dnm = 0
for(i in 1:length(df3$Freq)){
  if(df3$numberprob_cum[i]>0.5){
    dnm = df3$di[i-1] + (0.5-df3$numberprob_cum[i-1])/(df3$numberprob_cum[i]-df3$numberprob_cum[i-1])*(df3$di[i]-df3$di[i-1])
    break
  }
}
print(paste("dnm = ",dnm))
delta = dnm/(log(2))^(1/m)
print(paste("Scale Parameter= ",dnm/(log(2))^(1/m)))

ggplot(df4,aes(x=ln_di,y=rosinrammler_mass))+geom_point(col="black",size=2.5,pch=19,show.legend = TRUE)+geom_smooth(method = lm,se=FALSE,size=1,col='red',linetype='dashed',show.legend = TRUE,formula = y~x)+ggtitle("Rossin-Rammler Distribution for Cummulative Mass Probability")+labs(x="ln(D)",y="ln(-ln(1-F(D)))") 
regline = lm(df4$rosinrammler_mass~df4$ln_di,df4)
# y = mx + c
m = coef(regline)[2]
c = coef(regline)[1]

dmm = 0
for(i in 1:length(df3$Freq)){
  if(df3$massprob_cum[i]>0.5){
    dmm = df3$di[i-1] + (0.5-df3$massprob_cum[i-1])/(df3$massprob_cum[i]-df3$massprob_cum[i-1])*(df3$di[i]-df3$di[i-1])
    break
  }
}
print(paste("dmm = ",dmm))
print(paste("Shape parameter = ",m))
print(paste("Scale Parameter = ",dmm/(log(2))^(1/m)))

