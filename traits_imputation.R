#traits imputation
#in this file we get a dataset with species traits and we will conduct phylogenetic-based imputaion

#load required libraries
require(tidyverse)# for data processing and ggplots
require(V.PhyloMaker)# source of a phylogenetic megatree 
require(PVR)#for calculation of phylogenetic distance matrix eigenvalues
require(missForest)# for random forest-based imputation
require(reshape2)#for melting data.frames
require(ggthemes)#for nice plots layout


traits<-read.csv('traits_before.csv',sep=';') #.csv file with traits and taxonomic information (species/genus/family)
summary(traits)
#traits are columns 10-25

#phylogenetic tree-----
#using  V.phylomaker package we need to extract tree from a megatree, we are going to specify as an argument dataframe with three columns - family, genera and species. We use scenario.3 as default, recommended scenario in the V.phylomaker package
phylo<-phylo.maker(traits[,1:3])$scenario.3

plot(phylo, type='fan', cex=.6)
#quick check of cladograme

#phylo.imputations----------
# following Penone et al. (2014) Meth Ecol Evol 
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12232

deco<-PVRdecomp(phylo) #let's decompone the matrix of phylogenetic distances into eigenvalues
cumsum(deco@Eigen$values[1:100]/sum(deco@Eigen$values))#we need to kno how much of phylogenetic distance variability is explained by particular eigenvectors, Penone et al. used 15 which explained 60% of variability, here we can see that 16 first eigenvectors explain 60.16% of variability

#as some of variables were presented as factors, we need to code them as factors, i.e. Ellenberg's EIVs, pollination or flowering months
for(i in c(10:19)) traits[,i]<-factor(traits[,i], ordered = T)

set.seed(111) #to maintain replicability let's set seed of pseudorandom numbers generator
baza.miss<-missForest(cbind(traits[,c(10:25)],deco@Eigen$vectors[,1:16]), ntree=500, maxiter = 15)
#this may keep some time
#we bind our traits matrix with eignecectors, and then every single column containing missing values is predicted based on remaining columns. This literally relies on phylogenetic similarity of species and known traits values
#let's see how good is the fitness
baza.miss$OOBerror
#NRMSE is not bad - this is standardized RMSE for numeric traits
#PFC is fraction of falsce classified cases - here 40.5% is relatively low value, assuming multilevel traits, e.g. Ellenberg's indicators
# in case of prevalence of known trait values, only not-known are imputed


###comparison for appendix pictures----
#we need to bind original and predicted results
comparison<-rbind(traits[,c(10:25)],baza.miss$ximp[,1:16])
#only columns with traits
comparison$type<-c(rep(c('orig','imputed'),each=nrow(traits)))
#new column divideds original and imputed values

#we can compare density of distributions
ggplot(comparison, aes(x=SLA,col=type))+geom_density()
ggplot(comparison, aes(x=H,col=type))+geom_density()+scale_x_log10()

#even within the groups
ggplot(comparison, aes(x=liform,y=LDMC,col=type))+geom_boxplot()#+scale_y_log10()

#separately we have to show the continuous traits
#firstly, let's change factors into numeric backward
for(i in c(1:7)) comparison[,i]<-as.numeric(as.character(comparison[,i]))

app.imp.cont<-comparison[c(1:7,12:15,17)]%>%melt(id.vars='type',variable.name='trait',value.name='val')
#here we melted wide table into long table, to make easy facet-based ggplot
ggplot(app.imp.cont, aes(x=val,col=type))+facet_wrap(~trait,scales='free_x',ncol=4)+geom_density(aes(y=..scaled..))+theme_few()+theme(legend.position = 'bottom')+scale_color_manual(values=c('red','blue'))
#here we can see the overlap between predicted and real values and low errors in EIV-N and LDMC
#
ggplot(app.imp.cont, aes(x=val,col=type))+facet_wrap(~trait,scales='free_x',ncol=4)+geom_density(aes(y=..scaled..))+theme_few()+theme(legend.position = 'bottom')+scale_color_manual(values=c('red','blue'))+scale_x_log10()
#some variables are better to be shown with log scales

#and for factors:
#we melt data.frame with traits twice - first with number of all species with particular trait level, excluding NAs, second - with counts of all species with known information about trait level. This will show us how proportions of traits levels vary between original and imputed datasets
app.imp.fac<-comparison[,c(8:11,16,17)]%>%melt(id.vars='type',variable.name='trait',value.name='val')%>%filter(!is.na(val))%>%group_by(type,trait,val)%>%summarise(count=n())
#here we first melted data.frame, then we filter out NAs, groupbed by variables and summarized
app.imp.fcount<-comparison[,c(8:11,16,17)]%>%melt(id.vars='type',variable.name='trait',value.name='val')%>%filter(!is.na(val))%>%group_by(type,trait)%>%summarise(counttrait=n())
app.imp.fac<-left_join(app.imp.fac,app.imp.fcount)
#let's join these data.frames
app.imp.fac$perc<-app.imp.fac$count/app.imp.fac$counttrait
#And calculate the percentage fraction (perc)
ggplot(app.imp.fac,aes(y=perc,x=type,fill=val))+geom_col(position='stack',col='black')+facet_wrap(~trait, ncol=4)+theme_few()+scale_fill_manual(values=c('white','black','#1a9850','#bf5b17','#a6d96a','#66c2a5','#fee08b','#fdc086','#7fc97f','#386cb0','#beaed4','#f0027f','#666666','#2166ac','#d73027','#c51b7d','#ffff99'))+theme(legend.position = 'bottom')

#save results
#joint columns with no traits (1:8) from original database and imputed traits dataset (i.e. original values + NAs replaced by estimation)
write.table(cbind(traits[,1:9],baza.miss$ximp[,1:16]),'traits_new.csv',sep=';')
