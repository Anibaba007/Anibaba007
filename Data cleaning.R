library(reshape2)
library(dplyr)
library(ggplot2)
library(vegan)
library(stringr)
library(rgbif)
library(tidyverse)
library(ape)
library(V.PhyloMaker)
library(phytools)
library(PhyloMeasures)
library(ggtree)
library(treeio)
library(picante)



veg.ini <- read.csv("veg.csv", header = TRUE, check.names = FALSE)

# melt to dataframe to long format
veg.ini_2 <- reshape2::melt(veg.ini, id.vars ="species", na.rm = FALSE)

write.csv(veg.ini_2,"species1.csv") # save as cvs file

# GBIF authentication of species names

ns<-read.csv('gbif.csv',sep=';')

sp.ex<-name_backbone(name='Solidago gigantea') # view the complete of the species 
View(sp.ex)

# run the gbif command 
check_gbif_fun<-function(x) name_backbone(name = x, kingdom = 'plants')

gbifkey <- lapply(ns$species, check_gbif_fun)
for(i in 1:length(gbifkey))gbifkey[[i]]<-as.data.frame(gbifkey[[i]])
gbifkey <- bind_rows(gbifkey)
gbifkey$orig<-ns$species
write.table(gbifkey,'gbif_names.csv', sep = ',', row.names = TRUE)


# sum and group species in layer

species <- read.csv("species_data.csv")

species_merged <- species%>% group_by(plotid,species)%>%summarise(cover)

write.csv(species_merged,"species_sum.csv")

## Use dcast function to transpose the data into wide format 

species_sum <-read.csv("species_sum.csv")

species_trans <- dcast(species_sum, plotid~ species, 
      value.var="cover", fill=0,
      fun.aggregate = sum)
write.csv(species_trans,"community.csv")

species_layer <- read.csv("species_data.csv")

species_trans <- dcast(species_layer, plotid + layer~ species, 
                       value.var="cover", fill=0,
                       fun.aggregate = sum)

write.csv(species_layer,"community_layer.csv")


## generate a phylogeny tree for community data

comm.phylo <- read.csv("community_phylo_revised.csv")

## Create phylogenetic tree 

comm.phylo1 = phylo.maker(comm.phylo)

comm.tree = comm.phylo1$scenario.3

# Plot the tree
plotTree(comm.tree)
windows()
ggtree(comm.tree, branch.length='none', layout='circular') + geom_tiplab()


# export the tree for manual editing or use in other software 
write.tree(comm.tree, "commmuity_revised.tre")


# Import the tree back into R
# Or import a tree made with other software into R

tree <- read.tree("tree.tre")

ggplot(tree, aes(y)) + geom_tree() + theme_tree()

# To show all it internal node/tips 

ggtree(comm.tree) + geom_nodepoint(color="#b5e521", alpha=1/4, size=10) + geom_tippoint(color="#FDAC4F", shape=8, size=3)


# Displaying labels

p + geom_tiplab(size=2, color="purple")

ggtree(comm.tree) + geom_tiplab(as_ylab = FALSE, color = "firebrick")


ggtree(comm.tree, layout = "circular") + geom_tiplab(aes(angle=angle),size=2, color="blue")


#traits imputation
#in this file we get a dataset with species traits and we will conduct phylogenetic-based imputaion

#load required libraries
require(tidyverse)# for data processing and ggplots
require(V.PhyloMaker)# source of a phylogenetic megatree 
require(PVR)#for calculation of phylogenetic distance matrix eigenvalues
require(missForest)# for random forest-based imputation
require(reshape2)#for melting data.frames
require(ggthemes)#for nice plots layout


traits<-read.csv('traits.csv') #.csv file with traits and taxonomic information (species/genus/family)
summary(traits)
#traits are columns 10-25

?read.csv
#phylogenetic tree-----
#using  V.phylomaker package we need to extract tree from a megatree, we are going to specify as an argument dataframe with three columns - family, genera and species. We use scenario.3 as default, recommended scenario in the V.phylomaker package
phylo<-phylo.maker(traits[,1:3])$scenario.3

plot(phylo, type='fan', cex=.8)
#quick check of cladograme

#phylo.imputations----------
# following Penone et al. (2014) Meth Ecol Evol 
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12232

deco<-PVRdecomp(phylo) #let's decompose the matrix of phylogenetic distances into eigenvalues
cumsum(deco@Eigen$values[1:100]/sum(deco@Eigen$values))#we need to know how much of phylogenetic distance variability is explained by particular eigenvectors, Penone et al. used 15 which explained 60% of variability, here we can see that 16 first eigenvectors explain 60.16% of variability

#as some of variables were presented as character, we need to code them as factors, i.e. Ellenberg's EIVs, pollination or flowering months
for(i in c(14:16)) traits[,i]<-factor(traits[,i], ordered = T)


set.seed(111) #to maintain replicability let's set seed of pseudorandom numbers generator
base.miss<-missForest(cbind(traits[,c(4:21)],deco@Eigen$vectors[,1:16]), ntree=500, maxiter = 15)
traits$liform<-factor(traits$liform)
traits$Alien.status<-factor(traits$Alien.status)
traits$Invasive.status<-factor(traits$Invasive.status)
traits$Residence.time.status<-factor(traits$Residence.time.status)
traits$EIV.T<-as.numeric(traits$EIV.T)
traits$EIV.M<-as.numeric(traits$EIV.M)
traits$EIV.SR<-as.numeric(traits$EIV.SR)
traits$EIV.N<-as.numeric(traits$EIV.N)

#this may keep some time
#we bind our traits matrix with eignecectors, and then every single column containing missing values is predicted based on remaining columns. This literally relies on phylogenetic similarity of species and known traits values
#let's see how good is the fitness
base.miss$OOBerror
#NRMSE is not bad - this is standardized RMSE for numeric traits
#PFC is fraction of false classified cases - here 40.5% is relatively low value, assuming multilevel traits, e.g. Ellenberg's indicators
# in case of prevalence of known trait values, only not-known are imputed


###comparison for appendix pictures----
#we need to bind original and predicted results
comparison<-rbind(traits[,c(4:21)],base.miss$ximp[,1:18])
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
for(i in c(1:3,11:14)) comparison[,i]<-as.numeric(as.factor(comparison[,i]))

app.imp.cont<-comparison[c(1:18,19)]%>%melt(id.vars='type',variable.name='trait',value.name='val')
#here we melted wide table into long table, to make easy facet-based ggplot
ggplot(app.imp.cont, aes(x=val,col=type))+facet_wrap(~trait,scales='free_x',ncol=4)+geom_density(aes(y=..scaled..))+theme_few()+theme(legend.position = 'bottom')+scale_color_manual(values=c('red','blue'))
#here we can see the overlap between predicted and real values and low errors in EIV-N and LDMC
#
ggplot(app.imp.cont, aes(x=val,col=type))+facet_wrap(~trait,scales='free_x',ncol=4)+geom_density(aes(y=..scaled..))+theme_few()+theme(legend.position = 'bottom')+scale_color_manual(values=c('red','blue'))+scale_x_log10()
#some variables are better to be shown with log scales

#and for factors:
#we melt data.frame with traits twice - first with number of all species with particular trait level, excluding NAs, second - with counts of all species with known information about trait level. This will show us how proportions of traits levels vary between original and imputed datasets
app.imp.fac<-comparison[,c(1:3,9:14,19)]%>%melt(id.vars='type',variable.name='trait',value.name='val')%>%filter(!is.na(val))%>%group_by(type,trait,val)%>%summarise(count=n())
#here we first melted data.frame, then we filter out NAs, groupbed by variables and summarized
app.imp.fcount<-comparison[,c(1:3,9:14,19)]%>%melt(id.vars='type',variable.name='trait',value.name='val')%>%filter(!is.na(val))%>%group_by(type,trait)%>%summarise(counttrait=n())
app.imp.fac<-left_join(app.imp.fac,app.imp.fcount)
#let's join these data.frames
app.imp.fac$perc<-app.imp.fac$count/app.imp.fac$counttrait
#And calculate the percentage fraction (perc)
ggplot(app.imp.fac,aes(y=perc,x=type,fill=val))+geom_col(position='stack',col='black')+facet_wrap(~trait, ncol=4)+theme_few()+theme(legend.position = 'bottom')

#save results
#joint columns with no traits (1:8) from original database and imputed traits dataset (i.e. original values + NAs replaced by estimation)
write.table(cbind(traits[,1:3],base.miss$ximp[,1:18]),'traits_impute.csv', quote = FALSE, sep = ",")




