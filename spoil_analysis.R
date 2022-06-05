
# DCA analysis on spoil heaps data

# 1. Test how site differ with spoil heap claases and successional age

library(vegan)
library(tidyverse)
library(ggvegan)

comm_spoil <- read.csv("community_spoilheaps_0_revised_04_2022.csv", row.names = 1)

meta_spoil <- read.csv("community_meta_spoilheaps_0_revised_04_2022.csv")

str(meta_spoil$age)
# Verify the norm of the transformed row vectors
# Write a 1-line function that computes the norm of vector x
vec.norm <- function(x) sqrt(sum(x ^ 2))

## DCA on the log trannsformed

DCA <- decostand (comm_spoil, "hell") # log tranformed
DCA
DCA <- decorana(DCA)
DCA

## DCA on the log trannsformed and hellinger standardized data


# Since DCA does not show total initial value.
# The total inertia can be gotten by applying correspondence analysis (CA) on your data:

cca(log1p (comm_spoil))

# total inertia is 31.54, and if needed, variation captured by particular axes can be calculated as eigenvalue/total inertia (e.g., for the first axis, 1.00/31.54*100 = 3.17%


## Plot 

ordiplot (DCA, display = 'species', type = 'n') # Display = site; type = none
points (DCA, col = meta_spoil$class, pch = meta_spoil$class )

#plot empty space
plot(DCA, type='n')
#labels, shortened using make.cepnames function
orditorp(DCA, disp='spe',labels = make.cepnames(colnames(comm_spoil)))


# See symp.albu on right side. Which point is it?
#plot empty space
plot(DCA, type='n')
#labels, display=site names
orditorp(DCA, disp='sites',labels = rownames(comm_spoil))


# something _3_1_E, we can see site scores

View(scores(DCA)) # and sort by decreasing order:


# It is plot 20_3_1_E.

# See all columns values for that row:
  
comm_spoil[rownames(comm_spoil)=='20_3_1_E',]


# If I exclude that observation we have a better image:

comm_spoil_1<-comm_spoil[-which(rownames(comm_spoil)=='20_3_1_E'),]

meta_spoil_1<-meta_spoil[meta_spoil$plotid!='20_3_1_E',]


dca_hel <- decostand(comm_spoil_1, "hell")
dca_hel<-decorana(dca_hel)

plot(dca_hel)

orditorp(dca_hel, disp='spe',labels = make.cepnames(colnames(comm_spoil_1)))
points (dca_hel, col = meta_spoil_1$class, pch = meta_spoil_1$class )

# NMDS

NMDS <- metaMDS (comm_spoil_1)

NMDS
# Plot NMDS

ordiplot (NMDS, type = 't')
stressplot(NMDS)
plot (NMDS, display = 'sites', type = 't', main = 'Goodness of fit') # this function draws NMDS ordination diagram with sites
points (NMDS, display = 'sites', cex = goodness (NMDS)*200) # and this adds the points with size reflecting goodness of fit (bigger = worse fit)



vitava_fig <- ordiplot(dca_hel, display = "sites", type = "none")

points(vitava_fig,"sites", pch = 19, col = "green", select = meta_spoil$class == "1") # points for group 1

points(vitava_fig,"sites", pch = 20, col = "blue", select = meta_spoil$class == "2") # points for group 2

points(vitava_fig,"sites", pch = 25, col = "purple", select = meta_spoil$class == "3") # points for group 3


points(vitava_fig,"sites", pch = 23, col = "black", select = meta_spoil$class == "4") # points for group 4

points(vitava_fig,"sites", pch = 24, col = "pink", select = meta_spoil$class == "5") # points for group 4


plot(envfit(dca_hel, meta_spoil$age==("1")), text(meta_spoil, labs="1")) ## age 1
plot(envfit(dca_hel, meta_spoil$age== as.factor("2"))) ## age 2
plot(envfit(dca_hel, meta_spoil$age== as.factor("3"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("4"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("5"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("6"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("7"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("8"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("10"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("11"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("12"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("13"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("14"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("25"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("28"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("30"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("34"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("38"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("40"))) 
plot(envfit(dca_hel, meta_spoil$age== as.factor("56"))) 



dca_fit <-envfit(dca_hel~age,meta_spoil_1, perm = 999)

ordispider(dca_hel, meta_spoil_1$age, col="skyblue")
points(dca_hel, display = "sites", col = as.numeric(meta_spoil_1$age), pch=16)
plot(dca_fit, cex=1.2, axis=TRUE)



labels(dca_fit)
plot(dca_hel)
plot(dca_fit, labels=list(factors = paste("A", c(1,2,3,4))),
     bg = rgb(1,1,0,0.5))


??bg

?ordiarrows
?ordipointlabel 
?plot.decorana



