### install dependencies
pkg <- c('vegan','labdsv','ade4','ecodist','fso',
         'vegclust','ape','picante','mgcv','learnr')
has <- pkg %in% rownames(installed.packages())
if(any(!has)) install.packages(pkg[!has])

### specify path to wherever you saved the file (change as needed)
path_to_file <- 'C:/Users/user/Desktop/Plant Comm Eco/esa2021_0.0.1-20210722.tar.gz'
### install the main package
install.packages(path_to_file, repos=NULL, type='source')


require('esa2021')
learnr::run_tutorial('esa_tutorial', package='esa2021')

# Load package
require(ade4)
require(vegan)
require(labdsv)
require(ecodist)
require(vegclust)
require(ape)
require(picante)
require(mgcv)
require(fso)
require(ggplot2)
require(ggpubr)
# Useful functions to explore 
ls()
dim(spe)
NROW(spe)
NCOL(spe)
str(spe)
head(spe)
names(spe)
x <- spe[,11]
x
x <- spe[,'bolboschoenus_maritimus']
hist(x)
plot(x)
rm(x)

# Define objects and functions 

# function to pick an item from the vector by position
`item_picker` <- function(x, pos = NULL) {
  x[pos]
}
z <- 1:19               # assign a sequence of numbers to object `z`
item_picker(z, pos = 5) # pick the fifth element of `z`

# Load data

load('./Data/veg.rda') #Assumed to be in the working directory

xy <- veg$xy #Spatial data
spe <- veg$spe # Species data
env <- veg$env # environment data
tra <- veg$tra # traits data
phy <- veg$phy # phylogeny data
rm(veg) #cleanup
ls() #objects now in this local environment 

# Pre-Analysis

# Check Data Structure 
str(spe) # community (species) abundance matrix
str(env) # environment matrix
str(tra) # trait matrix
str(phy,1) # phylogeny
str(xy) # spatial coordinates


# Find and replace troublesome values 

### create example of matrix with some pesky values 

s <- spe
s[5,2]<- NA
s[6,4] <- (-999)

### Index the row and column of any NA values 
which(is.na(s), arr.ind = TRUE)

### Index the row and column of any negative values 
which(s < 0, arr.ind = TRUE)

### Replace value by logical test
s[is.na(s)] <- 777

head(s)

# DATA TRANSFORMATION

### load data

spe <- veg$spe # Species data
env <- veg$env # environment data
tra <- veg$tra # traits data

### basic transformations 

spet <- data.frame(log10(spe + 1))
envt <- data.frame(vegan::decostand(scale(env,center = F), 'range'))
trat <- data.frame(vegan::decostand(tra, 'range'))

### compare a few 

# Check if there is a linear relationship between actual data and transformed data
par(mfrow=c(1,1))

plot(tra$bfp,trat$bfp)
plot(env$k2o, envt$k2o)

### Abundances to presence/absence

s <- (spe>0) * 1 # from numeric to 0/1
s[1:5, 1:7] # peek at a few

range(s)

### Outliers 

### define outlier function

'outliers' <- function(x, mult=2, method='bray') {
  d <- as.matrix(vegan::vegdist(x, method = method, binary = F, diag = T, upper = T))
  diag(d) <- as.numeric(1)  # avoid zero-multiplication
  m <- apply(d, 2, mean)  # site means
  z <- scale(m)   # z-scores
  
  data.frame(mean_dist =m, z =z, is_outlier = abs(z) >= mult)
}

### try it
o <- outliers(spe, mult=2)
head(o,7)
which(o$is_outlier)

###  TEST VALIDITY OF SPECIES MATRIX
!anyNA(spe)   # expect TRUE, no missing value
all(rowSums(spe, na.rm = T) !=0)  # expect True, no empty sites
all(colSums(spe, na.rm = T) !=0)   # expect True, no empty species 

### VISUALIZE DATA

# Lets map the spatial coordinates

# spatial
plot(xy, pch=19, col='grey', xlab='Eastness', ylab='Northness')

# Plot the species abundance matrix as a heatmap.

# Species
?vegan::tabasco
(spe, col=get_palette(palette = "npg", 97))

Plot the soils (environmental) matrix as a heatmap.

# environment
vegan::tabasco(env, col=get_palette(palette = "npg", 97))

Plot the traits matrix as a heatmap.

# traits
vegan::tabasco(tra, col=get_palette(palette = "npg", 97))

Plot the phylogenetic tree.

# phylogeny
plot(phy, cex=0.6, no.margin=TRUE)


Gamma (regional) diversity
gamma <- sum(colSums(spe) > 0)
gamma

Alpha (per-site) diversity
alpha <- rowSums(spe > 0)
alpha    # within-site
avgalpha <- mean(rowSums(spe > 0))
avgalpha # average within-site

Beta (among-site) diversity: Whittaker's
gamma    <- sum(colSums(spe) > 0)
avgalpha <- mean(rowSums(spe > 0))
beta     <- gamma / avgalpha - 1
beta


### 1 -- proportion of zeros in the matrix
###   (independent of abundance)
eps      <- .Machine$double.eps # machine tolerance
propzero <- sum(spe < eps) / prod(dim(spe))
cat('Proportion of zeros in matrix:', propzero, '\n')

### 2 -- "dust bunny index" of McCune and Root (2015)
###   (integrates abundances)
dbi <- 1 - mean(as.matrix(vegan::decostand(spe, method='max')))
cat('Dust bunny index:', dbi, '\n')
d <- vegdist(spe, method='bray', binary=T)
str(d)
tabasco(as.matrix(d), col=get_palette(palette = "npg", 56)) 


### how many site-pairs share no species in common?
z <- vegan::no.shared(spe)
propnoshare <- sum(z) / length(z)
cat('Proportion of no-share sites:', propnoshare, '\n')














m1 <- metaMDS(D, k=2, try=50, trymax=51, trace=0)  # NMS
m2 <- cmdscale(D, k=2, add=T)               # PCoA
m3 <- prcomp(spe)                           # PCA
