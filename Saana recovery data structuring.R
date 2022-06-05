# Data structuring for Saana project
library(plyr)
library(tidyverse)
library(vegan)
library(lme4)
library(car)
library(FD)
library(lmerTest)
library(cowplot)

## Personal ggplot theme
theme_cw <- function () { 
  theme_bw(base_size=12) %+replace% 
    theme(
      panel.background = element_blank(), 
      plot.background = element_blank(), 
      axis.ticks = element_line(colour = "grey70", size = rel(0.5)),
      panel.grid.minor = element_blank(), 
      panel.grid.major.x = element_blank(),
      legend.background = element_blank(), 
      legend.key = element_blank(),
      strip.background = element_blank(), 
      strip.text=element_text(size=12),
      axis.text=element_text(size=12),
      complete = TRUE
    )
}


###  Ground cover data ------
ground.cover <- read.csv('Kilpisjarvi ground data.csv', stringsAsFactors = FALSE)
#View(ground.cover)
ground.cover$litter.depth.ave <- rowMeans(ground.cover[,c('litter.depth.1', 'litter.depth.2', 'litter.depth.3')])
ground.cover$height.ave <- rowMeans(ground.cover[,c('height.1', 'height.2', 'height.3')])
ground.cover$herbivore <- factor(ground.cover$herbivore)
ground.cover$fertilized <- factor(ground.cover$fertilized)
ground.cover$limed <- factor(ground.cover$limed)
ground.cover$block <- ground.cover$site %>% substring(1,1)
ground.cover$year <- '2019'

### Species cover 2019 ---------------
species.cover.2019 <- read.csv('Kilpisjarvi species data.csv', stringsAsFactors = FALSE) 
species.cover.2019[is.na(species.cover.2019)] <- 0 # replace blank NAs with 0
species.cover.2019.tidy <- species.cover.2019 %>% gather(key='plot', value='cover', -1) # tidy format
species.cover.2019.tidy$plot <- species.cover.2019.tidy$plot %>% substring(2) # remove leading X from plot names
species.cover.2019.tidy$year <- '2019'

non.species.columns <- c('site','plot','soil','herbivore','fertilized', 'limed',
                         'subplot', 'block','year')
plot.treatments <- ground.cover[,non.species.columns]
species.cover.2019.tidy <- species.cover.2019.tidy %>% left_join(plot.treatments, by=c('plot','year'))

# a version where all species are rows for using the vegan package
# species.cover.wide <- species.cover.2019.tidy %>% spread(species, cover)
# species.cover.wide$V1 <- NULL
# species.columns <- which(!(names(species.cover.wide) %in% non.species.columns))
# species.cover.wide$sp.richness <- species.cover.wide[,species.columns] %>% specnumber()
# species.cover.wide$sp.div <- species.cover.wide[,species.columns] %>% diversity(index='shannon')
# species.cover.wide$simpson <- species.cover.wide[,species.columns] %>% diversity(index='simpson')
# species.cover.wide$site.2 <- substr(species.cover.wide$site, 1, 1) # this is the same as block, why?

### Species cover 2004 ----------------------
species.cover.2004 <- read.csv('saana_vascular_cover_2004.csv', stringsAsFactors = FALSE)
species.cover.2004 <- species.cover.2004[,-c(59,60,61)]

species.cover.2004 <- species.cover.2004 %>% 
  dplyr::rename(block = site, plot = code, soil = casihabitat, herbivore = exclosure, 
         fertilized = manure, limed = liming)
species.cover.2004 <- species.cover.2004 %>% dplyr::rename(site = sitecode)

species.cover.2004$soil <- ifelse(species.cover.2004$soil == 'ca', 'fertile','infertile')
species.cover.2004$herbivore <- ifelse(species.cover.2004$herbivore == 'no', '1', '0')
species.cover.2004$fertilized <- ifelse(species.cover.2004$fertilized == 'yes', '1', '0')
species.cover.2004$limed <- ifelse(species.cover.2004$limed == 'yes', '1', '0')

species.cover.2004$plot <- species.cover.2004$plot %>%
  toupper() %>%
  gsub("EIM","",.) %>% 
  gsub(" ","",.)
names(species.cover.2004) <- tolower(names(species.cover.2004))
species.cover.2004$year <- as.character(species.cover.2004$year)
species.cover.2004$block <- as.character(species.cover.2004$block)

# shift to tidy format
non.species.columns.2 <- non.species.columns[!non.species.columns == 'subplot']
species.cover.2004.tidy <- species.cover.2004 %>% 
  gather(key = 'species', value = 'cover', -non.species.columns.2)

### Species cover 2010 ---------------------
species.cover.2010 <- read.csv('saana_vascular_cover_2010.csv', stringsAsFactors = FALSE)
species.cover.2010 <- species.cover.2010[,1:74]
species.cover.2010 <- species.cover.2010[,-c(8,10)]

names(species.cover.2010) <- gsub("first","", names(species.cover.2010))
species.cover.2010 <- species.cover.2010 %>% 
  dplyr::rename(block = site, plot = code, soil = casihabitat, herbivore = exclosure, 
         fertilized = manure, limed = liming)
species.cover.2010 <- species.cover.2010 %>% dplyr::rename(site = sitecode)

species.cover.2010$soil <- ifelse(species.cover.2010$soil == 'ca', 'fertile','infertile')
species.cover.2010$herbivore <- ifelse(species.cover.2010$herbivore == 'no', '1', '0')
species.cover.2010$fertilized <- ifelse(species.cover.2010$fertilized == 'yes', '1', '0')
species.cover.2010$limed <- ifelse(species.cover.2010$limed == 'yes', '1', '0')

species.cover.2010$plot <- species.cover.2010$plot %>%
  toupper() %>%
  gsub("EIM","",.) %>% 
  gsub(" ","",.)
names(species.cover.2010) <- tolower(names(species.cover.2010))
species.cover.2010$year <- as.character(species.cover.2010$year)
species.cover.2010$block <- as.character(species.cover.2010$block)

# tidy format
species.cover.2010.tidy <- species.cover.2010 %>% 
  gather(key = 'species', value = 'cover', -non.species.columns.2)

### Merge three years of species data --------------
species.cover.tidy <- full_join(species.cover.2004.tidy, species.cover.2010.tidy)
species.cover.tidy <- full_join(species.cover.tidy, species.cover.2019.tidy)
species.cover.tidy$block <- substr(species.cover.tidy$site, 1, 1)

## Cleaning up species code issues
## simple changes
species.cover.tidy$species[species.cover.tidy$species=='psealp'] <- 'psealb'
species.cover.tidy$species[species.cover.tidy$species=='solvig'] <- 'solvir'
species.cover.tidy$species[species.cover.tidy$species=='tolpus'] <- 'tofpus'

## combining the Cerastium alpinum subspecies
ceralp.subspecies <- c('ceralp1','ceralp2','ceralpalp','ceralpgla','ceralplan','cersp')
species.cover.tidy$species[species.cover.tidy$species %in% ceralp.subspecies] <- 'ceralp'

## NOTE: Lumping two sets of species into genus due to difficulty in telling them apart. 
## Not doing this at the earlier step just in case, but it is important here since we only 
## have trait data for one of each. 
## Combining Carex bigelowii with Carex vaginata, and Equisetum arvense with Equisetum pratense
species.cover.tidy$species[species.cover.tidy$species=='carbig'] <- 'carvag'
species.cover.tidy$species[species.cover.tidy$species=='equarv'] <- 'equpra'

species.cover.tidy.2 <- species.cover.tidy %>% 
  ddply(c('plot','year','site','block','soil','herbivore','fertilized','limed','species'),
        summarise,
        cover = sum(cover))

species.cover.wide <- species.cover.tidy.2 %>% spread(species, cover)
species.cover.wide$V1 <- NULL
species.cover.wide$subplot <- NULL
species.cover.wide[is.na(species.cover.wide)] <- 0
# double-check this later to make sure we're not just naming things differently
# especially the combined species
species.columns <- which(!(names(species.cover.wide) %in% non.species.columns))

## sort columns by year, site (which includes block and soil), fertilized, herbivore, limed
species.cover.wide <- species.cover.wide %>% arrange(year, site, fertilized, herbivore, limed)

## want to set up a similar dataframe for relative cover
species.cover.relative <- species.cover.wide
species.cover.relative$total <- rowSums(species.cover.relative[species.columns])
species.cover.relative[species.columns] <- species.cover.relative[species.columns]/species.cover.relative$total

### 2010 Litter depth data-------
litter.depth.2010 <- read.csv('litterdepth2010.csv', stringsAsFactors = FALSE)
litter.depth.2010 <- litter.depth.2010[,c(1:7, 13, 18, 20)]


litter.depth.2010 <- litter.depth.2010 %>% 
  dplyr::rename(block = site, plot = code, soil = casihabitat, herbivore = exclosure, 
         fertilized = manure, limed = liming, litter.cover = Litter, 
         litter.depth.ave = litterdepthmean)
litter.depth.2010 <- litter.depth.2010 %>% dplyr::rename(site = sitecode)

litter.depth.2010$soil <- ifelse(litter.depth.2010$soil == 'ca', 'fertile','infertile')
litter.depth.2010$herbivore <- ifelse(litter.depth.2010$herbivore == 'grazed', '1', '0')
litter.depth.2010$fertilized <- ifelse(litter.depth.2010$fertilized == 'fertilized', '1', '0')
litter.depth.2010$limed <- ifelse(litter.depth.2010$limed == 'limed', '1', '0')

litter.depth.2010$plot <- litter.depth.2010$plot %>%
  toupper() %>%
  gsub("EIM","",.) %>% 
  gsub(" ","",.)

litter.depth.2010$year <- as.character(litter.depth.2010$year)
litter.depth.2010$block <- substr(litter.depth.2010$site, 1, 1)

## add to ground cover data
ground.cover$year <- '2019'
ground.cover.simple <- ground.cover[,names(ground.cover) %in% names(litter.depth.2010)]
ground.cover.time <- rbind(ground.cover.simple, litter.depth.2010)

### Trait data ----------
trait.data <- read.csv('traitcalculations2011.csv', stringsAsFactors = FALSE)
trait.data.long <- trait.data %>% gather(key='species', value='value', -1)
trait.data.tidy <- trait.data.long %>% spread(key='trait', value='value')
trait.data.tidy$species <- tolower(trait.data.tidy$species)


## species present in both dataframe
species.only <- species.cover.relative[,species.columns] # currently relative cover, could also change it back
trait.data.2 <- filter(trait.data.tidy, species %in% names(species.only))
trait.data.3 <- trait.data.2[,c('CNratio','height','SLA')]
rownames(trait.data.3) <- trait.data.2$species

## species with non-zero cover in 2019
species.present <- names(which(colSums(species.only) > 0))
trait.data.4 <- trait.data.3[rownames(trait.data.3) %in% species.present,]
species.only.2 <- species.only[,rownames(trait.data.4)]

ft.cwd <- functcomp(trait.data.4, as.matrix(species.only.2))
ft.div <- dbFD(trait.data.4, as.matrix(species.only.2), corr='cailliez')

### Species functional groups
species.info <- read.csv('Kilpisjarvi species info.csv', stringsAsFactors = FALSE)
species.cover.tidy.3 <- species.cover.tidy.2 %>%
  left_join(species.info, by=c('species'='species.code')) %>%
  dplyr::rename(species.name = species.y)

#write.csv(trait.data.tidy, file='saana_trait_data.csv')
# manually added additional trait data from Elina's 2017 paper
trait.data.tidy.2 <- read.csv('saana_trait_data_2.csv', stringsAsFactors = FALSE)
trait.data.tidy.2 <- trait.data.tidy.2[,c('species','height','CNratio','SLA')]

trait.data.tidy.2$height.s <- scale(trait.data.tidy.2$height)
trait.data.tidy.2$CNratio.s <- scale(trait.data.tidy.2$CNratio)
trait.data.tidy.2$SLA.s <- scale(trait.data.tidy.2$SLA)

species.cover.traits <- species.cover.tidy.3 %>% left_join(trait.data.tidy.2, by='species')




### Soil data --------------
soil.data <- read.csv('Saana_soil_all_raw.csv', stringsAsFactors = FALSE)
soil.data$herbivore <- ifelse(soil.data$herbivore == 'grazed', '1', '0')
soil.data$fertilized <- ifelse(soil.data$fertilized == 'fertilized', '1', '0')
soil.data$limed <- ifelse(soil.data$limed == 'limed', '1', '0')

soil.data$plot <- soil.data$plot %>%
  toupper() %>%
  gsub("EIM","",.) %>% 
  gsub(" ","",.)

soil.data$site <- substr(soil.data$plot, 1, 2)

soil.data$c.n.ratio <- soil.data$c.perc / soil.data$n.perc
soil.data$nh4[soil.data$nh4 == "<0,0050"] <- 0
soil.data$nh4 <- as.numeric(soil.data$nh4)

### Plotting function -----
plot.saana <- function(df, y.id, y.se, y.lab) {
  y.id <- enquo(y.id)
  y.se <- enquo(y.se)
  ggplot(df, aes(x = soil, y = !!y.id, 
                 ymin = !!y.id - !!y.se,
                 ymax = !!y.id + !!y.se,
                 shape = herbivore, color = fertilized)) +
    geom_point(size=4, position=position_dodge(width=0.3)) +
    geom_errorbar(width=0.2, position=position_dodge(width=0.3)) +
    theme_cw() +
    xlab('Soil type') +
    ylab(y.lab) +
    scale_color_manual(values=c('seagreen','darkorange'), name='', 
                       labels=c('control','fertilized'),
                       guide = guide_legend(reverse=TRUE)) +
    scale_shape_discrete(name='',labels=c('excluded','herbivory'))
}

plot.saana.yr <- function(df, y.id, y.se, y.lab, l.pos = c(0.9, 0.9)) {
  y.id <- enquo(y.id)
  y.se <- enquo(y.se)
  ggplot(df, aes(x = year, y = !!y.id, 
                 ymin = !!y.id - !!y.se,
                 ymax = !!y.id + !!y.se,
                 shape = herbivore, color = fertilized)) +
    facet_grid(.~soil, labeller = as_labeller(soil.names)) +
    geom_point(size=4, position=position_dodge(width=0.3)) +
    geom_errorbar(width=0.2, position=position_dodge(width=0.3)) +
    theme_cw() +
    xlab('Year') +
    ylab(y.lab) +
    theme(legend.position = l.pos) +
    scale_color_manual(values=c('seagreen','darkorange'), name='', 
                       labels=c('control','fertilized'),
                       guide = guide_legend(reverse=TRUE)) +
    scale_shape_discrete(name='',labels=c('excluded','herbivory'))
}


plot.saana.simple <- function(df, y.id, y.se, y.lab, l.pos = c(0.9, 0.9)) {
  y.id <- enquo(y.id)
  y.se <- enquo(y.se)
  ggplot(df, aes(x = soil, y = !!y.id, 
                 ymin = !!y.id - !!y.se,
                 ymax = !!y.id + !!y.se,
                 color = fertilized)) +
    geom_point(size=4, position=position_dodge(width=0.3)) +
    geom_errorbar(width=0.2, position=position_dodge(width=0.3)) +
    theme_cw() +
    xlab('Soil type') +
    ylab(y.lab) +
    theme(legend.position = l.pos) +
    scale_color_manual(values=c('seagreen','darkorange'), name='', 
                       labels=c('control','fertilized'),
                       guide = guide_legend(reverse=TRUE))
}

plot.saana.yr.simple <- function(df, y.id, y.se, y.lab, l.pos = c(0.9, 0.9)) {
  y.id <- enquo(y.id)
  y.se <- enquo(y.se)
  ggplot(df, aes(x = year, y = !!y.id, 
                 ymin = !!y.id - !!y.se,
                 ymax = !!y.id + !!y.se,
                 color = fertilized)) +
    facet_grid(.~soil, labeller = as_labeller(soil.names)) +
    geom_point(size=4, position=position_dodge(width=0.3)) +
    geom_errorbar(width=0.2, position=position_dodge(width=0.3)) +
    theme_cw() +
    xlab('Year') +
    ylab(y.lab) +
    theme(legend.position = l.pos) +
    scale_color_manual(values=c('seagreen','darkorange'), name='', 
                       labels=c('control','fertilized'),
                       guide = guide_legend(reverse=TRUE))
}

plot.saana.limed <- function(df, y.id, y.se, y.lab) {
  y.id <- enquo(y.id)
  y.se <- enquo(y.se)
  ggplot(df, aes(x = limed, y = !!y.id, 
                 ymin = !!y.id - !!y.se,
                 ymax = !!y.id + !!y.se,
                 shape = herbivore, color = fertilized)) +
    facet_grid(.~soil, labeller = as_labeller(soil.names)) +
    geom_point(size=4, position=position_dodge(width=0.3)) +
    geom_errorbar(width=0.2, position=position_dodge(width=0.3)) +
    theme_cw() +
    scale_x_discrete(name = 'Liming', labels = c('Control', 'Limed')) +
    scale_color_manual(values=c('black','grey50'), name='', labels=c('control','fertilized')) +
    scale_shape_discrete(name='',labels=c('excluded','herbivory')) +
    ylab(y.lab)
}