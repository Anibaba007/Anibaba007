#gbif manipulation


nc<-read.csv('gbif_input_eg.csv',sep=';')
library(rgbif);library(tidyverse)

sp.ex<-name_backbone(name='Solidago gigantea')
View(sp.ex)

check_gbif_fun<-function(x) name_backbone(name = x, kingdom = 'plants')

gbifkey <- lapply(nc$species, check_gbif_fun)
for(i in 1:length(gbifkey))gbifkey[[i]]<-as.data.frame(gbifkey[[i]])
gbifkey <- bind_rows(gbifkey)
gbifkey$orig<-nc$species
write.table(gbifkey,'gbiofexample_output.csv',sep=';', row.names = F)
