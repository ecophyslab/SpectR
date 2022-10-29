#loading data
dir <- dirname(dirname(getwd()))

phys.data <- read.csv(paste0(dir,'/data/zoysia_physiology.csv')) # point to phys data
processed.spectra <- read.csv(paste0(dir,'/data/spectra/zoysia_japonica/2022_07_10/zoysia_processed_2022_07_10.csv'))

#Removing X from colnames
colnames(processed.spectra)[which(colnames(processed.spectra)=="X400"):length(processed.spectra)] <- gsub("X","",colnames(processed.spectra[which(colnames(processed.spectra)=="X400"):length(processed.spectra)]))

#standardize time format
processed.spectra$date <-  as.Date(processed.spectra$date)
phys.data$date <-  as.Date(phys.data$date)



#Merging
db <- merge( phys.data,processed.spectra,
             by.x = c( "plant_number", 'treatment', 'date','genotype'), 
             by.y = c( "plant_id", 'treatment', 'date','genotype'), 
             all = T) 




#### Example of spectral index calculation
require(hsdar)

test<-db[which(colnames(db)=="400"):length(db)]

test<-as.matrix(test)
wavelength <- as.numeric(colnames(test))
test_lib <- speclib(spectra = test,wavelength = as.numeric(colnames(test)))

db$PRI <- vegindex(test_lib, "PRI")
db$NDVI <- vegindex(test_lib, "NDVI")
db$Carter <- vegindex(test_lib, "Carter")


#### Ploting indexes vs physiology
require(ggplot2)
require(ggpubr)
theme_set(theme_minimal())

p1 <- ggplot(db,aes(x = mean_fvfm, y = PRI))+ geom_point(aes(x = mean_fvfm, y = PRI,color = treatment))+geom_smooth(aes(color = treatment, fill = treatment),method = "lm")+ scale_x_continuous(limits = c(0,.8))+ theme(legend.position = "top")

p2 <- ggplot(db,aes(x = mean_fvfm, y = NDVI))+ geom_point(aes(x = mean_fvfm, y = NDVI,color = treatment))+geom_smooth(aes(color = treatment, fill = treatment),method = "lm")+ scale_x_continuous(limits = c(0,.8))+ theme(legend.position = "none")

p3 <- ggplot(db,aes(x = mean_fvfm, y = Carter))+ geom_point(aes(x = mean_fvfm, y = Carter,color = treatment))+geom_smooth(aes(color = treatment, fill = treatment),method = "lm")+ scale_x_continuous(limits = c(0,.8))+ theme(legend.position = "none")


p4 <- ggplot(db,aes(x = sd_fvfm, y = PRI))+ geom_point(aes(x = sd_fvfm, y = PRI,color = treatment))+geom_smooth(aes(color = treatment, fill = treatment),method = "lm")+ scale_x_continuous(limits = c(0,.2))+ theme(legend.position = "none")

p5 <- ggplot(db,aes(x = sd_fvfm, y = NDVI))+ geom_point(aes(x = sd_fvfm, y = NDVI,color = treatment))+geom_smooth(aes(color = treatment, fill = treatment),method = "lm")+scale_x_continuous(limits = c(0,.2))+ theme(legend.position = "none")

p6 <- ggplot(db,aes(x = sd_fvfm, y = Carter))+ geom_point(aes(x = sd_fvfm, y = Carter,color = treatment))+geom_smooth(aes(color = treatment, fill = treatment),method = "lm")+ scale_x_continuous(limits = c(0,.2))+ theme(legend.position = "none")


ggarrange(p1,p4,p2,p5,p3,p6, ncol = 2,nrow = 3)