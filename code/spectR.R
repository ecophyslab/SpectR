#libraries
library(spectrolab)

#directory names

dir <- dirname(dirname(getwd()))
output_fn <- "zoysia_processed" #name of the output csv file

#folder names
spectra_fl <- 'spectra' #Root folder containing species folders
Species <- 'zoysia_japonica'#species folder containing sampling folders
Sampling <- '2022_07_10' #Sampling folder containing spectra files
svc_bn <- "Zoysia_2022710" #the basename under which the sig files are
meta_fn <-  paste0(dir,"/data/zoysia_spectra_meta.csv") #metadata filename

#settings
match <- c(990,1900) #sensor jumps
min_edge <- 400 #min wavelength to keep
max_edge <- 2400 #max wavelength to keep
res <- 1 #spectral resolution in nanometers
normalize <- FALSE #set to true if spectra should be normalized 



####The spectR starts####
spectra_meta <- read.csv(meta_fn) #File with metadata linking spectra files to corresponding Plant IDs and Treatments

#construction of full file names for spectra files
for (i in 1: length(spectra_meta$file_id)) {
  
  if (spectra_meta$file_id[i] < 10) {
    
    spectra_meta$Full_file_id[i] <- paste(spectra_meta$basename_id[i], '.000',spectra_meta$file_id[i],'.sig',sep = '')
    
  } else if (spectra_meta$file_id[i] >= 10 && spectra_meta$file_id[i] < 100) {
    
    spectra_meta$Full_file_id[i] <- paste(spectra_meta$basename_id[i], '.00',spectra_meta$file_id[i],'.sig',sep = '')
    
  } else if (spectra_meta$file_id[i] >= 100 && spectra_meta$file_id[i] < 1000) {
    
    spectra_meta$Full_file_id[i] <- paste(spectra_meta$basename_id[i], '.0',spectra_meta$file_id[i],'.sig',sep = '')
    
  } else if (spectra_meta$file_id[i] >= 1000) {
    
    spectra_meta$Full_file_id[i] <- paste(spectra_meta$basename_id[i], '.',spectra_meta$file_id[i],'.sig',sep = '')
  }
}

# path output analysis where csv file will be created
dirOut<-paste(dir,"/data/", spectra_fl, "/", Species,"/", Sampling,sep='')

# path reflectance files
dirspec<-paste(dir,"/data/", spectra_fl, "/", Species,"/", Sampling,sep='')



########################################
## 2. Start of loop for reading in files
########################################
lf = list.files(dirspec)[grep('.sig',list.files(dirspec))]

spec<- list()
spec_meta<-list()

for(i in 1:length(lf)){
  
  dirin = paste(dirspec, lf[i], sep = "/")
  
  ############################################
  # 3. Importing spectra
  ############################################
  spec[i]  = read_spectra(dirin, "sig",extract_metadata = T)
  spec_meta[[i]] <- meta(read_spectra(dirin, "sig",extract_metadata = T))
  
  
}

df <- data.frame(matrix(unlist(spec), nrow=length(lf), byrow=T),stringsAsFactors=FALSE)
colnames(df)<-colnames(as.data.frame(read_spectra(dirin, "sig")))[-1]

df$Full_file_id <- lf
df<-df[,c(ncol(df),1:(ncol(df)-1))] #placing file name at first column

spec <- as_spectra(df,name_idx = 1)
plot(spec[1:4])

spec<-match_sensors(spec,match) #eliminates jump due to sensor shifts
spec = spectrolab::resample(spec, seq(min_edge, max_edge, by = res) )
plot(spec[2:9])
spec_vn <- normalize(spec)
plot(spec_vn[2:9])

if (normalize ==T) {
  
  dati <- as.data.frame(spec_vn)
  plot(spec_vn,main = paste("All data normalized"))
} else  {
  
  dati <- as.data.frame(spec)
  plot(spec,main = paste("All data non-normalized"))
}

colnames(dati)[1]<- 'Full_file_id'
Mspectra<-merge(spectra_meta,dati, by = 'Full_file_id')  

Mspectra$file_id <- NULL #Removing meaningless variables
write.csv(Mspectra, paste(dirOut, paste(output_fn,"_",Sampling,".csv",sep=''), sep = "/"), row.names = FALSE)

spec_meta <- do.call("rbind", spec_meta)
write.csv(spec_meta, paste(dirOut, paste(output_fn,"_metadata_",Sampling,".csv",sep=''), sep = "/"), row.names = FALSE)

print("The spectR has performed its spooky magic...")




#### Example of spectral index calculation
library(hsdar)

test<-Mspectra[which(colnames(Mspectra)=="400"):length(Mspectra)]

test<-as.matrix(test)

wavelength <- as.numeric(colnames(test))

test_lib <- speclib(spectra = test,wavelength = as.numeric(colnames(test)))

allindices <- vegindex()

indexes <- vegindex(test_lib, allindices)

# ?soilindex

all_db <- cbind(Mspectra,indexes)