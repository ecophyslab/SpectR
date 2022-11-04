require(Rmisc)
require(spectrolab)
require(pls)
require(chillR)
require(paletteer)
require(caret)
require(ggpubr)

dir <-  dirname(dirname(getwd()))
fig_dir <-paste0(dir,"/outputs/") #where to store plots

######Load data
data_PLSR <- read.csv(paste0(dir,'/data/zoysia_physpec.csv')) #loading processed spec phys data

#play with sample size if you want
# data_PLSR <- data_PLSR[which(data_PLSR$genotype == "1410"),] 
# data_PLSR <- data_PLSR[sample(1:nrow(data_PLSR),floor(0.30*nrow(data_PLSR))),]
rownames(data_PLSR) <- as.character(seq(1,length(rownames(data_PLSR)),1))

validata <- read.csv(paste0(dir,'/data/zoysia_physpec_validata.csv'))  #loading independent dataset


#specify wavelengths an resolution
min_edge <- 400 #min wavelength to keep
max_edge <- 2400 #max wavelength to keep
res <- 1 #spectral resolution in nanometers

#specify number of iterations
iter <- 10
rmsep_plot <- F #set true to see component selection plots

#specify colors for treatments and reference treatment
colors <- paletteer_d("ggthemes::Tableau_20")[1:length(levels(as.factor(data_PLSR$treatment)))]
ref_trt <- "control" #choose a treatment that contains healthy looking spectra



#####PLSR start

#Removing X from colnames
colnames(data_PLSR)[which(colnames(data_PLSR)=="X400"):length(data_PLSR)] <- gsub("X","",colnames(data_PLSR[which(colnames(data_PLSR)=="X400"):length(data_PLSR)]))

colnames(validata)[which(colnames(validata)=="X400"):length(validata)] <- gsub("X","",colnames(validata[which(colnames(validata)=="X400"):length(validata)]))


#standardize time format
data_PLSR$date <-  as.Date(data_PLSR$date)
validata$date <-  as.Date(validata$date)

wvls <- as.character(seq(min_edge,max_edge,res)) #defining wavelengths and resolution to use
min_col <- which(colnames(data_PLSR)==wvls[1]) #defining first wavelength column in dataset
max_col <- which(colnames(data_PLSR)==wvls[length(wvls)]) #defining first wavelength column in dataset


fvfm_max <- max(data_PLSR$mean_fvfm, na.rm = T)
fvfm_min <- min(data_PLSR$mean_fvfm, na.rm = T)


options(scipen = 999) #makes R report decimals properly rather than using 1*e^x numbers
m <- list() #Model iterations are stored here
mods <- list() #prediction iterations are stored here
ncomp <- list() #number of components used in each iteration is stored here
R2<-vector() #R square values of each iteration are stored here
pv <-vector() #p values of each iteration are stored here
slope<-vector() #slope values of each iteration are stored here
Error<-vector() #RSME of each model iteration are stored here
Error_percent<-vector() #Relative RSME to allow comparisons among variables
bias <- vector() #Bias of each iteration is stored here
a<-1 #starts a counter for the loop
t<-list() #matrix to compile predicted vs measured values
var_imp <- list() #wavelength importance values of each iteration are stored here

#External dataset model testing 
mods_ex <- list() 
R2_ex<-vector() 
pv_ex <-vector() #p values of each iteration are stored here
slope_ex<-vector() 
Error_ex<-vector() 
Error_percent_ex<-vector() 
t_ex<-list()


for (i in 1:iter){
  set.seed(i) #This reinforces random sampling of the data for the k groups
  
  test_sample <-sample(1:nrow(data_PLSR),floor(0.20*nrow(data_PLSR))) #Sets a random 20% of observations aside to test the model after training
  

  train_sample <-setdiff(1:nrow(data_PLSR),test_sample)#grabs spectra for all the observations not present in the test group to train the model.
  
  
  train_set<-data_PLSR[train_sample,] #Uses the rownames selected in the previous step to gather the data_PLSR for training
  test_set<-data_PLSR[test_sample,] #Uses the rownames selected in the previous stem to gather the data_PLSR for testing
  
  
  fvfm_train<-plsr(train_set$mean_fvfm~value(as_spectra(train_set[min_col:max_col])),
                   ncomp=10,validation="LOO") #Makes models with up to 10 independent components (orthogonal combinations of the predictors) using the desired physiological variable as a response and all the spectra as predictors.
  
  #fvfm_train uncheck to see model object
  
  m[[a]] <- fvfm_train #Stores the resulting models

  
  ncomp[[a]] <- selectNcomp(fvfm_train, method = "onesigma", plot = rmsep_plot) #finds how many components are needed to minimize the RSME of the model. You want to minimize the components to avoid overfitting while maximizing the amount of error reduced. Set plot to TRUE the first time to find out the optimal number of components. Then to FALSE once you have set the correct number of components on the next line.

  mods[[a]]<-predict(fvfm_train,ncomp=5,newdata=value(as_spectra(test_set[min_col:max_col]))) #Uses the model to predict physiological variables for the test observations.Must adjust the # of components (ncomp).

  R2[a]<-summary(lm(test_set$mean_fvfm~mods[[a]]))$r.squared #saves the R square of the model
  
  pv[a]<-anova(lm(test_set$mean_fvfm~mods[[a]]))$`Pr(>F)`[1] #saves the p-value of the model
  
  slope[a]<-summary(lm(test_set$mean_fvfm~mods[[a]]))$coefficients[2] #saves the slope of the model

  Error[a]<-RMSEP(mods[[a]],test_set$mean_fvfm,na.rm = T)  #saves the RMSE of the prediction

  Error_percent[a] <- Error[a]/(max(train_set$mean_fvfm,na.rm = T)-min(train_set$mean_fvfm,na.rm = T))*100 #makes RMSE unitless to compare variables if needed.
  
  #Collecting the data_PLSR used and models so that I can plot the performance plots later on.
  t[[a]]<-as.data.frame(mods[[a]])
  t[[a]]$mean_fvfm <- test_set$mean_fvfm
  t[[a]]$treatment <- test_set$treatment
  t[[a]]$row <- rownames(test_set)
  
  bias[[a]] <- mean(t[[a]]$mean_fvfm, na.rm = T)-mean(t[[a]][,1], na.rm = T) #Calculating bias
  
  
  var_imp[a]<-varImp(fvfm_train)  #Collecting wavelength importance across model iterations

  a <- a+1  
  
    if (i == iter) {
    print('Done!')
    
    } else {
    print(paste('Round',i,sep=' ')) #Tells you which round of the iteration you are entering
    }
  }


fvfm_C <- mean(reshape2::melt(ncomp)$value)  # calculates the average number of components used across the 100 iterations

max(reshape2::melt(ncomp)$value) # tells you how many components where used in the most complex model across the 100 iterations
plot(density(reshape2::melt(ncomp)$value),main = 'Distribution of components') #shows the distribution of component values. The narrower and the closer to one the better.

max(R2) # tells you the highest R square found across the 100 iterations
mean(R2) # tells you the average R square found across the 100 iterations

plot(density(R2),main = 'Distribution of R2 values')  #shows the distribution of R square values. The narrower and the higher the better.

#Stores the diagnostic plots to see if the model is well constructed
pdf(paste0(fig_dir,"mean_fvfm_diagnostics.pdf"),width = 15,height = 15) #creating a pdf

layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))
plot(density(reshape2::melt(ncomp)$value),main = 'Distribution of components')
plot(density(R2),main = 'Distribution of R2 values')

dev.off()

fvfm_pv <- median(pv) #stores the median p-value found across the 100 iterations. I like using this value better than the mean.

fvfm_R2 <- median(R2) #  #stores the median R square found across the 100 iterations. I like using this value better than the mean.

fvfm_slope <- median(slope) #  #stores the median slope found across the 100 iterations. I like using this value better than the mean.

median(Error)# #reports the median RSMEP found across the 100 iterations. I like using this value better than the mean.

fvfm_E <- median(Error_percent)# #stores the median percent RSMEP found across the 100 iterations. I like using this value better than the mean.


#Variable importance
Importance <- as.data.frame(matrix(ncol=length(min_col:max_col),nrow=iter)) #preps wavelength importance data for plotting

colnames(Importance) <- colnames(test_set[min_col:max_col]) #brings in the wavelength names

for (i in 1:length(var_imp)) {
  
  Importance[i,]<-unlist(var_imp[i])

}

Mean_imp <- colMeans(Importance, na.rm = T) #calculates the average importance of each wavelength across iterations.

new_Mean_imp <- c(rep(NA,as.numeric(colnames(test_set)[min_col])-1), Mean_imp)



#Plot average variable importance through the iterations


pdf(paste0(fig_dir,"mean_fvfm_importance.pdf"),width = 15,height = 15) #creating a pdf

    layout(matrix(c(1,1), 1, 2, byrow = TRUE))
    
    plot(abs(new_Mean_imp),type="n",xlim = c(min_edge,max_edge),xaxt = "n", ylab = "", xlab = "",xaxt="n",yaxt="n", col="darkgrey", col.axis    ="#666666", col.lab="darkgrey", main = "Predicting Fv/Fm all treatments",cex.main= 2)
    axis(side = 2,col="darkgrey", col.axis="#666666", col.lab="darkgrey", 
         font.axis=13, font=1,cex.axis=2, cex.lab=2)
    #abline(0, 0)
    lines(abs(new_Mean_imp),col = 'black',lwd = 2)
    par(new = TRUE)
    plot(as_spectra(test_set[which(test_set$treatment==ref_trt & !is.na(test_set[min_col])),][1,min_col:max_col],name_idx = 1),col = 'darkgrey',lwd = 1, lty = 1,yaxt = "n",ylab = "Importance", xlab = "Wavelength    (nm)",xaxt="n",bg="transparent",cex.lab=2) #M0
    axis(side=1, col="darkgrey", col.axis="#666666", col.lab="darkgrey", 
         font.axis=13, font=1,at=seq(0, 2500, by=400),cex.axis=2, cex.lab=2)
    box(col="darkgrey")

dev.off()


#predicted vs Measured plot preping
for (i in 1:length(t)) {
  t[[i]]$cat <- i
}

#preping pred vs meas data by treatment
D <- cbind(do.call(rbind,t))
colnames(D) <- c("predicted_mean_fvfm","mean_fvfm","treatment","row","cat")
D$cat <- as.factor(D$cat)
D <-D[complete.cases(D),]
D$row <- as.factor(D$row)

D_points <- aggregate(D,by = list(D$row),FUN = "mean",na.rm=T)
D_points_CI<-group.CI(predicted_mean_fvfm~row,D)
D_points_CI$mean_fvfm <-D_points$mean_fvfm
D_points_CI$row <-as.numeric(as.character(D_points_CI$row))
D_points_CI <-D_points_CI[order(D_points_CI$row),]
D_points_CI$treatment <- data_PLSR$treatment[D_points_CI$row]


#Across treatment
pdf(paste0(fig_dir,"mean_fvfm_test_all.pdf"),width = 7,height = 7)
  ggplot(D_points_CI, aes(y = mean_fvfm,x = predicted_mean_fvfm.mean,color=treatment))+
    geom_point(aes(y=mean_fvfm,x=predicted_mean_fvfm.mean),alpha = 1,size = 2.5)+
    geom_errorbar( aes(xmin = predicted_mean_fvfm.lower, xmax =  predicted_mean_fvfm.upper),width = 0.001) +
    geom_line(aes(y=mean_fvfm,x=predicted_mean_fvfm.mean),stat="smooth",method = "lm", formula = y ~ x,size = 1.3,linetype ="solid"  ,alpha = 1,colour ="darkred")+
    geom_abline(intercept = 0, slope = 1)+
    scale_y_continuous(n.breaks = 10,limits = c(min(D_points_CI$predicted_mean_fvfm.lower[!is.nan(D_points_CI$predicted_mean_fvfm.lower)]),max(D_points_CI$predicted_mean_fvfm.upper[!is.nan(D_points_CI$predicted_mean_fvfm.upper)])))+
    scale_x_continuous(n.breaks = 10,limits = c(min(D_points_CI$predicted_mean_fvfm.lower[!is.nan(D_points_CI$predicted_mean_fvfm.lower)]),max(D_points_CI$predicted_mean_fvfm.upper[!is.nan(D_points_CI$predicted_mean_fvfm.upper)])))+
    scale_color_manual(values=colors)+
    theme_minimal()+
    theme(panel.border = element_rect(colour = "grey", fill=NA, size=1))+
  ylab("Measured Fv/Fm")+xlab("Predicted Fv/Fm")+
  theme(aspect.ratio= 1,
        legend.position = "top",
        plot.title = element_text(size="10", color = "#666666", hjust = .05, vjust = -50)) +
  ggtitle(paste(
    paste("Median p-value",ifelse(median(pv)> 0.001,paste("=",round(median(pv),3)),paste("< 0.001"))),"\n",
    paste("Median R2 =",round(median(R2),2)),"\n",
    paste("Median slope =",round(median(slope),2)),"\n",
    paste("Median % RMSEP =",round(median(Error_percent),2)),"\n",
    paste("Median bias =",round(median(bias),4)),"\n",
    paste("# Components =",round(mean(reshape2::melt(ncomp)$value))))) 
  # expand_limits(x = 0, y = 0)


dev.off()

#You can use the following lines to use an independent dataset not used in training nor testing to validate  your model.
fvfm_validated<-predict(fvfm_train,ncomp=round(mean(reshape2::melt(ncomp)$value)),newdata=value(as_spectra(validata[which(colnames(validata)==wvls[1]):which(colnames(validata)==wvls[length(wvls)])])))
fvfm_validated <- as.data.frame(fvfm_validated)


data_compare <- validata[0:(which(colnames(validata)==wvls[1])-1)] #adjust columns with phys variables as needed
data_compare$predicted_mean_fvfm <- fvfm_validated$`train_set$mean_fvfm.5 comps` #must adjust the # of components

R2_ex<-summary(lm(data_compare$mean_fvfm~data_compare$predicted_mean_fvfm))$r.squared #saves the R square of the model

pv_ex<-anova(lm(data_compare$mean_fvfm~data_compare$predicted_mean_fvfm))$`Pr(>F)`[1] #saves the p-value of the model

slope_ex<-summary(lm(data_compare$mean_fvfm~data_compare$predicted_mean_fvfm))$coefficients[2] #saves the slope of the model

Error_ex<-RMSEP(data_compare$predicted_mean_fvfm,data_compare$mean_fvfm,na.rm = T)  #saves the RMSE of the prediction

Error_percent_ex <- Error_ex/(max(data_compare$mean_fvfm,na.rm = T)-min(data_compare$mean_fvfm,na.rm = T))*100 #makes RMSE unitless to compare variables if needed.

ggplot(data_compare, aes(y = mean_fvfm,x = predicted_mean_fvfm,color=treatment))+
  geom_point(aes(y=mean_fvfm,x=predicted_mean_fvfm),alpha = 1,size = 2.5)+
  # geom_errorbar( aes(xmin = predicted_mean_fvfm.lower, xmax =  predicted_mean_fvfm.upper),width = 0.001) +
  geom_line(aes(y=mean_fvfm,x=predicted_mean_fvfm),stat="smooth",method = "lm", formula = y ~ x,size = 1.3,linetype ="solid"  ,alpha = 1,colour ="darkred")+
  geom_abline(intercept = 0, slope = 1)+
  scale_y_continuous(n.breaks = 10,limits = c(min(data_compare$predicted_mean_fvfm[!is.nan(data_compare$predicted_mean_fvfm)]),max(data_compare$predicted_mean_fvfm[!is.nan(data_compare$predicted_mean_fvfm)])))+
  scale_x_continuous(n.breaks = 10,limits = c(min(data_compare$predicted_mean_fvfm[!is.nan(data_compare$predicted_mean_fvfm)]),max(data_compare$predicted_mean_fvfm[!is.nan(data_compare$predicted_mean_fvfm)])))+
  scale_color_manual(values=colors)+
  theme_minimal()+
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=1))+
  ylab("Measured Fv/Fm")+xlab("Predicted Fv/Fm")+
  theme(aspect.ratio= 1,
        legend.position = "top",
        plot.title = element_text(size="10", color = "#666666", hjust = .05, vjust = -50)) +
  ggtitle(paste(
    paste("Median p-value",ifelse(median(pv_ex)> 0.001,paste("=",round(median(pv_ex),3)),paste("< 0.001"))),"\n",
    paste("Median R2 =",round(median(R2_ex),2)),"\n",
    paste("Median slope =",round(median(slope_ex),2)),"\n",
    paste("Median % RMSEP =",round(median(Error_percent_ex),2)),"\n",
    paste("Median bias =",round(median(bias),4)),"\n",
    paste("# Components =",round(mean(reshape2::melt(ncomp)$value))))) 
# expand_limits(x = 0, y = 0)



#fvfm over time real  
RealvsTime <- ggplot(data_compare, aes(y = mean_fvfm,x = date,color=treatment))+
  geom_point(aes(y=mean_fvfm,x=date),alpha = 1,size = 2.5)+
  geom_smooth(aes(y=mean_fvfm,x=date,color=treatment,fill=treatment),stat="smooth",method = "loess", formula = y ~ x,size = 1.3,linetype ="solid"  ,alpha = 0.3)+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  theme_minimal()+
  ylab("Measured Fv/Fm")+xlab("Date of experiment")


#fvfm  over time spectra
SpectralvsTime <- ggplot(data_compare, aes(y = predicted_mean_fvfm,x = date,color=treatment))+
  geom_point(aes(y=predicted_mean_fvfm,x=date),alpha = 1,size = 2.5)+
  geom_smooth(aes(y=predicted_mean_fvfm,x=date,color=treatment,fill=treatment),stat="smooth",method = "loess", formula = y ~ x,size = 1.3,linetype ="solid"  ,alpha = 0.3)+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  theme_minimal()+
  ylab("Reflectance-predicted Fv/Fm")+xlab("Date of experiment")


pdf(paste0(fig_dir,"mean_fvfm_time_spectral.pdf"),width = 16,height = 12)

ggarrange(RealvsTime,SpectralvsTime,ncol = 1)

dev.off()
