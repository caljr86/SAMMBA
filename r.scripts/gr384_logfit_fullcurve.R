###########################################################################
#
#   This script is a big loop with three parts:
#   1 - Create gc384 file with raw particle analysis
#   2 - Create gc384_apw: an aggregation of all gametophytes per well (APW: area per well)
#   3 - Create plateplot: gr_fit is created based on gc384_apw, with growth rates from higher slope of 6 windows 
#   4 - 
#
###########################################################################

#Load libraries
library(growthrates);library(ggplot2);library(dplyr);library(combinat)

for(plt in c(1)){
  start_time <- Sys.time()
plate <- paste0("plate", plt)
setwd("~/Dropbox/1_posdoc/1_CCMAR/experiments/384WP/myrsini")
setwd(paste0(getwd(), "/", plate, "/"))
output <- paste0(getwd(), "/")
design <- read.table(file = "plate_design.csv", header = T, sep = ",")
nrow(design)

gc384 <- NULL
for(f in list.dirs(recursive = F)){
date <- unlist(strsplit(f, "/"))[2]

for(i in list.files(paste0(f, "/segmented/results/"))){
  # Skip to the next folder if no files are found
  if (length(f) == 0) {
    #cat("Skipping empty folder:", folder, "\n")
    next  # Go to the next iteration of the loop
  }
  area <- read.table(file = paste0(f, "/segmented/results/", i), header =T, sep = ",")
  area <- area[,1:2]
  colnames(area) <- c("gametophyte", "area")
  photo <- as.character(unlist(strsplit(i, "_segmented"))[1])
  photo <- as.numeric(unlist(strsplit(photo, "_"))[4])
  well1 <- as.character(unlist(strsplit(i, "_"))[2])
  well2 <- design[design$PHOTO == photo, ][,3]
  
    if(well1 != well2){
    print("Warning: Well divergence")
    }
  metadata <- design[design$WELL == well1,][,c(1,2,5,6)]
  metadata2 <- metadata[rep(1, nrow(area)),]
  plate_df <- rep(plate, nrow(area))
  date1 <- rep(date, times = nrow(area))
  img <- cbind(photo, well2, area, plate_df, metadata2, date1)
  colnames(img) <- c("photo", "well", "gametophyte", "area", "plate", "row", "col", "temperature", "sex", "date" )
  gc384 <- rbind(gc384, img)
}
}

gc384$gametophyte <- as.numeric(gc384$gametophyte)
gc384$col <- as.numeric(gc384$col)
gc384$date <- as.Date(as.character(gc384$date), format="%y%m%d")

#Reorder wells sequence
gc384 <- gc384[order(gc384$row, gc384$col),]
gc384$well <- factor(gc384$well, levels = unique(paste0(gc384$row, gc384$col)))

head(gc384)
nrow(gc384)
str(gc384)



#######################################
####### Total Area per well ###########
#######################################

gc384_apw <- aggregate(area ~ well + date + temperature + plate + row + col + photo, data = gc384, sum)
gc384_apw$counts <- aggregate(gametophyte ~ well+date+ temperature+plate+row+col+photo, data = gc384, length)$gametophyte
day0 <- unique(na.omit(gc384_apw$date))[1]
gc384_apw$day <- as.numeric(gc384_apw$date - as.Date(day0))
gc384_apw <- gc384_apw[order(gc384_apw$row, gc384_apw$col),]

if (!exists("gc384_apw_list")) {
  gc384_apw_list <- list()
}

gc384_apw_list[[unique(gc384_apw$plate)]] <- gc384_apw

#####################################
######### Quality control 1 ###########
#####################################

######### This QC Checks if any well has missing days
#First, get total number of analyzed days
total_days <- length(unique(gc384_apw$day))
# Create a frequency table oer well
well_freq <- setNames(as.data.frame(table(gc384_apw$well)), c("well", "day"))
head(well_freq)
#Get wells that have less than total days
miss_wells <- well_freq[well_freq$day < total_days,]
head(miss_wells)

#Range of missing photos
range(design[design$WELL %in% miss_wells$well,]$PHOTO)

######## This QC checks for outliers

#Loop to plot well by well
outliers <- NULL
for(w in unique(gc384_apw$well)){
  subset_df <- subset(gc384_apw, well == w)
  area <- boxplot.stats(subset_df$area)$out
  #Conditional to get only subsets that have outliers
  if (length(area) > 0) {
  day <- subset_df$day[subset_df$area %in% area]
    #Avoid days from extremities
  day <- day[day > 0 & day < max(day)]
  if (length(day) > 0) {
  #Get wells respective to these days
  well <- rep(w, length(day))
  df <- cbind(well, day)
  #Create data frame
  outliers <- as.data.frame(rbind(outliers,df))
  }
  }
}

outliers <- unique(outliers$well)

#Range of missing photos
design[design$WELL %in% outliers,]$PHOTO

#Use this range on FIJI to check why these wells are not in the analysis.
#Use bash script rsync_croped.sh to download images from source and macro croped_stack.ijm to load them on FIJI.

######################################################
############## Growth Curve Plate Plot ###############
######################################################

#Set plate plot area
par(mfrow = c(16,21), # vector for c(rows, cols)
    mar = c(4, 4, 2, 2), #vector for numner of lines of margins in subplots
    mai = c(0.03,0,0,0.03), #Vector for margin size in inches in subplots
    oma = c(5, 5, 2, 2), # Vector for outer margins for axis labels in bigger plot
    mgp = c(3,0,0)) #c(axis title, tick labels, tick marks)


# Create storage for growth curves fit parameters
gr_fit <- NULL

#Loop to plot well by well
for(w in unique(gc384_apw$well)){
  subset_df <- subset(gc384_apw, well == w)
  
  #Calculate y-axis limits
  y0 <- round(log(min(gc384_apw$area)*0.9),1)
  ymax <- round(log(max(gc384_apw$area)*1.1),1)
  
  df2 <- subset_df[1,][-c(2,8)]
  
  if(length(subset_df$day) >=6){
  #Use function from growthrates package to calculate growth parameters using a specific window (h) of 6 consecutive days for each well
  fit <- fit_easylinear(subset_df$day, subset_df$area, h = 5)

  coefs_fit <- as.data.frame(as.list(c(coef(fit)[3], fit@rsquared, regression = "exponential")))
  df2 <- cbind(df2, coefs_fit)
  # Add to gr_fit
  gr_fit <- rbind(gr_fit, df2)
  
  #Calculate fit coordinates
  x_range <- subset_df$day[fit@ndx]
  #y0 * exp(mumax * time)
  y_range <-  log(coef(fit)[2] * exp(coef(fit)[3] * x_range))
  #Plot growth curve without fit
  plot(log(subset_df$area) ~ subset_df$day, ylim = c(y0, ymax), yaxt = "n", xaxt = "n", 
       #type = "l", 
       axes = T, xlab = "", ylab = "", pch = 16, col = "black", cex = 1)
  #Add linear regression line to the window used to calculate growth parameters
  lines(x_range, y_range, , col = "red", lwd = 1)
  
   } else if(length(subset_df$day) > 1 & length(subset_df$day) < 6){
     
     fit <- summary(lm(log(subset_df$area) ~ subset_df$day))
     coefs_fit <- data.frame(mumax = fit$coefficients[2,1], 
                             r2 = fit$adj.r.squared,
                             regression = "linear")
     df2 <- cbind(df2, coefs_fit)
     # Add to gr_fit
     gr_fit <- rbind(gr_fit, df2)
     
     #Calculate fit coordinates
     x_range <- subset_df$day
     #y0 * exp(mumax * time)
     y_range <-  fit$coefficients[2,1]* x_range + fit$coefficients[1,1]
     #Plot growth curve without fit
     plot(log(subset_df$area) ~ subset_df$day, ylim = c(y0, ymax), yaxt = "n", xaxt = "n", 
          #type = "l", 
          axes = T, xlab = "", ylab = "", pch = 16, col = "black", cex = 1)
     #Add linear regression line to the window used to calculate growth parameters
     lines(x_range, y_range, , col = "red", lwd = 1)
     }

  
  #Plate plot settings
  
#Conditional function to create y-axis labels only in the leftmost col
  if (w %in% gc384_apw$well[gc384_apw$col == 1]){ 
  axis(2, at = c(y0,10,11.8), labels = c(y0-1,10,12), cex.axis = 1, tick = F, hadj = 1.2 , las = 1)}
  
#Conditional function to create x-axis labels only in the bottom row
  max_day <- max(gc384_apw$day)
  if (w %in% gc384_apw$well[gc384_apw$row == "P"]){ 
    axis(1, at = c(par("usr")[1]+0.5, par("usr")[2]/2, par("usr")[2]*0.9),
         labels = c(0,max_day/2, max_day), cex.axis = 1, tick = F, padj = 0)}
  
#Conditional function to create top labels for plate col names
  if (w %in% gc384_apw$well[gc384_apw$row == "A"]){ 
    mtext(substr(w, 2,3), side = 3, line = 0, cex = 1)}
  
#Conditional function to create side labels for plate row names
  if (w %in% gc384_apw$well[gc384_apw$col == 24]){ 
    mtext(substr(w, 1,1), side = 4, cex.axis = 1, line = 1, las = 1)}
  }

# Reset par() to apply global labels to x and y axis
par(fig = c(0, 1, 0, 1), 
    #oma = c(2, 0, 0, 2), 
    mar = c(6, 4, 4, 2) + 0.1,
    #mai = c(1,1,0,0),
    new = TRUE)

# Empty plot to hold labels
plot(NA, xlim = c(1, 24), ylim = c(2, 17), type = "n", axes = F, xlab = "", ylab = "")

# Add y-axis title
mtext(expression("Ln[Area]("*mu*"m²)"), side = 2, line = 6, cex = 1.5, las = 3)

# Add x-axis title
mtext("Time(DAF)", side = 1, line = 8, cex = 1.5, las = 1)


dev.print(jpeg, file = paste0(output, plate, "_gc_all_logfit.jpg"), width = 6000, 
          height = 4000, res = 300, unit = "px") 

###################################
#### Process gr_fit
##################################

# Configure growth parameters data frame
gr_fit <- as.data.frame(gr_fit)
gr_fit$mumax <- as.numeric(gr_fit$mumax)
str(gr_fit)

gr_fit$well <- as.character(gr_fit$well)
str(gr_fit)
head(gr_fit)

##############################
# Create position parameter
################################

#margin_rows <- c("A","B","O", "P")  # First and last row
margin_rows <- c("A","P")  # First and last row
#margin_cols <- c(1,2,23, 24)  # First and last col
margin_cols <- c(1, 24)  # First and last col

position <- ifelse(gr_fit$row %in% margin_rows | gr_fit$col %in% margin_cols, "Edges", "Center")
length(position)

gr_fit$position <- position

#Create colors based on position
colors <- ifelse(gr_fit$position == "Center", "red","blue")


#######################################  Check details per well comparisons ###########


x <- "H9"
dff <- subset(gc384_apw, well == x)

if(length(subset_df$day) >=6){
fit <- fit_easylinear(dff$day, dff$area, h = 5)
#fit <- fit_growthmodel(FUN = grow_gompertz, p = p, dff$day, dff$area)
#fit <- fit_growthmodel(FUN = grow_logistic, p = p, dff$day, dff$area)

plot(fit)
#coef(fit)

} else {

  fit <- summary(lm(log(dff$area) ~ dff$day))
  coefs_fit <- data.frame(mumax = fit$coefficients[2,1], 
                          r2 = fit$adj.r.squared,
                          regression = "linear")
  #Calculate fit coordinates
  x_range <- dff$day
  #y0 * exp(mumax * time)
  y_range <-  fit$coefficients[2,1]* x_range + fit$coefficients[1,1]
  #Plot growth curve without fit
  plot(log(dff$area) ~ dff$day, ylim = c(y0, ymax), yaxt = "n", xaxt = "n", 
       #type = "l", 
       axes = T, xlab = "", ylab = "", pch = 16, col = "black", cex = 1)
  #Add linear regression line to the window used to calculate growth parameters
  lines(x_range, y_range, , col = "red", lwd = 1)
    
}

dev.print(jpeg, file = paste0(output, "/well_",x, "_growthcurve.jpg"), width = 1000, 
          height = 1000, res = 300, unit = "px") 

# Plot several wells
wells <- c("O5", "P16", "F15", "P9")
dff3 <- NULL
for(w in wells){

  dff2 <- subset(gc384_apw, well == w)
  dff3 <- rbind(dff3, dff2)
}

# Assign colors to treatments

colors <- rainbow(length(wells))
colors_match <- colors[as.numeric(factor(dff3$well, levels = wells))]

plot(dff3$area ~ dff3$day, col = colors_match, main = "Well comparisons", xlab = "Time(day)", ylab = expression("Ln[Area]("*mu*"m²)"), 
     pch = 19, 
     #xlim = c(0,25), ylim = c(log(min(df$area))*0.9, log(max(df$area))*1.1), 
     cex = 1)
legend("topleft", legend = unique(dff3$well), col = colors, pch = 16)

dev.print(jpeg, file = paste0(output, "wells_", paste(wells, collapse = "_"), "_gc.jpg"), width = 5000, 
          height = 2500, res = 300, unit = "px") 


####################################################
#### Create list with all data from different plates
####################################################

if (!exists("gr_fit_list")) {
  gr_fit_list <- list()
}

gr_fit_list[[unique(gr_fit$plate)]] <- gr_fit

}


end_time <- Sys.time()
total_time <- end_time - start_time
print(total_time)
names(gr_fit_list)


###########################################################################
########## Analysis of variance to check which sex grows faster ###########
##########################################################################

#First create a new df with only mumax and sexes
male_female <- rbind(gr_fit_list[["plate56"]], gr_fit_list[["plate58"]])
#Add positions to compare too

male_female$sex <- factor(male_female$sex, levels = c("female", "male"))
male_female$position <- factor(male_female$position, levels = c("Center", "Edges"))
male_female$treatment <- factor(male_female$treatment, levels = "nomat", "mat")

head(male_female)
str(male_female)
male_female$mumax <- as.numeric(male_female$mumax)

# Shapiro normality test
shapiro.test(male_female$mumax)
hist(male_female$mumax)

wt <- wilcox.test(male_female$mumax ~ male_female$sex)
format(wt$p.value, scientific = TRUE, digits = 22)
library(summarytools)
 stby(data = male_female$mumax, INDICES = male_female$sex, FUN = descr)

#######################################
#############################################
# Compare edge effects on males and females
##########################################
###########################################

male_female$group <-  interaction(male_female$sex, male_female$position,  sep = ":", lex.order = T)

male_female <- male_female[order(male_female$group), ]
head(male_female)

############################
# Make box plot 
#########################

cols <- c("red", "blue", "red", "blue")
par(mar=c(5, 10, 3, 1))
bp <- boxplot(mumax ~ group, data = male_female,
              boxwex = 0.4, 
              ylim = c(min(male_female$mumax)*0.8, max(male_female$mumax)*1.2),
              xlab = "", ylab = "",  outline = F,
              xaxt = "n", yaxt = "n", 
              #main = "Growth Rate Edge Effects",
              col = cols)

tab_group <- as.numeric(table(male_female$group))
xpos <- rep(which(levels(male_female$group) %in% bp$names), times = tab_group)

points(jitter(xpos, amount = 0.2), male_female$mumax, col = rgb(0, 0, 0, alpha = 0.2), pch = 16, cex = 1)

mtext(expression("Specific Growth Rate ("*day^{-1}*")"), side = 2, line = 7, cex = 2)

# To add the label of x axis 
yscale <-  seq(par("yaxp")[1], par("yaxp")[2], length.out = par("yaxp")[3]+1)

#If removing mat80
axis(1, at = c(1.5, 3.5), labels = c("Females", "Males") , tick = FALSE , cex.axis = 4, padj = 1)
axis(2, at = yscale, labels = yscale, las = 1, hadj = 1.5, cex.axis = 2)

# Divide panels within plt area
ylims <- par("usr")[3:4]  # returns c(ymin, ymax)
# Vertical gray line
segments(y0 = ylims[1], y1 = ylims[2], x0 = 2.5, col = "gray")

library(dunn.test)
#Perform Dunn Test
dt <- dunn.test(male_female$mumax, male_female$group, method = "BH")
pvals <- dt$P.adjusted
names <- gsub(" ", "", dt$comparisons)
names(pvals) <- names
pvals

dt_sex <- dunn.test(male_female_center$mumax, male_female_center$sex)
pvals_sex <- dt_sex$P.adjusted
names <- gsub(" ", "", dt_sex$comparisons)
names(pvals_sex) <- names
pvals_sex

#Get stats letters and export pvalues data
library(multcompView)
letters <- multcompLetters(pvals)
letters <- data.frame(letters = as.character(letters$Letters),
                      treatment = rownames(letters$LetterMatrix))
letters$treatment <- factor(letters$treatment, bp$names)
letters <- letters[order(letters$treatment),]

pvals <- as.data.frame(pvals)
write.table(pvals, file = paste0("pvals_SGR_MxF.csv"), sep = ",")

#add stats letters
text(c(1:ncol(bp$stats)), bp$stats[5,]*1.1, labels = letters$letters, cex = 3)

legend(x = 3.5, y = 0.23, legend = c("Center", "Edges"), 
       col=cols[c(1,2)], y.intersp = 0.4,
       pch = 15, bty = "n", pt.cex = 3, cex = 3,  horiz = F, inset = c(-0.02, 0))

dev.print(jpeg, file = paste0("boxplot_SGR_FxM.jpg"), width = 3000, 
          height = 3000, res = 300, unit = "px") 

#######################################
## Descriptive stats ##############3
#################################
library(summarytools);

desc_stats <- NULL
for(v in c("mumax")){
  df1 <- stby(data = male_female$mumax, INDICES = male_female$group, FUN = descr)
  for(g in unique(male_female$group)){
    df2 <- c(g,unlist(df1[g])[c(1,2,10)])
    desc_stats <- rbind(desc_stats,df2)
  }
}
colnames(desc_stats) <- c("group","mean", "sd", "var")
desc_stats <- as.data.frame(desc_stats)

#write.table(desc_stats, file = paste0("desc_stats_mumax_sexcols.csv"), sep = ",", row.names = F)
write.table(desc_stats, file = paste0("desc_stats_mumax_edeffect.csv"), sep = ",", row.names = F)

#####################################
############ Correlation Growth rate and initial area
#####################################

setwd("~/Dropbox/1_posdoc/1_CCMAR/experiments/384isol/growthrate")

par(mfrow = c(2,2),
    mar = c(0,0.5,0,0),
    oma = c(1,2,2,3)
    
)

for(plt3 in c(56,58)){
  
  platedf3 <- paste0("plate", plt3)
  df3 <- gc384_apw_list[[platedf3]]
  df3_day0 <- subset(df3, day == 0)
  df3_mumax <- gr_fit_list[[platedf3]]
  
model <- lm(df3_mumax$mumax ~ df3_day0$area)
sum_model <- summary(model)

colors <- ifelse(df3_mumax$position == "Edges", "red","blue")

# Get the plot limits
x_lim <- range(df3_day0$area)  # x-axis limits based on the data
y_lim <- range(df3_mumax$mumax)  # y-axis limits based on the data

par(mar = c(4,4,4,6), xpd = T)
plot(df3_mumax$mumax ~ df3_day0$area, type = "n", 
     main = paste0("SDGR x Area(Day0) in ", unique(df3_mumax$sex)), 
     xlab = expression("Area Day0("*mu*"m²)"),
     ylab = expression("SDGR(day"^{-1}*")"),
     xlim = c(x_lim[1]*0.8, x_lim[2*1.2]),
ylim = c(y_lim[1]*0.8, y_lim[2]*1.2)
)

text(df3_day0$area, df3_mumax$mumax, labels = df3_mumax$well,
     col = col[factor(df3_mumax$position, levels = c("Edges", "Center"))], cex = 0.6)

# Calculate y values for the x-axis limits
x_range <- par("usr")[1:2]  # Current x-axis limits
y_range <- coef(model)[1] + coef(model)[2] * x_range  # y = intercept + slope * x

# Ensure the line stays within the plotting area
lines(x_range, y_range, col = "red", lwd = 2)

xlim <- par("usr")[c(1,2)]
ylim <- par("usr")[c(3,4)]

text(xlim[2]*0.85, ylim[2]*0.98, paste0("R² =", round(sum_model$r.squared,4)), cex = 1.5)
pvalue <- signif(sum_model$coefficients[2,4],4)
text(xlim[2]*0.85, ylim[2]*0.95, paste0("pvalue =", pvalue), cex = 1.5)

############ Calculate evaporation correlation ##################

evap_rate <- read.table(file = "../evaporation/evaprate_ser.csv", sep = ",", header = T)
class(evap_rate)
evap_rate_nomat <- subset(evap_rate, treatment == "nomat")    

evap_gr <- merge(evap_rate_nomat, df3_mumax, by = "well")    
str(evap_gr)
plot(evap_gr$mumax ~ evap_gr$rate, ylim = c(y_lim[1]*0.8, y_lim[2]*1.2),
     main = paste0("SGR x SER in ", unique(df3_mumax$sex)), 
     xlab = expression("SDER(day"^{-1}*")"), 
     ylab = expression("SDGR(day"^{-1}*")"), 
     pch = 16, cex = 0.7, col = colors)


model <- lm(df3_mumax$mumax ~ evap_rate_nomat$rate)
sum_model <- summary(model)

# Get the plot limits
x_limits <- range(evap_rate_nomat$rate)  # x-axis limits based on the data
y_limits <- range(df3_mumax$mumax)  # y-axis limits based on the data

# Calculate y values for the x-axis limits
x_range <- par("usr")[1:2]  # Current x-axis limits
y_range <- coef(model)[1] + coef(model)[2] * x_range  # y = intercept + slope * x

# Ensure the line stays within the plotting area
lines(x_range, y_range, col = "red", lwd = 2)

xlim <- par("usr")[c(1,2)]
ylim <- par("usr")[c(3,4)]

text(xlim[2]*0.9, ylim[2]*0.95, paste0("R² =", round(sum_model$r.squared,4)), cex = 1.5)
pvalue <- signif(sum_model$coefficients[2,4],4)
text(xlim[2]*0.9, ylim[2]*0.90, paste0("pvalue =", pvalue), cex = 1.5)

}

# Reset par() to apply global labels to x and y axis
par(fig = c(0, 1, 0, 1), 
    #oma = c(2, 0, 0, 2), 
    mar = c(6, 4, 4, 2) + 0.1,
    #mai = c(1,1,0,0),
    new = TRUE)
plot(NA, xlim = c(0, 5), ylim = c(0, 5), type = "n", axes = F, xlab = "", ylab = "")

legend(4.9, 5.7, legend = c("Edges", "Center"), col = col, pch = 16, bty = "n", inset = c(0,0), pt.cex = 1.5)

text(c(-0.38,2.5), c(5.4,5.4,2.5,2.5), c("A", "B", "C", "D"), cex = 4)

dev.print(jpeg, file = "Correlations_A0_SDERs.jpg", width = 5000, 
          height = 3500, res = 300, unit = "px") 

###########################################################
########## BOX plot per row or per col #################
###########################################################

setwd("~/Dropbox/1_posdoc/1_CCMAR/experiments/384isol/growthrate")

#Set graphical parameters
par(mfrow = c(2,2),
    mar = c(1,1,3,0),
    oma = c(1,1,1,1)
    )

#Loop to make boxplot
for(plt3 in c(56,58)){

  platedf3 <- paste0("plate", plt3)
  df3_mumax <- gr_fit_list[[platedf3]]
  
##########################################
#Growthrate by cols
#######################################

  #First, compare column groups 

library(dunn.test)
#Perform Dunn Test
dt_col <- dunn.test(df3_mumax$mumax, df3_mumax$col, method = "BH")
pvals <- dt_col$P.adjusted
names <- gsub(" ", "", dt_col$comparisons)
names(pvals) <- names
pvals_df <- as.data.frame(pvals)
write.table(pvals_df, file = paste0(platedf3, "_pvals_SGR_percol.csv"), sep = ",")

#Get stats letters and export pvalues data
library(multcompView)
letters_col <- multcompLetters(pvals, threshold = 0.05,)
letters_col <- data.frame(letters = as.character(letters_col$Letters),
                      col = rownames(letters_col$LetterMatrix))
letters_col$col <- as.numeric(letters_col$col)
letters_col <- letters_col[order(letters_col$col),]


#BOx plot

bp <- boxplot(df3_mumax$mumax ~ df3_mumax$col, main = paste0(unique(df3_mumax$sex), " SDGR. Test:  ", unique(df3_mumax$treatment)) ,
        xlab = "cols", 
        ylab = expression("SGR(day"^{-1}*")"), outline = F, cex.lab = 2,
        ylim = c(-0.1, 0.5))

#add stats letters
text(c(1:ncol(bp$stats)), bp$stats[5,]*1.1, labels = letters_col$letters, cex = 1.5)
colors <- ifelse(df3_mumax$position == "Center", "red","blue")

points(jitter(df3_mumax$col, amount = 0.2), df3_mumax$mumax, col = colors, pch = 16, cex = 1)

#axis(1, at = 1:24, labels = unique(gr_fit$col), cex.axis = 2)
#ylabels <- par("yaxp")
#ylabels = seq(y_par[1], y_par[2], length.out = y_par[3])
#yat <- round(seq(par("usr")[3],par("usr")[4], length.out = 4),2)
#axis(2, yat, labels = ylabels)
}

dev.print(jpeg, file = paste0("gr_percol_perplate_all.jpg"), width = 5000, 
          height = 3000, res = 300, unit = "px") 

###############################################################
#Check if there area column effects on sex differences
####################################################################

male_female$group2 <- interaction(male_female$sex, male_female$col,  sep = ":", lex.order = F)

male_female_center <- male_female
male_female_center <- subset(male_female, position == "Center")
#male_female_center <- subset(male_female_center, (row != "B" | row != "O") & (col != 2 | col != 23))
dt_sexcol <- dunn.test(male_female_center$mumax, male_female_center$group2, method = "BH", kw = F, table = F)
pvals <- dt_sexcol$P.adjusted
names <- gsub(" ", "", dt_sexcol$comparisons)
names(pvals) <- names
pvals_df <- as.data.frame(pvals)
write.table(pvals_df, file = paste0("pvals_SGR_sexcol.csv"), sep = ",")

#Get stats letters and export pvalues data
library(multcompView)
library(dunn.test)
letters_sexcol <- multcompLetters(pvals, threshold = 0.05,)
letters_sexcol <- data.frame(letters = as.character(letters_sexcol$Letters),
                          group2 = rownames(letters_sexcol$LetterMatrix))

color_sexcol <- rep(c("#EFEF00", "#BF7B00" ), times = 24) # "#7F5200" "#BF7B00"  "#DF9000" "#7F7F00" "#BFBF00" "#DFDF00" "#EFEF00"
par(mar = c(7, 12, 4.1, 3.1))

bp_sexcol <- boxplot(male_female_center$mumax ~ male_female_center$group2, #main = "Sex-Biased SDGR", 
                     xlab = "", 
                     ylab = "", outline = F, cex.lab = 2, xaxt = "n", yaxt = "n", col = color_sexcol,
                     ylim = c(0, 0.2))

ylabels <- par("yaxp")
ylabels = seq(ylabels[1], ylabels[2], length.out = ylabels[3]+1)

axis(1, at = seq(1.5,ncol(bp_sexcol$stats), by = 2), labels = 1:24, cex.axis = 2)
mtext("Columns", side = 1, line = 4, cex = 3, las = 1)

axis(2, ylabels, labels = ylabels, las = 1, cex.axis = 2, hadj = 1.5)
mtext(expression("SGR(day"^{-1}*")"), side = 2, line = 7, cex = 4, las = 3)

letters_sexcol$group2 <- factor(letters_sexcol$group2, levels = bp_sexcol$names)
letters_sexcol <- letters_sexcol[order(letters_sexcol$group2),]
#add stats letters
text(c(1:ncol(bp_sexcol$stats)), bp_sexcol$stats[5,]*1.05, labels = letters_sexcol$letters, cex = 3.5)

legend(x=par("usr")[2]*0.78, y = par("usr")[4]*1.05, col = color_sexcol, legend = c("Females", "Males"), pch = 15, 
       bty = "n", cex = 3, y.intersp = 0.4, x.intersp = 0.2)


dev.print(jpeg, file = paste0(output, "gr_sexcol_noedges2.jpg"), width = 8000, 
          height = 3000, res = 300, unit = "px") 

########################################
#Initial area by cols
##########################################

boxplot(df3_mumax$y0 ~ df3_mumax$col, main = "Area at Day0 per col (n=16)", xlab = "cols", ylab = expression("Ln[Area]("*mu*"m²)"), outline = F,
        ylim = c(min(df3_mumax$y0)*0.8, max(df3_mumax$y0)*1.1))
points(jitter(df3_mumax$col, amount = 0.2), df3_mumax$y0, col = "red", pch = 16,cex = 0.5)
dev.print(jpeg, file = paste0(output, "/A0_percol_all.jpg"), width = 2500, 
          height = 2000, res = 300, unit = "px") 

dunn_col_area <- dunn.test(df3_mumax$y0 , df3_mumax$col, method = "bonferroni")
dunn_col_pv_area <- dunn_col_area$P.adjusted
dunn_col_area$comparisons[which(dunn_col_pv_area  <= 0.05)]

#################################
#Initial counts by cols
###################################

boxplot(df3_day0$counts ~ df3_day0$col, main = "Counts at Day0 per col (n=16)", xlab = "cols", ylab = "Counts", outline = F,
        ylim = c(min(df3_day0$counts)*0.8, max(df3_day0$counts)*1.1))
points(jitter(df3_day0$col, amount = 0.2), df3_day0$counts, col = "red", pch = 16,cex = 0.5)
dev.print(jpeg, file = paste0(output, "/A0_percol_all.jpg"), width = 2500, 
          height = 2000, res = 300, unit = "px") 

dunn_col_area <- dunn.test(df3_day0$counts , df3_day0$col, method = "bonferroni")
dunn_col_pv_area <- dunn_col_area$P.adjusted
dunn_col_area$comparisons[which(dunn_col_pv_area  <= 0.05)]

# Growth rate by rows
boxplot(gr_fit$mumax ~ gr_fit$row, main = "Growth rate per row (n=24)", xlab = "Rows", ylab = "Specific Growth Rate", outline = F, 
        ylim = c(min(gr_fit$mumax)*0.8, max(gr_fit$mumax)*1.1))
points(jitter(as.numeric(as.factor(gr_fit$row)), amount = 0.2), gr_fit$mumax, col = "red", pch = 16,cex = 0.5)
dev.print(jpeg, file = paste0(output, "/gr_perrow_all.jpg"), width = 2500, 
          height = 2000, res = 300, unit = "px") 

dunn_row <- dunn.test(gr_fit$mumax , gr_fit$row, method = "bonferroni")
dunn_row_pv <- dunn_row$P.adjusted
dunn_row$comparisons[which(dunn_row_pv  <= 0.05)]



#remove 1st layers
gr_fit_noedge <- subset(gr_fit, (row != "A" & row != "P") & (col != 1 & col != 24))
boxplot(gr_fit_noedge$mumax ~ gr_fit_noedge$row, main = "Growth rate per row (n=22)", xlab = "Rows", ylab = "Specific Growth Rate")
boxplot(gr_fit_noedge$mumax ~ gr_fit_noedge$col, main = "Growth rate per col (n=14)", xlab = "cols", ylab = "Specific Growth Rate", outline = F)
points(jitter(gr_fit_noedge$col-1, amount = 0.05), gr_fit_noedge$mumax, col = "red", pch = 16,cex = 0.5)

#remove 2nd layers
gr_fit_noedge2 <- subset(gr_fit_noedge, (row != "B" & row != "O") & (col != 2 & col != 23))
boxplot(gr_fit_noedge2$mumax ~ gr_fit_noedge2$row, main = "Growth rate per row (n=20)", xlab = "Rows", ylab = "Specific Growth Rate", 
        ylim = c(min(gr_fit_noedge2$mumax)*0.8, max(gr_fit_noedge2$mumax)*1.1), 
        outline = F)
points(jitter(as.numeric(as.factor(gr_fit_noedge2$row)), amount = 0.05), gr_fit_noedge2$mumax, col = "red", pch = 16,cex = 0.5)

boxplot(gr_fit_noedge2$mumax ~ gr_fit_noedge2$col, main = "Growth rate per col (n=12)", xlab = "cols", ylab = "Specific Growth Rate",
        ylim = c(min(gr_fit_noedge2$mumax)*0.8, max(gr_fit_noedge$mumax)*1.1),
        outline = F)
points(jitter(gr_fit_noedge2$col-2, amount = 0.05), gr_fit_noedge2$mumax, col = "red", pch = 16,cex = 0.5)
dev.print(jpeg, file = paste0(output, "/growthcurves_2layers.jpg"), width = 2500, 
          height = 2000, res = 300, unit = "px") 

kruskal.test(gr_fit_noedge2$mumax ~ gr_fit_noedge2$col)
dunn_2layers <- dunn.test(gr_fit_noedge2$mumax , gr_fit_noedge2$col, method = "bonferroni")
dunn_2layers_pv <- dunn_2layers$P.adjusted
dunn_2layers$comparisons[which(dunn_2layers_pv <= 0.05)]


############################################################
#Compare upper half X lower half #########################
#########################################################

gr_fit_up <- gr_fit[gr_fit$row %in% c(LETTERS[1:8]),]
gr_fit_down <- gr_fit[gr_fit$row %in% c(LETTERS[9:16]),]


par(mfrow = c(2,1))
boxplot(gr_fit_up$mumax ~ gr_fit_up$col, main = "GR/col (n=8) - Upper Half", xlab = "cols", ylab = "Specific Growth Rate", 
        ylim = c(min(gr_fit_up$mumax)*0.8,max(gr_fit_up$mumax)*1.2),
        outline = FALSE)
# Add jittered points
points(jitter(gr_fit_up$col, amount = 0.2), gr_fit_up$mumax, col = "red", pch = 16,cex = 0.5)

boxplot(gr_fit_down$mumax ~ gr_fit_down$col, main = "GR/col (n=8) - Lower Half", xlab = "cols", ylab = "Specific Growth Rate", 
        ylim = c(min(gr_fit_up$mumax)*0.8,max(gr_fit_up$mumax)*1.2),
        outline = FALSE)
# Add jittered points
points(jitter(gr_fit_down$col, amount = 0.2), gr_fit_down$mumax, col = "red", pch = 16, cex = 0.5)

kruskal.test(gr_fit_up$mumax ~ gr_fit_up$col)
dunn_up <- dunn.test(gr_fit_up$mumax , gr_fit_up$col, method = "bonferroni")
dunn_up_pv <- dunn_up$P.adjusted

kruskal.test(gr_fit_down$mumax ~ gr_fit_down$col)
dunn_down <- dunn.test(gr_fit_down$mumax , gr_fit_down$col, method = "bonferroni")
dunn_down_pv <- dunn_down$P.adjusted
dunn_down$comparisons[which(dunn_down_pv <= 0.05)]

dev.print(jpeg, file = paste0(output, "/growthcurve_up_down.jpg"), width = 5000, 
          height = 4000, res = 300, unit = "px") 

############################################################
############ Analyze 48x8 array ###########################
##########################################################
male_female$half <- ifelse(male_female$row %in% c(LETTERS[1:8]), "upper", "lower")
  
male_female$group3 <- interaction(male_female$sex, male_female$half, male_female$col,  sep = ":", lex.order = F)

male_female_center <- subset(male_female, position == "Center")
#male_female_center <- subset(male_female_center, (row != "B" | row != "O") & (col != 2 | col != 23))
dt_sexcol_half <- dunn.test(male_female$mumax, male_female$group3, method = "BH", kw = F, table = F)
pvals <- dt_sexcol_half$P.adjusted
names <- gsub(" ", "", dt_sexcol_half$comparisons)
names(pvals) <- names
pvals_df2 <- as.data.frame(pvals)
write.table(pvals_df2, file = paste0("pvals_SGR_sexcol_halves.csv"), sep = ",")

#Get stats letters and export pvalues data
library(multcompView)
letters_sexcol_half <- multcompLetters(pvals, threshold = 0.05,)
letters_sexcol_half <- data.frame(letters = as.character(letters_sexcol_half$Letters),
                             group3 = rownames(letters_sexcol_half$LetterMatrix))

color_sexcol <- rep(c("#EFEF00", "#BF7B00" ), times = 24) # "#7F5200" "#BF7B00"  "#DF9000" "#7F7F00" "#BFBF00" "#DFDF00" "#EFEF00"
par(mar = c(5.1, 6.1, 4.1, 3.1))

bp_sexcol <- boxplot(male_female_center$mumax ~ male_female_center$group2, main = "Sex-Biased SDGR", xlab = "cols", 
                     ylab = expression("SGR(day"^{-1}*")"), outline = F, cex.lab = 2, xaxt = "n", col = color_sexcol,
                     ylim = c(0.08, 0.16))

axis(1, at = seq(1.5,ncol(bp_sexcol$stats), by = 2), labels = 1:24, cex.axis = 1, hadj = 0.1)

letters_sexcol$group2 <- factor(letters_sexcol$group2, levels = bp_sexcol$names)
letters_sexcol <- letters_sexcol[order(letters_sexcol$group2),]
#add stats letters
text(c(1:ncol(bp_sexcol$stats)), bp_sexcol$stats[5,]*1.05, labels = letters_sexcol$letters, cex = 0.8)

legend(x=par("usr")[2]*0.89, y = par("usr")[4]*1.02, col = color_sexcol, levels(male_female_center$sex), pch = 15, 
       bty = "n", cex = 2, y.intersp = 0.4, x.intersp = 0.2)


dev.print(jpeg, file = paste0("gr_sexcol_noedges.jpg"), width = 8000, 
          height = 3000, res = 300, unit = "px") 


#Compare left side X right side #######################
gr_fit_lh <- gr_fit[gr_fit$col %in% c(1:12),]
gr_fit_rh <- gr_fit[gr_fit$col %in% c(13:24),]

par(mfrow = c(1,2))
boxplot(gr_fit_lh$mumax ~ gr_fit_lh$col, main = "GR/row (n=12) - Left Half", xlab = "cols", ylab = "Specific Growth Rate", 
        ylim = c(min(gr_fit_lh$mumax)*0.9, max(gr_fit_lh$mumax)*1.1),
        outline = FALSE)
# Add jittered points
points(jitter(gr_fit_lh$col, amount = 0.2), gr_fit_lh$mumax, col = "red", pch = 16, cex = 0.5)

boxplot(gr_fit_rh$mumax ~ gr_fit_rh$col, main = "GR/row (n=12) - Right Half", xlab = "cols", ylab = "Specific Growth Rate", 
        ylim = c(min(gr_fit_lh$mumax)*0.9, max(gr_fit_lh$mumax)*1.1),
        outline = FALSE)
# Add jittered points
points(jitter(gr_fit_rh$col-12, amount = 0.2), gr_fit_rh$mumax, col = "red", pch = 16, cex = 0.5)

#Stats on right
kruskal.test(gr_fit_rh$mumax ~ gr_fit_rh$col)
dunn_right <- dunn.test(gr_fit_rh$mumax , gr_fit_rh$col, method = "bonferroni")
dunn_right_pv <- dunn_right$P.adjusted
dunn_right$comparisons[which(dunn_right_pv <= 0.05)]

#Stats on left
kruskal.test(gr_fit_lh$mumax ~ gr_fit_lh$col)
dunn_left <- dunn.test(gr_fit_lh$mumax , gr_fit_lh$col, method = "bonferroni")
dunn_left_pv <- dunn_left$P.adjusted
dunn_left$comparisons[which(dunn_left_pv <= 0.05)]

dev.print(jpeg, file = paste0(output, "/growthcurve_right_left.jpg"), width = 5000, 
          height = 4000, res = 300, unit = "px") 

########## Margins X Edges ################

#Define margin wells (edges) and center wells
margin_rows <- c("A", "P")  # First and last row
margin_cols <- c(1, 24)  # First and last col

gr_fit_bp <- gr_fit %>%
  mutate(
    position = case_when(
      row %in% margin_rows | col %in% margin_cols ~ "Margin",  # Edge wells
      TRUE ~ "Center"  # Everything else
    )
  )

# Perform statistical test (t-test)
t_test_result <- t.test(mumax ~ position, data = gr_fit_bp)
print(t_test_result)

library(ggplot2)

# Plot growth rate distributions
ggplot(gr_fit_bp, aes(x = position, y = mumax, fill = position)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) + ylim(min(gr_fit_bp$mumax)*0.8, max(gr_fit_bp$mumax)*1.2) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(title = "Growth Rate Comparison: Center vs. Margin",
       x = "Well Position", y = "Maximum Growth Rate (µmax)") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal()

dev.print(jpeg, file = "center_vs_margin_growthrate.jpg", width = 2000, 
          height = 1000, res = 300, unit = "px") 

#######################################################
################## Heatmap Plate Plot##################
#######################################################

# Load necessary package
library(dplyr)

# Set up 384-well plate layout with space for the legend at the bottom
par(mar = c(2, 4, 2, 2), oma = c(5, 4, 2, 2))  # Increased bottom margin

# Compute the mean and standard deviation of `mumax`
mu_avg <- mean(gr_fit$mumax, na.rm = TRUE)
mu_sd <- sd(gr_fit$mumax, na.rm = TRUE)

# Define color palette (10 shades from blue to white to red)
color_palette <- colorRampPalette(c("blue", "white", "red"))(10)

# Categorize wells based on deviation from mean
gr_fit$deviation <- (gr_fit$mumax - mu_avg) / mu_sd  # Z-score (standardized difference)

# Assign colors based on deviation range
gr_fit$color <- color_palette[
  cut(gr_fit$deviation, breaks = 10, labels = FALSE, include.lowest = TRUE)
]

# Extract row (A-P) and col (1-24) from well names
gr_fit$row <- as.numeric(factor(substr(gr_fit$well, 1, 1), levels = rev(LETTERS[1:16])))
gr_fit$col <- as.numeric(gsub("[A-Z]+", "", gr_fit$well))

# Set up layout for main plot and legend
layout(matrix(c(1, 2), nrow = 2), heights = c(0.85, 0.15))  # Main plot takes up 85%, legend 15%

# Main plot (384-well plate)
par(mar = c(2, 4, 2, 2))  # Adjust margins for the main plot
plot(NA, xlim = c(0.5, 24.5), ylim = c(0.5, 16.5), type = "n", axes = FALSE, xlab = "", ylab = "")
for (i in 1:nrow(gr_fit)) {
  well <- gr_fit[i, ]
  rect(well$col - 0.5, well$row - 0.5, well$col + 0.5, well$row + 0.5, col = well$color, border = "black")
}
mtext(1:24, side = 3, at = 1:24, line = 0, cex = 1.5)
mtext(rev(LETTERS[1:16]), side = 2, at = 1:16, line = -1, cex = 1.5, las = 1)

# Legend plot
par(mar = c(2, 4, 1, 2))  # Adjust margins for the legend plot
plot(NA, xlim = c(0, 10), ylim = c(0, 0), type = "n", axes = FALSE, xlab = "", ylab = "")
legend_x_positions <- seq(0.5, 9.5, length.out = 10)
for (i in 1:10) {
  rect(legend_x_positions[i] - 0.5, 0, legend_x_positions[i] + 0.5, 1, 
       col = color_palette[i], border = "black")
}
legend_labels <- round(seq(min(gr_fit$deviation, na.rm = TRUE), 
                           max(gr_fit$deviation, na.rm = TRUE), length.out = 10), 2)
mtext(legend_labels, side = 1, at = seq(0.5, 9.5, length.out = 10), line = 0, cex = 2)
mtext("Deviation from Mean (Z-score)", side = 1, line = 3, cex = 1.5, font = 2)

      dev.print(jpeg, file = paste0(plate, "_heatmap_grfit.jpg"), width = 5000, 
                height = 3000, res = 300, unit = "px") 
     
      
 
################## Correlation Evaporation - Growth rate #################################      
  evap_rate <- read.table(file = "../evaporation/evaprate_ser.csv", sep = ",", header = T)
      class(evap_rate)
      evap_rate_nomat <- subset(evap_rate, treatment == "nomat")    
      nrow(evap_rate_nomat)
  
  evap_gr <- merge(evap_rate_nomat, gr_fit, by = "well")    
  head(evap_gr)
  str(evap_gr)
  
  colors <- ifelse(gr_fit_list[[1]]$position == "Center", "red","blue")
            
  plot(evap_gr$rate ~ evap_gr$mumax, #ylim = c(0.09, 0.1), 
       main = paste0("SGR x SER ", platedf3), xlab = "SGR", ylab = "SER", pch = 16, cex = 0.7, col = colors)
     
  
  model <- lm(evap_rate_nomat$rate ~ gr_fit_list[[1]]$mumax)
  sum_model <- summary(model)
  
  # Get the plot limits
  x_limits <- range(gr_fit$mumax)  # x-axis limits based on the data
  y_limits <- range(evap_rate_nomat$rate)  # y-axis limits based on the data
  
  # Calculate y values for the x-axis limits
  x_range <- par("usr")[1:2]  # Current x-axis limits
  y_range <- coef(model)[1] + coef(model)[2] * x_range  # y = intercept + slope * x
  
  # Ensure the line stays within the plotting area
  lines(x_range, y_range, col = "red", lwd = 2)
  
  text(max(gr_fit$mumax)*0.45, max(evap_rate_nomat$rate)*0.9, paste0("R² =", round(sum_model$r.squared,3)))
  pvalue <- signif(sum_model$coefficients[2,4],4)
  text(max(gr_fit$mumax)*0.45, max(evap_rate_nomat$rate)*0.8, paste0("pvalue =", pvalue))
     
     # Add legend
     legend("topright", legend = c("Center", "Edges"), col = c("blue", "red"), pch = 16, bty = "n")
     
     dev.print(jpeg, file = "scatter_rate_edgeff.jpg", width = 5000, 
               height = 3000, res = 300, unit = "px") 
     
     
  