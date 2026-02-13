#   This script is a big loop with seven parts:
#   1 - Create gc384 file with raw particle analysis
#   2 - Create gc384_apw: an aggregation of all gametophytes per well (APW: area per well)
#   3 - Create plate plot: gr_fit is created based on gc384_apw, with growth rates from higher slope of 6 windows 
#   4 - Create plate plot with heat maps
#   5 - Create box plots with growth rates per strain
#   6 - Create tables with statistical value
#   7 - Create tables with survival data
#
###########################################################################
############################################################################

#Set time 0
start_time <- Sys.time()

#Load libraries
required_packages <- c("growthrates", "ggplot2", "dplyr", "combinat", "summarytools", "car",  "cld", "DescTools")
lapply(required_packages, library, character.only = TRUE)

#Loop to merge all tables well by well for one or several plates
for (plt in c(1:20)) { #Define here the plates that you want to analyze. Separated by commas. If using letters use quotes: "A", or "3A".
  
  #Reset main directory for each plate
  setwd("~/Dropbox/1_posdoc/1_CCMAR/experiments/384WP/brigitta")
 
  #Create plate object and set specific working directory (wd)
plate <- paste0("plate", plt)

#Set working directory for each plate
setwd(paste0(getwd(), "/", plate))
output <- paste0(getwd(), "/")

#Set plate design objects
design <- read.table(file = "plate_design.csv", header = T, sep = ",")

# Nested loop to merge all tables per plate
gc384 <- NULL
for (f in list.dirs(recursive = F)) {
date <- unlist(strsplit(f, "/"))[2]
      for (i in list.files(paste0(f, "/segmented/results/"))) {
      # Skip to the next folder if no files are found
      if (length(f) == 0) {
       #cat("Skipping empty folder:", folder, "\n")
      next  # Go to the next iteration of the loop
                }
       area <- read.table(file = paste0(f, "/segmented/results/", i), header = T, sep = ",")
       area <- area[,1:2]
        colnames(area) <- c("gametophyte", "area")
  
        well1 <- as.character(unlist(strsplit(i, "_"))[3])
        
        #Select the columns that are important to plot
  metadata <- design[design$WELL == well1,][,c(1:length(names(design)))]
  metadata2 <- metadata[rep(1, nrow(area)),]
  date1 <- rep(date, times = nrow(area))
  
  #Create df with all gametophytes area in a well
  img <- cbind(date1, metadata2, area)
  
  #Rename columns
   colnames(img) <- c("date", 
                     names(metadata2), # Variable columns. Depends on experiments
                     "gametophyte", "area") #Measurements column. May vary too.
  colnames(img)[colnames(img) %in% c("ROW", "COLUMN", "WELL", "PHOTO")] <- c("row", "col", "well", "photo")
  gc384 <- rbind(gc384, img)
}
}

gc384$date <- as.Date(as.character(gc384$date), format = "%y%m%d")

#Reorder wells sequence
gc384 <- gc384[order(gc384$row, gc384$col),]
gc384$well <- factor(gc384$well, levels = unique(paste0(gc384$row, gc384$col)))

head(gc384)
nrow(gc384)
str(gc384)

#######################################
####### Total Area per well ###########
#######################################

#Aggregate gametophytes areas
#Reset gc384_apw from previous iterations
gc384_apw <- NULL

#Sum (aggregate) all gametophyte areas per well and reorganize table
gc384_apw <- aggregate(area ~ date + photo + well + row + col + #Standard columns
                         strain + sex + plate + temperature, #Variable columns
                       data = gc384, sum)

#Aggregate gametophyte counts 
gc384_apw$counts <- aggregate(gametophyte ~ date + photo + well + row + col + #Standard columns
                                strain + sex + plate + temperature, #variable columns
                              data = gc384, length)$gametophyte

# Create "day" column subtracting from initial date.
day0 <- min(na.omit(gc384_apw$date))
gc384_apw$day <- as.numeric(gc384_apw$date - as.Date(day0))
gc384_apw <- gc384_apw[order(gc384_apw$row, gc384_apw$col),]

# Create list with data frames from all plates
if (!exists("gc384_apw_list")) {
  gc384_apw_list <- list()
}

plate_name <- paste0("plate", as.character(unique(gc384_apw$plate)))
gc384_apw_list[[plate_name]] <- gc384_apw

#####################################
######### Quality controls ##########
#####################################

######### This QC Checks if any well has missing time points

  #First, get total number of analyzed time points
  total_days <- length(unique(gc384_apw$day))
  
# Create a frequency table per well
  well_freq <- setNames(as.data.frame(table(gc384_apw$well)), c("well", "time.points"))
  well_freq$plate <- rep(plate, times = nrow(well_freq))
  
  #Get wells that have less than total days
  miss_wells <- well_freq[well_freq$time.points < total_days,]
  #head(miss_wells)
  
  if (nrow(miss_wells) > 0) {
    
  if (!exists("miss_wells_list", inherits = F)) {
    miss_wells_list <- list()
  }
  
  plate_miss <- as.character(unique(miss_wells$plate))
  miss_wells_list[[plate_miss]] <- miss_wells
  
  }
  
#Use bash script rsync_croped.sh to download images from source and macro croped_stack.ijm to load them on FIJI.

######################################################
############## Growth Curve Plate Plot ###############
######################################################

#Add empty wells to main data frame. This is needed as the two outer rows and columns are excluded to minimize edge effects.

#Create hypothetical object with all possible wells
rows <- LETTERS[1:16]
columns <- 1:24
wells <- as.vector(outer(rows, columns, paste0))

#Check which wells are empty in the current plate analyzed
empty_wells <- setdiff(wells, unique(gc384_apw$well))

#Create empty df with same format of gc384_apw
empty_wells_df <- data.frame(matrix(ncol = ncol(gc384_apw), nrow = 0))
colnames(empty_wells_df) <- colnames(gc384_apw)

#Fill empty df with wells absent from current plate
for (w in 1:length(empty_wells)) {
  empty_wells_df[w,3] <- empty_wells[w]
}

#Adjust column properties
empty_wells_df$col <- as.numeric(substr(empty_wells_df$well,2,3))
empty_wells_df$row <- factor(substr(empty_wells_df$well,1,1), levels = LETTERS[1:16])
empty_wells_df$day <- 0
empty_wells_df$area <- NA

# Merge empty well df with gc384_apw
gc384_apw <- rbind(gc384_apw, empty_wells_df)
gc384_apw$row <- factor(gc384_apw$row, levels = LETTERS[1:16])

#tail(gc384_apw$row)

#Reorder df order based in row and column order and reclassify well column as factor
gc384_apw <- gc384_apw[order(gc384_apw$row, gc384_apw$col),]
#tail(gc384_apw)
gc384_apw$well <- factor(gc384_apw$well, levels = unique(gc384_apw$well))

# Growth curves plate plot
#Set plate plot area
par(mfrow = c(16,24), # vector for c(rows, cols)
    mar = c(4, 4, 2, 2), #vector for number of lines of margins in subplots
    mai = c(0.03,0,0,0.03), #Vector for margin size in inches in subplots
    oma = c(5, 5, 2, 2), # Vector for outer margins for axis labels in bigger plot
    mgp = c(3,0,0)) #c(axis title, tick labels, tick marks)

# Create storage for growth curves fit parameters
gr_fit_simple <- NULL

#Loop to plot well by well
for (w in unique(gc384_apw$well)) {
  #Calculate y-axis limits
  y0 <- round(log(min(gc384_apw$area + 1, na.rm = T)*0.9),1)
  ymax <- round(log(max(gc384_apw$area, na.rm = T)*1.1),1)
  
    #Create subset for growth curve plot and growth rate calculations  
    subset_df <- subset(gc384_apw, well == w)
    #Subset without area, counts and day for the gr_fit object
    df2 <- subset_df[1,1:(ncol(subset_df) - 3)]
    
    if (length(subset_df$day) > 2) {
      
    fit <- summary(lm(log(subset_df$area) ~ subset_df$day))
    coefs_fit <- data.frame(mumax = fit$coefficients[2,1], 
                            r2 = fit$adj.r.squared,
                            regression = "linear")
    
    df2 <- cbind(df2, coefs_fit)
    
    # Add to gr_fit
    gr_fit_simple <- rbind(df2, gr_fit_simple)
    }
    
    plot(log(subset_df$area) ~ subset_df$day, ylim = c(y0, ymax), yaxt = "n", xaxt = "n", 
         axes = T, xlab = "", ylab = "", pch = 16, col = "black", cex = 1)

    #Conditional function to create y-axis labels only in the leftmost column
    if (w %in% gc384_apw$well[gc384_apw$col == 1]) { 
      axis(2, at = c(par("usr")[3]*1.3, mean(c(par("usr")[3], par("usr")[4])), par("usr")[4])*0.9,
           labels = c(y0 - 1,ymax/2,ymax), cex.axis = 1, tick = F, hadj = 1.2 , las = 1)
    }
    
    #Conditional function to create x-axis labels only in the bottom row
    max_day <- max(gc384_apw$day, na.rm = T)
    if (w %in% gc384_apw$well[gc384_apw$row == "P"]) { 
      axis(1, at = c(par("usr")[1]*0.97, mean(c(par("usr")[1],par("usr")[2])), par("usr")[2])*0.85,
           labels = c(0, max_day/2, max_day), cex.axis = 1, tick = F, padj = 0)
    }
    
    #Conditional function to create top labels for plate col names
    if (w %in% gc384_apw$well[gc384_apw$row == "A"]) { 
      mtext(substr(w, 2,3), side = 3, line = 0, cex = 1)
    }
    
    #Conditional function to create side labels for plate row names
    if (w %in% gc384_apw$well[gc384_apw$col == 24]) { 
      mtext(substr(w, 1,1), side = 4, cex.axis = 1, line = 1, las = 1)
    }
}
    # Reset par() to apply global labels to x and y axis to whole plate plot
    par(fig = c(0, 1, 0, 1), 
        mar = c(6, 4, 4, 2) + 0.1,
        new = TRUE)
    
    # Empty plot to hold labels
    plot(NA, xlim = c(1, 24), ylim = c(2, 17), type = "n", axes = F, xlab = "", ylab = "")
    
    # Add y-axis title
    mtext(expression("Ln[Area("*mu*"m²)]"), side = 2, line = 6, cex = 1.5, las = 3)
    
    # Add x-axis title
    mtext("Time(DAF)", side = 1, line = 8, cex = 1.5, las = 1)
    
    #Save plate plot in its specific plate and correct identification
    dev.print(jpeg, file = paste0(output, plate, "_gc_all_nofit.jpg"), width = 6000, 
              height = 4000, res = 300, unit = "px") 

    #Wells to check. Wells not empty
    wne <- setdiff(gc384_apw$well, empty_wells)
    #Get wells that are not empty where growth rate couldn't be calculated
    wtc <- setdiff(wne, gr_fit_simple$well)
    
    warning(
      paste0("Please Check areas in wells: ", paste(wtc, collapse = ", "), " from plate: ", plt)
    )
     
    
    # Growth curves plate plot with regression
    #Set plate plot area
    par(mfrow = c(16,24), # vector for c(rows, cols)
        mar = c(4, 4, 2, 2), #vector for number of lines of margins in subplots
        mai = c(0.03,0,0,0.03), #Vector for margin size in inches in subplots
        oma = c(5, 5, 2, 2), # Vector for outer margins for axis labels in bigger plot
        mgp = c(3,0,0)) #c(axis title, tick labels, tick marks)
    
# Recalculate growth rates and plot the regression lines considering exponential phase of growth and regressions
# Create storage for growth curves fit parameters
gr_fit <- NULL

#Loop to plot well by well
for (w in unique(gc384_apw$well)) {
  #Calculate y-axis limits
  y0 <- round(log(min(gc384_apw$area + 1, na.rm = T)*0.9), 1)
  ymax <- round(log(max(gc384_apw$area, na.rm = T)*1.1), 1)
  
  #Create subset for growth curve plot and growth rate calculations  
  subset_df <- subset(gc384_apw, well == w)
  #Subset without area, counts and day for the gr_fit object
  df2 <- subset_df[1,1:(ncol(subset_df) - 3)]
    
     #Conditional plot for empty wells
    if (sum(is.na(subset_df$area)) >= 1) {
    
    # Empty placeholder plot for excluded wells
    range_x <- range(gc384_apw$day, na.rm = T)
    range_y <- round(log(range(gc384_apw$area, na.rm = TRUE)))
    
    plot(NA, xlim = range_x, ylim = range_y,
         yaxt = "n", xaxt = "n", 
         axes = T, type = "n",
         xlab = "", ylab = "")
    
    ##################################################################################################################
    # Conditional plot if the wells has only one time point. Necessary because its not possible to calculate growth rate
  } else if (length(subset_df$day) == 1) {
    plot(log(subset_df$area) ~ subset_df$day, ylim = c(y0, ymax), yaxt = "n", xaxt = "n", 
         axes = T, xlab = "", ylab = "", pch = 16, col = "black", cex = 1)
    
    ##################################################################################################
    # Conditional plot if wells have enough days. This can be changed accordingly.
   
  } else if (length(subset_df$day) > 1 & length(subset_df$day) < 5) {
  # First, calculate growth rate considering all time points. The first goal is to check whether growth is positive or negative.
    #This will influence the next steps
    
    fit <- summary(lm(log(subset_df$area) ~ subset_df$day))
    coefs_fit <- data.frame(mumax = fit$coefficients[2,1], 
                            r2 = fit$adj.r.squared,
                            regression = "linear")
    df2 <- cbind(df2, coefs_fit)
    
    # Add to gr_fit
    gr_fit <- rbind(df2, gr_fit)
    
    #Calculate fit coordinates
    x_range <- subset_df$day
    y_range <-  fit$coefficients[2,1] * x_range + fit$coefficients[1,1]
    
    #Plot growth curve without fit
    plot(log(subset_df$area) ~ subset_df$day, ylim = c(y0, ymax), yaxt = "n", xaxt = "n", 
         axes = T, xlab = "", ylab = "", pch = 16, col = "black", cex = 1)
    
    #Add linear regression line to the window used to calculate growth parameters
    lines(x_range, y_range, , col = "red", lwd = 1)
    
    #Plot to differentiate increasing or decreasing growth. If increasing, use fit_easylinear()
  } else if (length(subset_df$day) > 5) {
    fit <- summary(lm(log(subset_df$area) ~ subset_df$day))
    pre_coeffs <- data.frame(mumax = fit$coefficients[2,1], 
                            r2 = fit$adj.r.squared,
                            regression = "linear")
    
    if (pre_coeffs$mumax > 0) {
    
    ################################################
  #Use function from growthrates package to calculate growth parameters using a specific window (h) of consecutive days for each well
    # Here h should be changed accordingly to get better comparisons.
    ##################################################
  fit <- fit_easylinear(subset_df$day, subset_df$area, h = 5)
  
  #Create object with growth parameters (mumax (growth rate), R² and the type of regression))
  coefs_fit <- as.data.frame(as.list(c(coef(fit)[3], fit@rsquared, regression = "exponential")))
  df2 <- cbind(df2, coefs_fit)
  
  # Add to gr_fit
  gr_fit <- rbind(df2, gr_fit)
  
  #Calculate fit coordinates
  x_range <- subset_df$day[fit@ndx]
  y_range <-  log(coef(fit)[2] * exp(coef(fit)[3] * x_range))
  
  #Plot growth curve without fit
  plot(log(subset_df$area) ~ subset_df$day, ylim = c(y0, ymax), yaxt = "n", xaxt = "n", 
       axes = T, xlab = "", ylab = "", pch = 16, col = "black", cex = 1)
  
  #Add linear regression line to the window used to calculate growth parameters
  lines(x_range, y_range, , col = "red", lwd = 1)
    }
  ################################
  # Linear plot
  
   if (pre_coeffs$mumax < 0) {
     
     fit <- summary(lm(log(subset_df$area) ~ subset_df$day))
     coefs_fit <- data.frame(mumax = fit$coefficients[2,1], 
                             r2 = fit$adj.r.squared,
                             regression = "linear")
     df2 <- cbind(df2, coefs_fit)
     
     # Add to gr_fit
     gr_fit <- rbind(df2, gr_fit)
     
     #Calculate fit coordinates
     x_range <- subset_df$day
     y_range <-  fit$coefficients[2,1] * x_range + fit$coefficients[1,1]
     
     #Plot growth curve without fit
     plot(log(subset_df$area) ~ subset_df$day, ylim = c(y0, ymax), yaxt = "n", xaxt = "n", 
          axes = T, xlab = "", ylab = "", pch = 16, col = "black", cex = 1)
     
     #Add linear regression line to the window used to calculate growth parameters
     lines(x_range, y_range, , col = "red", lwd = 1)
     
     }

  }
  #Plate plot settings
  
  #Conditional function to create y-axis labels only in the leftmost column
  if (w %in% gc384_apw$well[gc384_apw$col == 1]) { 
    axis(2, at = c(par("usr")[3]*1.3, mean(c(par("usr")[3], par("usr")[4])), par("usr")[4])*0.9,
         labels = c(y0 - 1,ymax/2,ymax), cex.axis = 1, tick = F, hadj = 1.2 , las = 1)
  }
  
  #Conditional function to create x-axis labels only in the bottom row
  max_day <- max(gc384_apw$day, na.rm = T)
  if (w %in% gc384_apw$well[gc384_apw$row == "P"]) { 
    axis(1, at = c(par("usr")[1]*0.97, mean(c(par("usr")[1],par("usr")[2])), par("usr")[2])*0.85,
         labels = c(0, max_day/2, max_day), cex.axis = 1, tick = F, padj = 0)
  }
  
#Conditional function to create top labels for plate col names
  if (w %in% gc384_apw$well[gc384_apw$row == "A"]) { 
    mtext(substr(w, 2,3), side = 3, line = 0, cex = 1)
    }
  
#Conditional function to create side labels for plate row names
  if (w %in% gc384_apw$well[gc384_apw$col == 24]) { 
    mtext(substr(w, 1,1), side = 4, cex.axis = 1, line = 1, las = 1)
  }
  
  } # End of plate plot loop

# Reset par() to apply global labels to x and y axis to whole plate plot
par(fig = c(0, 1, 0, 1), 
    mar = c(6, 4, 4, 2) + 0.1,
    new = TRUE)

# Empty plot to hold labels
plot(NA, xlim = c(1, 24), ylim = c(2, 17), type = "n", axes = F, xlab = "", ylab = "")

# Add y-axis title
mtext(expression("Ln[Area("*mu*"m²)]"), side = 2, line = 6, cex = 1.5, las = 3)

# Add x-axis title
mtext("Time(DAF)", side = 1, line = 8, cex = 1.5, las = 1)

#Save plate plot in its specific plate and correct identification
dev.print(jpeg, file = paste0(output, plate, "_gc_all_logfit.jpg"), width = 6000, 
          height = 4000, res = 300, unit = "px") 

# Set gr_fit as data frame, check its structure and change class of well column
gr_fit <- as.data.frame(gr_fit)
gr_fit$well <- as.character(gr_fit$well)
str(gr_fit)
head(gr_fit)


if (!exists("gr_fit_list")) {
  gr_fit_list <- list()
}

plate_gr <- as.character(unique(gr_fit$plate))
gr_fit_list[[plate_gr]] <- gr_fit

#######################################  Check details per well comparisons ###########

graphics.off()

x <- "H12"
dff <- subset(gc384_apw, well == x)

  #Plot growth curve without fit
  plot(dff$area ~ dff$day, #ylim = c(y0, ymax), #yaxt = "n", xaxt = "n", 
       #type = "l", 
       axes = T, pch = 16, col = "black", cex = 1,
       main = paste0("Growth Curve well: ", x),
       xlab = "Time(day)", ylab = expression("Aarea("*mu*"m²)"))


dev.print(jpeg, file = paste0(output, "/well_",x, "_growthcurve.jpg"), width = 1000, 
          height = 1000, res = 300, unit = "px") 

# Plot several wells

par(xpd = T,
    mar = c(5.1, 4.1, 4.1, 4.1))

wells <- c("H13", "H14", "I13", "I14")
dff3 <- NULL

for (w in wells) {
  dff2 <- subset(gc384_apw, well == w)
  dff3 <- rbind(dff3, dff2)
}

# Assign colors to treatments

colors <- rainbow(length(wells))
colors_match <- colors[as.numeric(factor(dff3$well, levels = wells))]

plot(dff3$area ~ dff3$day, col = colors_match, main = "Well comparisons", xlab = "Time(day)", ylab = expression("APW("*mu*"m²)"), 
     pch = 19, cex = 1)
legend(x = max(dff3$day) + 1, y = max(dff3$area)*1.1, legend = unique(dff3$well), col = colors, pch = 16, bty = "n", title = "Well")

dev.print(jpeg, file = paste0(output, "wells_", paste(wells, collapse = "_"), "_gc.jpg"), width = 1000, 
          height = 1000, res = 200, unit = "px", ) 

}

####################################################
#### Create list with all data from different plates
####################################################

names(gr_fit_list)
names(gc384_apw_list)

### End of First loop ############
#######################################################################################################################################################################

# In this part you can compare the strains or treatments by choosing the plate to compare

#Which plates do you want to compare?

#Create sequence from 1 to 11 by two
#Odd plates
comp_plates <- as.character(seq(1, 20, 1) )
#Even plates
#comp_plates <- as.character(seq(2, 12, 2))

#Create an empty object to assemble all data in the end
comparison <- NULL

#Loop to create a new df with only the plates you want to compare
for (p in comp_plates) {
comparison <- rbind(comparison, gr_fit_list[[p]])
}

unique(comparison$plate)
head(comparison)

# Set variables as factors for proper ordering on plots

comparison$mumax <- as.numeric(comparison$mumax)
comparison$temperature <- factor(comparison$temperature)
#Remove NAs
comparison <- comparison[!(is.na(comparison$mumax)),]
head(comparison)
nrow(comparison)
str(comparison)
comparison$r2 <- as.numeric(comparison$r2)
sum(is.na(comparison$r2))

#Plot raw data per plate
colors <- rainbow(length(unique(comparison$plate)), alpha = 0.8)
colors_match <- colors[comparison$plate]

length(colors_match)
length(comparison$r2)

output2 <- "~/Dropbox/1_posdoc/1_CCMAR/experiments/384WP/brigitta/"

par(mfrow = c(1,2), mai = c(1,1,1,1),
    mar = c(5.1, 4.1, 1, 8), xpd = T)
plot(comparison$r2, col = colors_match, main = "R2 per plate", xlab = "Sample Index", ylab = "R squared")
plot(comparison$mumax, col = colors_match, main = "Umax per plate", xlab = "Sample Index", ylab = "SGR")

legend(par("usr")[2]*0.85, par("usr")[4], legend = paste0("Plate", unique(comparison$plate)), x.intersp = 0.2,
       col = colors, pch = 16, bty = "n", pt.cex = 1.5)


dev.print(jpeg, file = paste0(output2,"platewise_rsquare_mumax_all.jpg"), width = 3000, 
          height = 1500, res = 200, unit = "px") 

#Plot raw data per temperature
colors2 <- c("blue", "red")
colors_temp <- colors2[comparison$temperature]

par(mfrow = c(1,2), mai = c(1,1,1,1),
    mar = c(5.1, 4.1, 1, 8), xpd = T)
plot(comparison$r2, col = colors_temp, main = "R2 per plate", xlab = "Sample Index", ylab = "R squared")
plot(comparison$mumax, col = colors_temp, main = "Umax per plate", xlab = "Sample Index", ylab = "SGR")

legend(par("usr")[2]*0.85, par("usr")[4], legend = paste0(levels(comparison$temperature), "ºC"), x.intersp = 0.2,
       col = colors2, pch = 16, bty = "n", pt.cex = 1.5)

dev.print(jpeg, file = paste0(output2,"tempwise_rsquare_mumax_all.jpg"), width = 3000, 
          height = 1500, res = 200, unit = "px") 


# Check number of replicates

comp_reps <- as.data.frame(table(comparison$strain))

#Strains with more than 24 replicates
comp_reps_m24 <- subset(comp_reps, Freq > 48)

comparison_m24 <- comparison[comparison$strain %in% comp_reps_m24$Var1,]
unique(comparison_m24$plate)

#Filter some data

summary(comparison$mumax)
#Filter by linearity of the regression
comparison <- subset(comparison, r2 >= 0.5)
comparison <- comparison[!(is.na(comparison$r2)),]

# Shapiro normality test. If pvalue is smaller than 0.05, the data is not normally distributed. However, visual inspection of Q–Q plots showed approximate normality
# in the central distribution (red line), with deviations in the tails. Given the large sample size (n > 1200), parametric analyses can be considered appropriate.

par(mfrow = c(1,2))
hist(comparison$mumax, probability = T)
lines(density(comparison$mumax),col = "red")
qqPlot(comparison$mumax)
ks.test(comparison$mumax, "pnorm", mean = mean(comparison$mumax), sd = sd(comparison$mumax))

# Check normality of the variance test
model_mumax <- lm(mumax ~ strain*temperature, data = comparison)
summary(model_mumax)

par(mfrow = c(2,2))
plot(model_mumax)
vif(model_mumax, type = "predictor")
ks.test(residuals(model_mumax), "pnorm", mean = mean(comparison$mumax), sd = sd(comparison$mumax))

cutoff <- 4/((nrow(comparison) - length(model_mumax$coefficients) - 2))

leveneTest(residuals(model_mumax) ~ comparison$strain)
dev.off()
plot(model_mumax, which = 5, cook.levels = cutoff)
plot(model_mumax, which = 4)

outliers <- outlierTest(model_mumax)
class(outliers)

comparison_outliers <- comparison[rownames(comparison) %in% names(outliers$rstudent),]
head(comparison_outliers)
nrow(comparison_outliers)


#Identify outliers

comparison$out_studentized <- rstudent(model_mumax)
comparison$cooks_distance <- cooks.distance(model_mumax)

resid_thresh <- 3
cook_thresh  <- 4 / nrow(comparison)

outlier_table <- comparison %>%
  dplyr::mutate(
    outlier_residual = abs(out_studentized) > resid_thresh,
    outlier_cook     = cooks_distance > cook_thresh
  ) %>%
  dplyr::filter(outlier_residual | outlier_cook) %>%
  dplyr::select(
    strain,
    mumax,
    well,
    plate,
    out_studentized,
    cooks_distance,
    outlier_residual,
    outlier_cook
  )

outlier_table
nrow(outlier_table)

comparison <- comparison[!(rownames(comparison) %in% rownames(outlier_table)),]
nrow(comparison)

#######################################
#############################################
# Plot growth curves per strain and calculate survival rates
##########################################
###########################################

# This creates a new column (variable) to compare multiple parameters
#comparison$group <-  interaction(comparison$treatment, comparison$volume, comparison$strain,  sep = ":", lex.order = T)

#calculate variation of variations

desc_stats_area <- NULL
surv_all <- NULL 

for (p in unique(names(gc384_apw_list))) {
  #Get Initial and final areas. Make growth curve per strain
  #Set plate plot area
  par(mfrow = c(3, 4), # vector for c(rows, cols)
      mar = c(4, 4, 2, 2), #vector for number of lines of margins in subplots
      mai = c(0.03,0,0,0.03), #Vector for margin size in inches in subplots
      oma = c(5, 5, 2, 2), # Vector for outer margins for axis labels in bigger plot
      mgp = c(3,0,0)) #c(axis title, tick labels, tick marks)
  
  df1 <- gc384_apw_list[[p]]
  for (q in unique(df1$strain)) {
  
    #Create object per strain
    df2 <- subset(df1, strain == q)
    
    #Plot growth curves
    
    # Growth curves plate plot with regression
   
    #Define the mininum and maximum of Y-axis
    y0 <- round(log(min(df1$area + 1, na.rm = T) * 0.9), 1)
    ymax <- round(log(max(df1$area, na.rm = T)) * 1.1, 1)
    
    #Plot area per day
    plot(log(df2$area) ~ df2$day, ylim = c(y0, ymax), yaxt = "n", xaxt = "n", 
         axes = T, xlab = "", ylab = "", pch = 16, col = "black", 
         cex = 1)
    
    #Conditional function to create y-axis labels only in the leftmost column
    if (q %in% unique(df1$strain)[c(1,5,9)]) { 
      axis(2, at = c(par("usr")[3]*1.3, mean(c(par("usr")[3], par("usr")[4])), par("usr")[4])*0.9,
           labels = c(y0 - 1, ymax/2, ymax), cex.axis = 1, tick = F, hadj = 1.2 , las = 1)
    }
    
    #Conditional function to create x-axis labels only in the bottom row
    #max_day <- max(df2$day, na.rm = T)
    #if (q %in% unique(df1$strain)[7:10]){ 
    #  axis(1, at = c(par("usr")[1]*0.85, mean(c(par("usr")[1],par("usr")[2])), par("usr")[2])*0.95,
    #       labels = c(0, max_day/2, max_day), cex.axis = 1, tick = F, padj = 0)
    #}
    days <- unique(df2$day)
    max_day <- max(df2$day, na.rm = T)
    if (q %in% unique(df1$strain)[7:10]) { 
      axis(1, at = days,
           labels = days,cex.axis = 1, tick = F, padj = 0)
    }
    
    # Add title
    mtext(q, side = 3, line = -2, cex = 1)
    
    ######################## Survivorship only for plates cultured at 26C ################
    
   
   # if(p %in% c("plate2", "plate4","plate6","plate8", "plate10", "plate12")){
      for (w in unique(df2$well)) {
        
        #Define a new object per well. To calculate the survivor rate
        df3 <- subset(df2, well == w)
        #p <- unique(df3$plate)
       
         if (length(df3$day) > 1) {
          
          #Calculate final area
          f_area <- df3$area[which.max(df3$day)]
          #Calculate initial area
          i_area <- df3$area[which.min(df3$day)]
          
          #Calculate survivor rate ( simple percentage)
        surv1 <- (1 - (i_area - f_area)/i_area
                         ) * 100
                 
        #Create a new data frame with all the specific data from each loop
        surv2 <- data.frame(plate = p,
                           strain = q,
                           temperature = unique(df3$temperature),
                           well = w,
                           day = max(df3$day),
                           survival = surv1)
        
        #Add the new data from the new loop to one data frame
        surv_all <- rbind(surv_all, surv2)
        }
      }
    #}
    
               }
  
  
  # Reset par() to apply global labels to x and y axis to whole plate plot
  par(fig = c(0, 1, 0, 1), 
      mar = c(6, 4, 4, 2) + 0.1,
      new = TRUE)
  
  # Empty plot to hold labels
  plot(NA, xlim = c(1, 4), ylim = c(1, 3), type = "n", axes = F, xlab = "", ylab = "")
  
  # Add title
  mtext(p, side = 3, line = 4.5, cex = 1)
  
  # Add y-axis title
  mtext(expression("Ln[Area("*mu*"m²)]"), side = 2, line = 6, cex = 1.5, las = 3)
  
  # Add x-axis title
  mtext("Time(DAF)", side = 1, line = 8, cex = 1.5, las = 1)
  
  #Save plate plot in its specific plate and correct identification
  dev.print(jpeg, file = paste0(output2, p, "/", p, "_gc_allstrains.jpg"), width = 3000, 
            height = 2000, res = 300, unit = "px") 
  
}
  
######## Survival Rate Plot ##########################################################
  
write.table(surv_all, file = paste0(output2, "survival_rates.csv"), sep = ",", row.names = F)
#First exclude all survival lower than 0% and higher than 100%.

nrow(surv_all)
surv_all$strain <- factor(surv_all$strain)
length(levels(surv_all$strain))

surv_all_filt <- subset(surv_all, survival >= 0 & survival <= 100)
#surv_all_filt <- surv_all
nrow(surv_all_filt)
length(unique(surv_all_filt$strain))

#Second, reorder strains by decreasing survival rate
surv_all_filt$strain <- reorder(surv_all_filt$strain, surv_all_filt$survival, median)
surv_all_filt$day <- as.factor(surv_all_filt$day)
surv_all_filt$temperature <- factor(surv_all_filt$temperature)

length(unique(surv_all_filt$strain))
table(surv_all_filt$strain)

strain_levels_temp <- NULL

for (temp in levels(surv_all_filt$temperature)) {
surv_all_filt_temp <- subset(surv_all_filt, temperature == temp)

strain_levels <- levels(surv_all_filt_temp$strain)
strain_levels_df <- cbind(strain_levels, rep(temp, times = length(strain_levels)))
strain_levels_temp <- rbind(strain_levels_temp, strain_levels_df)

  par(mar = c(6, 10, 3, 8), xpd = T)
  bp <- boxplot(survival ~ strain, data = surv_all_filt_temp,
                boxwex = 0.4, horizontal = T,
                ylim = c(0, 130),
                xlab = "", ylab = "",  outline = F, 
                xaxt = "n", yaxt = "n", 
                #col = cols,
                main = "Survival Rate")
  
  #tab_group <- length(unique(comparison$antibiotic))
  xpos <- match(surv_all_filt_temp$strain, bp$names)
  
  colors <- rainbow(length(unique(surv_all_filt_temp$day)), alpha = 0.8)
  colors_match <- colors[surv_all_filt_temp$day]
  
  points(jitter(surv_all_filt_temp$survival, amount = 0.2), xpos, col = rgb(1, 0, 0, alpha = 0.3), pch = 16, cex = 0.5) #may be wrong
  #points(jitter(surv_all_filt$survival, amount = 0.2), xpos, col = colors_match, pch = 16, cex = 0.5) #may be wrong
  
  mtext("Gametophyte Survival(%)", side = 1, line = 3, cex = 2)
  
  # To add the label of x axis 
  yscale <-  seq(par("xaxp")[1], par("xaxp")[2], length.out = par("xaxp")[3] + 1)
  xaxis <- par("xaxp")
  #n_treatments <- levels(comparison$treatment)
  n_comp <- unique(bp$names)
  #n_vols <- unique(comparison$volume)
  
  #Add xlabels
  axis(2, at = 1:length(bp$names), labels = FALSE, tck = -0.01)
  text(y = 1:length(bp$names), x = -8, labels = n_comp, cex = 1, adj = 1)
  axis(1, at = yscale, labels = yscale, las = 1, cex.axis = 2, padj = 0.8)
  
  #legend(x = par("usr")[2]*0.9, y = par("usr")[4]*1.1, legend = levels(surv_all_filt$day), col = colors, pch = 16, bty = "n", title = "Maximum Day")
  
  ######### Get pvalues and get the Compact Letter Display
  #“All pairwise comparisons among the 53 strains (n = 1378) were conducted with multiple-testing correction,
  #and results were summarized using a compact letter display following Piepho (2004).”
  

  #Calculates all pairwise calculations
  choose(length(unique(surv_all_filt$strain)), 2)
  
  #Perform Dunn Test
   dt <- DunnTest(surv_all_filt$survival, surv_all_filt$strain)
  
  #write.table(dt, file = paste0(output2, "survival_pvalues_", temp, "C.csv"), sep = ",")
  
  #Get letters
  cld_res <- make_cld(dt)
  
  letters <- as.data.frame(cld_res)
  
  letters$group <- factor(letters$group, levels(surv_all_filt$strain))
  letters <- letters[order(letters$group),]
  
    write.table(letters, file = paste0(output2, "survival_letters_", temp, "C.csv"), sep = ",")
  
  #add stats letters
  text(bp$stats[5,] + 1, c(1:ncol(bp$stats)), labels = letters$cld, cex = 1, adj = 0)
  
  
  dev.print(jpeg, file = paste0(output2, "survival_allstrains_", temp, "C.jpg"), width = 2500, 
            height = 4000, res = 300, unit = "px") 
}
  #######################################################################
  #Descriptive statistics for survival rates
  
  
  desc_stats_surv <- NULL
  df5 <- stby(data = surv_all$survival, INDICES = surv_all$strain, FUN = descr)
  
  for (g in unique(surv_all_filt$strain)) {
    df6 <- c(g, unlist(df5[g])[c(1,2,10,15)])
    desc_stats_surv <- rbind(desc_stats_surv, df6)
  }
  
  colnames(desc_stats_surv) <- c("group","mean", "sd", "var", "replicates")
  desc_stats_surv <- as.data.frame(desc_stats_surv)
  
  write.table(desc_stats_surv, file = paste0(output2,"survival_allstrain_stats.csv"), sep = ",", row.names = F)
  
############################################################################################################################
####################### plot growth rates and do the statistics ################################
  ###################################################################################
  gr_fit_all <- do.call(rbind, gr_fit_list)
  gr_fit_all <- gr_fit_all[!(is.na(gr_fit_all$mumax)),]
  gr_fit_all$mumax <- as.numeric(gr_fit_all$mumax)
  gr_fit_all$strain <- factor(gr_fit_all$strain, levels = levels(surv_all$strain))
  gr_fit_all$group <- interaction(gr_fit_all$strain, gr_fit_all$temperature,  sep = ":", lex.order = F)
  
  gr_fit_all$strain <- reorder(gr_fit_all$strain, gr_fit_all$mumax, median, decreasing = T)
  
  nrow(gr_fit_all)
  str(gr_fit_all)
  
  table(gr_fit_all$strain)
  
  # Descriptive statistics for Growth rate (mumax)
  
  desc_stats_gr <- NULL
  df4 <- stby(data = gr_fit_all$mumax, INDICES = gr_fit_all$group, FUN = descr)
  
  for (g in unique(gr_fit_all$group)) {
    df7 <- c(g, unlist(df4[g])[c(1,2,10,15)])
    desc_stats_gr <- rbind(desc_stats_gr, df7)
  }
  
  colnames(desc_stats_gr) <- c("group", "mean.SGR", "sd", "var", "replicates")
  desc_stats_gr <- as.data.frame(desc_stats_gr)
  
  head(desc_stats_gr)
  
  write.table(desc_stats_gr, file = paste0(output2, "SGR_allstrains_descStats.csv"), sep = ",", row.names = F)

  for (temp in c(16,26)) {
    
    gr_fit_all_temp <- subset(gr_fit_all, temperature == temp)

    par(mar = c(6, 10, 3, 8), xpd = T)
bp <- boxplot(mumax ~ strain, data = gr_fit_all_temp,
              #boxwex = 0.4, 
              horizontal = T,
              ylim = c(-1, 1),
              xlab = "", ylab = "",  outline = F, 
              xaxt = "n", yaxt = "n", 
              #col = cols,
              main = paste0("Growth Rate temperature: ", temp, "ºC")
              )

#tab_group <- length(unique(comparison$antibiotic))
xpos <- match(gr_fit_all_temp$strain, bp$names)

colors <- rainbow(length(unique(gr_fit_all_temp$day)), alpha = 0.8)
colors_match <- colors[gr_fit_all_temp$day]

points(jitter(gr_fit_all_temp$mumax, amount = max(gr_fit_all_temp$mumax)/30), xpos, col = rgb(1, 0, 0, alpha = 0.3), pch = 16, cex = 0.5) #may be wrong
#points(jitter(surv_all_filt$survival, amount = 0.2), xpos, col = colors_match, pch = 16, cex = 0.5) #may be wrong

mtext("Specific Growth Rate (day-1)", side = 1, line = 4, cex = 2)

# To add the label of x axis 
yscale <-  seq(par("xaxp")[1], par("xaxp")[2], length.out = par("xaxp")[3]*1.1)
xaxis <- par("xaxp")
#n_treatments <- levels(comparison$treatment)
n_comp <- unique(bp$names)
#n_vols <- unique(comparison$volume)

#Add xlabels
axis(2, at = 1:length(bp$names), labels = FALSE, tck = -0.01)
text(y = 1:length(bp$names), x = par("usr")[1]*1.05, labels = n_comp, cex = 1, adj = 1)
axis(1, at = yscale, labels = yscale, las = 1, cex.axis = 2, padj = 1)

#legend(x = par("usr")[2]*0.9, y = par("usr")[4]*1.1, legend = levels(surv_all_filt$day), col = colors, pch = 16, bty = "n", title = "Maximum Day")

######### Get pvalues and get the Compact Letter Display
#“All pairwise comparisons among the 53 strains (n = 1378) were conducted with multiple-testing correction,
#and results were summarized using a compact letter display following Piepho (2004).”


#Calculates all pairwise calculations
choose(length(unique(gr_fit_all_temp$strain)), 2)

#Perform Dunn Test
library("cld")
library("DescTools")

dt <- DunnTest(gr_fit_all_temp$mumax, gr_fit_all_temp$strain)

cld_res <- make_cld(dt)

letters <- as.data.frame(cld_res)

letters$group <- factor(letters$group, levels(gr_fit_all_temp$strain))
letters <- letters[order(letters$group),]

write.table(letters, file = paste0(output2, "SGR_pvals_SGR.csv"), sep = ",", row.names = F)

#add stats letters
text(bp$stats[5,] + 0.1, c(1:ncol(bp$stats)), labels = letters$cld, cex = 1, adj = 0)


dev.print(jpeg, file = paste0(output2, "SGR_allstrains_", temp, "C.jpg"), width = 2500, 
          height = 4000, res = 300, unit = "px") 
  }
  
  ################## Mixed model #################
  library(lme4)
  library(lmerTest)
  m <- lmer(mumax ~ strain * temperature + (1|well), data = comparison)
  class(m)
  
  library(emmeans)
  
  emm <- emmeans(m, ~ strain | temperature)
  pairs(emm, adjust = "tukey", infer = TRUE)

  emm_df <- as.data.frame(emm)
  
  emm_df$group <- cut(
    emm_df$emmean,
    breaks = quantile(emm_df$emmean, probs = c(0, .25, .75, 1)),
    labels = c("Low", "Medium", "High"),
    include.lowest = TRUE, na.rm = T
  )

  ggplot(emm_df, aes(x = emmean, y = reorder(strain, emmean))) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL)) +
    facet_wrap(~ temperature)
  #############################

desc_stats_surv <- NULL
df5 <- stby(data = surv_all$survival, INDICES = surv_all$strain, FUN = descr)

for (g in unique(surv_all_filt$strain)) {
  df6 <- c(g, unlist(df5[g])[c(1,2,10,15)])
  desc_stats_surv <- rbind(desc_stats_surv, df6)
}

colnames(desc_stats_surv) <- c("group","mean", "sd", "var", "replicates")
desc_stats_surv <- as.data.frame(desc_stats_surv)

write.table(desc_stats_surv, file = paste0(output2,"SGR_allstrain_stats.csv"), sep = ",", row.names = F)

end_time <- Sys.time()
total_time <- end_time - start_time
print(total_time)


#######################################
## Descriptive stats ##############3
#################################
library(summarytools);

desc_stats <- NULL
for (v in c("mumax")) {
  df1 <- stby(data = comparison$mumax, INDICES = comparison$group, FUN = descr)
  for (g in unique(comparison$group)) {
    df2 <- c(g,unlist(df1[g])[c(1,2,10)])
    desc_stats <- rbind(desc_stats,df2)
  }
}
colnames(desc_stats) <- c("group","mean", "sd", "var")
desc_stats <- as.data.frame(desc_stats)

#write.table(desc_stats, file = paste0("desc_stats_mumax_sexcols.csv"), sep = ",", row.names = F)
write.table(desc_stats, file = paste0("desc_stats_mumax_edeffect.csv"), sep = ",", row.names = F)


     
     
  
