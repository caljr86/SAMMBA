
setwd("/home/cicero/Dropbox/1_posdoc/1_CCMAR/experiments/384isol/growthrate/")

library(growthrates)
plate <- "plate58 - modeltest"
setwd(paste0(getwd(), "/", plate, "/"))
output <- paste0(getwd(), "/")
design <- read.table(file = "plate_design.csv", header = T, sep = ",")
nrow(design)

gc384 <- NULL
for(f in list.dirs(recursive = F)){
  filter <- unlist(strsplit(f, "/"))[2]
  
  for(i in list.files(paste0(f, "/segmented_labkit/results/"))){
    # Skip to the next folder if no files are found
    if (length(f) == 0) {
      #cat("Skipping empty folder:", folder, "\n")
      next  # Go to the next iteration of the loop
    }
    area <- read.table(file = paste0(f, "/segmented_labkit/results/", i), header =T, sep = ",")
    area <- area[,1:2]
    colnames(area) <- c("gametophyte", "area")
    photo <- as.character(unlist(strsplit(i, "_segmented"))[1])
    well <- as.character(unlist(strsplit(i, "_"))[3])
    well2 <- rep(well, times = nrow(area))
    metadata <- design[design$WELL == well,][,c(1,2,5)]
    metadata2 <- metadata[rep(1, nrow(area)),]
    date1 <- rep(filter, times = nrow(area))
    img <- cbind(photo, well2, area, metadata2, date1)
    colnames(img) <- c("photo", "well", "gametophyte", "area", "row", "column", "plate", "filter" )
    gc384 <- rbind(gc384, img)
  }
}

str(gc384)
head(gc384)
gc384$gametophyte <- as.numeric(gc384$gametophyte)
gc384$column <- as.numeric(gc384$column)
summary(gc384)
nrow(gc384)

#Reorder wells sequence
gc384$well <- factor(gc384$well, levels = unique(paste0(gc384$row, gc384$column)))
gc384 <- gc384[order(gc384$row, gc384$column),]

summary(gc384)

quantile <- as.numeric(quantile(gc384$area, probs = 0.24, na.rm = FALSE,
                     names = TRUE, type = 7, digits = 2))

gc384 <- gc384[gc384$area > quantile,]
head(gc384)
nrow(gc384)
length(unique(gc384$well))

############ Correlation DAPI x PI
pi <- gc384[gc384$filter == "pi",]
dapi <- gc384[gc384$filter == "dapi",]
t = "area"
x <- pi[, t]         # shorter
y <- dapi[, t]       # longer

range_x <- mean(range(x))
  range_y <- mean(range(y))
  min_x <- min(x)
  max_x <- max(x)
  min_y <- min(y)
  max_y <- max(y)
  
  #Kolmogorov-Smirnoff test
  ks.test(x, y)
   
  #Mann-Whitney test
    wilcox.test(x,y)
  
# DTW
    
    # Install if needed
    #install.packages("dtw")
    
    # Load library
    library(dtw)
    
    y_interp <- approx(1:length(y), y, n = length(x))$y
        # Perform DTW
    alignment <- dtw(y_interp, x, keep = TRUE)
    
 matched_df <- data.frame(
   df1_index = alignment$index1,
   df1_value = y_interp[alignment$index1],
   df2_index = alignment$index2,
   df2_value = x[alignment$index2]
 )
 
 plot(matched_df$df1_value ~ matched_df$df2_value) 
 
 model <- lm(matched_df$df1_value ~ matched_df$df2_value)
 
  sum_model <- summary(model)
  col <- c("blue", "red")
  plot(pi[,t] ~ dapi[,t],  type = "n", xlim = c(min_x*0.8,max_x*1.1), ylim = c(min_y*0.8, max_y*1.1),
       main = paste0("Total ", t, " per well"), xlab ="", ylab = "")
  
  points(jitter(dapi[,t], amount = min_x*0.1), jitter(pi[,t], amount = min_y*0.1), pch = 19)
  
  if (t == "area"){ 
    mtext(expression("DAPI Area("*mu*"m²)"), side = 1, line = 2, cex = 0.9, font = 2)
    mtext(expression("PI Area("*mu*"m²)"), side = 2, line = 2, cex = 0.9, font = 2)}
  
  if (t == "counts"){ 
    mtext("DAPI Counts", side = 1, line = 2, cex = 0.9, font = 2)
    mtext("PI Counts", side = 2, line = 2, cex = 0.9, font = 2)}
  
  
  # Get the plot limits
  x_limits <- range(dapi[,t])  # x-axis limits based on the data
  y_limits <- range(pi[,t])  # y-axis limits based on the data
  
  # Calculate y values for the x-axis limits
  x_range <- par("usr")[1:2]  # Current x-axis limits
  y_range <- coef(model)[1] + coef(model)[2] * x_range  # y = intercept + slope * x
  
   # Ensure the line stays within the plotting area
  lines(x_range, y_range, col = "red", lwd = 2)
  
  text(mean(c(range_x, min_x)), mean(c(max_y, range_y)) , paste0("R² =", round(sum_model$r.squared,4)))
  pvalue <- signif(sum_model$coefficients[2,4],4)
  text(mean(c(range_x, min_x)), mean(c(max_y, range_y))*0.95, paste0("pvalue =", pvalue))
  

dev.print(jpeg, file = paste0(output, plate, "_filter_optmz_per_gam_", t, ".jpg"), width = 5000, 
          height = 2000, res = 300, unit = "px") 

####### Total Area per well ###########

gc384_apw <- aggregate(area ~ well + filter, data = gc384, sum)
gc384_apw$counts <- aggregate(gametophyte ~ well+filter, data = gc384, length)$gametophyte
gc384_apw$row <- substr(gc384_apw$well,1,1)
gc384_apw$col <- as.numeric(gsub("[A-Z]", "", gc384_apw$well))
gc384_apw <- gc384_apw[order(gc384_apw$row, gc384_apw$col),]
gc384_apw$filter <- factor(gc384_apw$filter, levels = c("dapi", "pi"))


head(gc384_apw)
nrow(gc384_apw)
str(gc384_apw)
summary(gc384_apw)

############ Correlation DAPI x PI
pi <- gc384_apw[gc384_apw$filter == "pi",]
dapi <- gc384_apw[gc384_apw$filter == "dapi",]

par(mfrow = c(1,2),
    mar = c(4,5,4,6), xpd = T)
for(t in c("counts", "area")){
  range_x <- mean(range(dapi[,t]))
  range_y <- mean(range(pi[,t]))
  min_x <- min(dapi[,t])
  max_x <- max(dapi[,t])
  min_y <- min(pi[,t])
  max_y <- max(pi[,t])
  model <- lm(pi[,t] ~ dapi[,t])
  
sum_model <- summary(model)

col <- c("blue", "red")
plot(pi[,t] ~ dapi[,t],  type = "n", xlim = c(min_x*0.8,max_x*1.1), ylim = c(min_y*0.8, max_y*1.1),
     main = paste0("Total ", t, " per well"), xlab ="", ylab = "")

points(jitter(dapi[,t], amount = min_x*0.1), jitter(pi[,t], amount = min_y*0.1), pch = 19)

if (t == "area"){ 
  mtext(expression("DAPI Area("*mu*"m²)"), side = 1, line = 2, cex = 0.9, font = 2)
  mtext(expression("PI Area("*mu*"m²)"), side = 2, line = 2, cex = 0.9, font = 2)}

if (t == "counts"){ 
  mtext("DAPI Counts", side = 1, line = 2, cex = 0.9, font = 2)
  mtext("PI Counts", side = 2, line = 2, cex = 0.9, font = 2)}
  
  
# Get the plot limits
x_limits <- range(dapi[,t])  # x-axis limits based on the data
y_limits <- range(pi[,t])  # y-axis limits based on the data

# Calculate y values for the x-axis limits
x_range <- par("usr")[1:2]  # Current x-axis limits
y_range <- coef(model)[1] + coef(model)[2] * x_range  # y = intercept + slope * x

# Ensure the line stays within the plotting area
lines(x_range, y_range, col = "red", lwd = 2)

text(mean(c(range_x, min_x)), mean(c(max_y, range_y)) , paste0("R² =", round(sum_model$r.squared,4)))
pvalue <- signif(sum_model$coefficients[2,4],4)
text(mean(c(range_x, min_x)), mean(c(max_y, range_y))*0.95, paste0("pvalue =", pvalue))

}

dev.print(jpeg, file = paste0(output, plate, "_filter_optmz_per_well_", t, ".jpg"), width = 5000, 
          height = 2000, res = 300, unit = "px") 

#####################################################
########## Correlation LabKit x Manual ###############
#####################################################

setwd("/media/cicero/T7/Imagens/experiments_kelps/ccmar/384WP/growthcurve")
plate <- "plate58 - modeltest"
setwd(paste0(getwd(), "/", plate, "/"))
output <- paste0(getwd(), "/")
model_optm <- read.table(file = "pi/manual/model_optmization.csv", header = T, row.names = NULL, sep = ",")
head(model_optm)

wilcox <- wilcox.test(model_optm$Area ~ model_optm$analysis)
wilcox


manual <- subset(model_optm, analysis == "manual")
labkit <- subset(model_optm, analysis == "labkit")

nrow(manual)
nrow(labkit)

model <- lm(labkit$Area ~ manual$Area)

sum_model <- summary(model)

# Get the plot limits
x_limits <- range(manual$Area)  # x-axis limits based on the data
x_min <- min(manual$Area)
x_max <- max(manual$Area)
y_limits <- range(labkit$Area)  # y-axis limits based on the data
y_min <- min(labkit$Area)
y_max <- max(labkit$Area)

par(mar = c(5.1, 7, 4.1, 2.1))

#Plot regression with labkit on y-axis as manual is the reference
plot(labkit$Area ~ manual$Area,  xlim = c(x_min*0.8, x_max*1.1), ylim = c(y_min*0.8, y_max*1.1), pch = 19,
     main = paste0("Linear regression"), xlab = "" , ylab = "")

# Add y-axis title
mtext(expression("Labkit Area("*mu*"m²)"), side = 2, line = 3, cex = 3)
mtext(expression("Manual Area("*mu*"m²)"), side = 1, line = 4, cex = 3)

# Calculate y values for the x-axis limits
x_range <- c(0, x_max)  # Current x-axis limits
y_range <- coef(model)[1] + coef(model)[2] * x_range  # y = intercept + slope * x

# Ensure the line stays within the plotting area
lines(x_range, y_range, col = "red", lwd = 2)

text(mean(c(x_limits, x_min)), mean(c(y_max, y_range))*1.2 , paste0("R² =", round(sum_model$r.squared,4)), cex = 3)
pvalue <- signif(sum_model$coefficients[2,4],4)
text(mean(c(x_range, x_min))*1.5, mean(c(y_max, y_range))*1.5, paste0("pvalue =", pvalue), cex = 3)

dev.print(jpeg, file = paste0(output, "_model_optmz_per_gametophyte.jpg"), width = 2500, 
          height = 2000, res = 300, unit = "px") 

#### GGPLOT ###


df <- data.frame(manual = manual$Area,
                 labkit = labkit$Area)

# ---- 1. Linear regression plot with CI ----
p1 <- ggplot(df, aes(x = manual, y = labkit)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(title = "Linear regression",
       x = expression("Manual Area("*mu*"m²)"), y = expression("Labkit Area("*mu*"m²)")) +
  theme_minimal()

# ---- 2. Bland–Altman plot ----
# Differences vs mean of methods
df$mean <- rowMeans(df)
df$diff <- df$manual - df$labkit

mean_diff <- mean(df$diff)                  # bias
sd_diff <- sd(df$diff)                      # standard deviation of differences
loa_upper <- mean_diff + 1.96 * sd_diff     # upper limit of agreement
loa_lower <- mean_diff - 1.96 * sd_diff     # lower limit of agreement

# Build annotation text
annot_text <- paste0("Bias = ", round(mean_diff, 2),
                     "\nLoA = ", round(loa_lower, 2),
                     " to ", round(loa_upper, 2))

ggplot(df, aes(x = mean, y = diff)) +
  geom_point() +
  geom_hline(yintercept = mean_diff, color = "red", linetype = "dashed") +
  geom_hline(yintercept = loa_upper, color = "darkgreen", linetype = "dotted") +
  geom_hline(yintercept = loa_lower, color = "darkgreen", linetype = "dotted") +
  annotate("text", x = max(df$mean), y = loa_upper,
           label = annot_text, hjust = 1, vjust = 1.2, size = 3.5) +
  labs(title = "Bland–Altman Plot",
       x = "Mean of methods (Manual & LabKit)",
       y = "Difference (Manual - LabKit)") +
  theme_minimal()


dev.print(jpeg, file = paste0(output, "_model_optmz_bland.altman_plot.jpg"), width = 2500, 
          height = 2000, res = 300, unit = "px") 
