## Dilution-to-extiction analysis in a larger scale: 196 wells per sex

setwd("~/Dropbox/1_posdoc/1_CCMAR/experiments/384isol")
test <- "dilution"
plate <- "plate32"
date <- "250506"

output <- paste0(getwd(), "/", test, "/", plate, "/", date, "/")

setwd(paste0(getwd(), "/", test, "/", plate, "/", date, "/segmented_labkit/results/"))

design <- read.table(file = "../../../plate_design.csv", header = T, sep = ",")


seg_res_raw <- NULL
for(i in list.files()){
  area <- read.table(file = i, header =T, sep = ",")
  area <- area[,1:2]
  colnames(area) <- c("gametophyte", "area")
  photo <- as.character(unlist(strsplit(i, "_segmented"))[1])
  well <- as.character(unlist(strsplit(i, "_"))[2])
  #well <- design[design$photo == photo, ][,3]
  well2 <- rep(well, times = nrow(area))
  metadata <- design[design$WELL == well,][,c(1:2,5:7)]
  metadata2 <- metadata[rep(1, nrow(area)),]
  date1 <- rep(date, times = nrow(area))
  img <- cbind(photo, well2, area, metadata2, date1)
  colnames(img) <- c("photo", "well", "gametophyte", "area", "row", "column", "species", "treatment", "sex", "date" )
  seg_res_raw <- rbind(seg_res_raw, img)
}

str(seg_res_raw)
head(seg_res_raw)
seg_res_raw$gametophyte <- as.numeric(seg_res_raw$gametophyte)
seg_res_raw$column <- as.numeric(seg_res_raw$column)

#Check ranges. Well A2 was contaminated
seg_res_raw2 <- subset(seg_res_raw, well != "A2" & well != "M4")
                       
seg_res_raw2[seg_res_raw2$area %in% range(seg_res_raw2$area),]



######### GM dilution #########

cpw_df <- NULL

for(w in unique(seg_res_raw$well)){
  cpw <- subset(seg_res_raw, well == w)
  counts <- nrow(cpw)
  area <- sum(cpw$area)
  sex <- unique(cpw$sex)
  df <- cbind(w, counts, area, sex)
  colnames(df) <- c("well", "counts", "area", "sex")
  cpw_df <- as.data.frame(rbind(cpw_df, df))
}

nrow(cpw_df)
head(cpw_df)
cpw_df$counts

cpw_df$area <- as.numeric(cpw_df$area)
cpw_df$counts <- as.numeric(cpw_df$counts)
cpw_df$apc <- cpw_df$area/cpw_df$counts 

#Check counts per wells
cpw_df[cpw_df$counts >= 4,]

#Some wells had contamination, so I had to filter 3 wells
cpw_df_norm <- cpw_df[cpw_df$counts <= 3,]
nrow(cpw_df_norm)

#Number of wells with males
nrow(subset(cpw_df_norm, sex == "male"))
#Number of male singlets
nrow(subset(cpw_df_norm, counts == 1 & sex == "male"))
#Wells with male singlets
subset(cpw_df_norm, counts == 1 & sex == "male")$well
#Total number of detected male fragments
sum(subset(cpw_df_norm, sex == "male")$counts)

#Number of wells with females
nrow(subset(cpw_df_norm, sex == "female"))
#Number of female singlets
nrow(subset(cpw_df_norm, counts == 1 & sex == "female"))
#Wells with female singlets
subset(cpw_df_norm, counts == 1 & sex == "female")$well
#Total number of detected female fragments
sum(subset(cpw_df_norm, sex == "female")$counts)

cpw_df_norm$counts

#######Other than singlets
# DOublets
nrow(subset(cpw_df_norm, counts == 2 & sex == "male"))
subset(cpw_df_norm, counts == 2 & sex == "male")$well
nrow(subset(cpw_df_norm, counts == 2 & sex == "female"))
#Triplets
nrow(subset(cpw_df_norm, counts == 3 & sex == "male"))
nrow(subset(cpw_df_norm, counts == 3 & sex == "female"))

###################################################

singlets <- NULL

for(s in unique(seg_res_raw$sex)){
  df1 <- subset(seg_res_raw, sex == s)
    cpw <- table(df1$well)
    ww1 <- sum(cpw == 1)
    df <- c(ww1, "singlet", s)
    singlets <- as.data.frame(rbind(singlets, df))
  }

head(singlets)

doublets <- NULL

for(s in unique(seg_res_raw$sex)){
  df1 <- subset(seg_res_raw, sex == s)
  cpw <- table(df1$well)
  ww1 <- sum(cpw == 2)
  df <- c(ww1, "doublet", s)
  doublets <- as.data.frame(rbind(doublets, df))
}

head(doublets)

triplets <- NULL

for(s in unique(seg_res_raw$sex)){
  df1 <- subset(seg_res_raw, sex == s)
  cpw <- table(df1$well)
  ww1 <- sum(cpw == 3)
  df <- c(ww1, "triplet", s)
  triplets <- as.data.frame(rbind(triplets, df))
}

head(triplets)

multiplets <- NULL

for(s in unique(seg_res_raw$sex)){
  df1 <- subset(seg_res_raw, sex == s)
  cpw <- table(df1$well)
  ww1 <- sum(cpw >= 4)
  df <- c(ww1, "multiplet", s)
  multiplets <- as.data.frame(rbind(multiplets, df))
}

head(multiplets)

isolation <- as.data.frame(rbind(singlets, doublets, triplets, multiplets))
colnames(isolation) <- c("counts", "level", "sex")
isolation$counts <- as.numeric(isolation$counts)
str(isolation)
head(isolation)

zeros_males <- 196 - sum(subset(isolation, sex == "male")$counts)
zeros_males <- c(zeros_males, "empty", "male")
zeros_females <- 196 - sum(subset(isolation, sex == "female")$counts)
zeros_females <- c(zeros_females, "empty", "female")

isolation <- rbind(isolation, zeros_males, zeros_females)
isolation$level <- factor(isolation$level, levels = c("singlet", "doublet", "triplet", "multiplet", "empty"))
isolation$counts <- as.numeric(isolation$counts)
isolation
str(isolation)

# PLot isolation
dev.off()
col <- gray(c(2,4,6,8,10)/10)

par(mar = c(5, 5, 6, 8),
    oma = c(0,0,0,0), xpd = TRUE)

  bp <- barplot(counts ~ level + sex, data = isolation,  col = col, 
                ylab = "Number of Wells", xaxt = "n", 
                yaxt = "n", xlab = "Sex", ylim = c(min(isolation$counts), max(isolation$counts)*1.1))
  #x axis
  axis(1, at = bp, labels = F)
  text(x = bp, y = -5, labels = unique(isolation$sex), xpd = T)
  
  #y axis
  axis(2, at = seq(0,196,12), labels = F)
  text(x = -0.05, y = seq(0,196,28), #srt = 360, #xpd = T,
       labels = seq(0,196,28))
  
  #Title
  title(main = "Large-Scale DTE isolation")
  

# Add legend in the upper-right outer margin
par(fig = c(0, 1, 0, 1), new = TRUE)
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")  # empty plot
legend("topright", legend = levels(isolation$level),
       pt.bg = col, col = "black", pch = 21, bty = "n", title = "Well isolation level", cex = 1,
       inset = c(-0.35, -0.1), y.intersp = 0.5)

dev.print(jpeg, file = paste0(output, "Large_scale_DTE_isolation.jpg"), width = 2000, 
          height = 2500, res = 300, unit = "px") 

write.table(isolation, file = paste0(output, "large_scale_isolation.csv"), sep = ",", col.names = T, row.names = F)
