setwd("~/Dropbox/1_posdoc/1_CCMAR/experiments/384WP")
test <- "dilution"
plate <- "plate44"
date <- "250424"

output <- paste0(getwd(), "/", test, "/", plate, "/", date, "/")

setwd(paste0(getwd(), "/", test, "/", plate, "/", date, "/segmented_labkit/results/"))

design <- read.table(file = "../../../plate_design.csv", header = T, sep = ",")

head(design)

seg_res_raw <- NULL
for (i in list.files()) {
  area <- read.table(file = i, header = T, sep = ",")
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
  colnames(img) <- c("photo", "well", "gametophyte", "area", "row", "column", "species", "strain",  "treatment", "date" )
  seg_res_raw <- rbind(seg_res_raw, img)
}

str(seg_res_raw)
head(seg_res_raw)
seg_res_raw$gametophyte <- as.numeric(seg_res_raw$gametophyte)
seg_res_raw$column <- as.numeric(seg_res_raw$column)
#seg_res_raw[[date]] <- list(seg_res_raw)

unique(seg_res_raw$species)

#### Counts per well


cpw_df <- NULL

for (w in unique(seg_res_raw$well)) {
  cpw <- subset(seg_res_raw, well == w)
  counts <- nrow(cpw)
  area <- sum(cpw$area)
  species <- unique(cpw$species)
  df <- cbind(w, counts, area, species)
  colnames(df) <- c("well", "counts", "area", "species")
  cpw_df <- as.data.frame(rbind(cpw_df, df))
}

nrow(cpw_df)
head(cpw_df)
unique(cpw_df$species)


cpw_df$area <- as.numeric(cpw_df$area)
cpw_df$counts <- as.numeric(cpw_df$counts)
cpw_df$apc <- cpw_df$area/cpw_df$counts 

table(cpw_df[cpw_df$counts == 1,]$species)
table(cpw_df[cpw_df$counts <= 2,]$species)
table(cpw_df[cpw_df$counts <= 3,]$species)
table(cpw_df[cpw_df$counts <= 4,]$species)
table(cpw_df[cpw_df$counts <= 5,]$species)

#Number of wells with males
hydrolithon <- subset(cpw_df, species == "Hydrolithon")
nrow(hydrolithon)
hydrolithon_singlets <- subset(hydrolithon, counts == 1)
nrow(hydrolithon_singlets)

############ Singlets: Well with only 1 GM ##################

singlets <- NULL

for (d in unique(seg_res_raw$strain)) {
  df <- subset(seg_res_raw, strain == d)
  cpw <- table(df$well)
  species <- unique(df$species)
  ww1 <- as.numeric(sum(cpw == 1))
  df2 <- c(d, ww1, species)
  singlets <- as.data.frame(rbind(singlets, df2))
}

head(singlets)
colnames(singlets) <- c("strain", "counts", "species")
singlets$species[singlets$species == "hydrolithon"] <- "Hydrolithon sp"

singlets$counts <- as.numeric(singlets$counts)
singlets <- singlets[order(singlets$counts),]
singlets$species <- factor(c("Hydrolithon sp", "Erythrotrichia sp", "Hydrolithon sp",  "Ulva sp", 
                             "Hydrolithon sp", "Halymenia sp","Halymenia sp","Halymenia sp"),
                           levels = c("Erythrotrichia sp", "Halymenia sp", "Hydrolithon sp", "Ulva sp"))
singlets_rng <- singlets

#Merge with singlets from Phyllariopsis
singlets_all <- rbind(singlets_rng, singlets_phylla)
levels(singlets_all$species)
sum(singlets_all$counts)

#counts_matrix <- xtabs(counts ~ strain + species, data = singlets_all)

reds1 <- colorRampPalette(c("#EE7777", "#cc3333", "#AA2222"))(3)  # lighter red to deep red
reds2 <- colorRampPalette(c("#EE7777", "#cc3333", "#AA2222"))(4)
other <- c("#339976", "#b68600")

col <- c(reds1, other[1], reds2, other[2])
names(col) <- singlets_all$strain
#"#EE7777" "#E26060" "#D74949" "#CC3333" "#C02D2D" "#B52727" "#AA2222" "#339976" "#b68600"
par(mar = c(4,8,2,5),
    xpd = T)
bp <- barplot(singlets_all$counts, names.arg = singlets_all$species, col = col, main = "Reds and Greens", horiz = T, yaxt = "n",
             xlab = "Singlet count", ylab = "", xlim = c(0, max(singlets_all$counts)*1.1))

#bp <- barplot(counts_matrix, col = col[rownames(counts_matrix)], names.arg = colnames(counts_matrix), 
 #       main = "Strain Isolation", horiz = T, yaxt = "n",  xlab = "Singlet count", ylab = "", xlim = c(0, 70))

# Add a legend outside the plot
legend("topright", inset = c(-0.2, 0), legend = singlets_all$strain, fill = col, bty = "n")

axis(2, at = bp, labels = F)
text(y = bp, x = -2, labels = singlets_all$species, adj = 1, offset = 2, xpd = T, font = 3)

dev.print(jpeg, file = paste0(output, date, "_singlets_all.jpg"), width = 2000, 
          height = 1500, res = 300, unit = "px")  

####

subset(seg_res_raw)

