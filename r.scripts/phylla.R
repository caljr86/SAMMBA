setwd("~/Dropbox/1_posdoc/1_CCMAR/experiments/384WP")
test <- "dilution"
plate <- "plate31"
date <- "250424"

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
  metadata <- design[design$WELL == well,][,c(1:2,5,6)]
  metadata2 <- metadata[rep(1, nrow(area)),]
  date1 <- rep(date, times = nrow(area))
  img <- cbind(photo, well2, area, metadata2, date1)
  colnames(img) <- c("photo", "well", "gametophyte", "area", "row", "column", "species", "treatment", "date" )
  seg_res_raw <- rbind(seg_res_raw, img)
}

str(seg_res_raw)
head(seg_res_raw)
seg_res_raw$gametophyte <- as.numeric(seg_res_raw$gametophyte)
seg_res_raw$column <- as.numeric(seg_res_raw$column)
seg_res_raw[[date]] <- list(seg_res_raw)

############ Singlets: Well with only 1 GM ##################

cpw_df <- NULL

for(w in unique(seg_res_raw$well)){
  cpw <- subset(seg_res_raw, well == w)
  counts <- nrow(cpw)
  area <- sum(cpw$area)
  species <- unique(cpw$species)
  img <- unique(cpw$photo)
  df <- cbind(w, img, counts, area, species)
  colnames(df) <- c("well", "photo", "counts", "area", "species")
  cpw_df <- as.data.frame(rbind(cpw_df, df))
}

cpw_df[cpw_df$counts == 1 & cpw_df$species == "Phylla",]

nrow(subset(seg_res_raw, species == "Phylla"))

singlets <- NULL

for(d in unique(seg_res_raw$species)){
  df <- subset(seg_res_raw, species == d)
  cpw <- table(df$well)
  ww1 <- as.numeric(sum(cpw == 1))
  df2 <- c(d, ww1)
  singlets <- as.data.frame(rbind(singlets, df2))
}

singlets$strain <- seq(1,nrow(singlets),1)
colnames(singlets) <- c("species", "counts", "strain")
singlets <- singlets[,c("strain", "counts", "species")]
singlets$counts <- as.numeric(singlets$counts)
singlets <- singlets[order(singlets$counts),]
singlets$species <- factor(c("Laminaria ochroleuca", "P. purpurascens") ,
                           levels = c("Laminaria ochroleuca", "P. purpurascens"))

singlets_phylla <- droplevels(subset(singlets, species == "P. purpurascens"))
levels(singlets_phylla$species)
singlets_phylla$strain <- "PP01OR"

