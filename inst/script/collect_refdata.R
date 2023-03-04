# create temporary directory and download the scRNA-seq data:
if(dir.exists("temp_folder")==FALSE) {
    dir.create(path = "temp_folder")
}

# download data (last tested on 02. March 2023) 
utils::download.file(
    url = "http://50.255.35.37:8080/downloads/human_v2.0.zip",
    destfile = "temp_folder/human_v2.0.zip")

utils::download.file(
    url = "http://50.255.35.37:8080/downloads/mouse_v1.0.zip",
    destfile = "temp_folder/mouse_v1.0.zip")

# unzip
utils::unzip(zipfile = "temp_folder/human_v2.0.zip", exdir = "temp_folder/")
utils::unzip(zipfile = "temp_folder/mouse_v1.0.zip", exdir = "temp_folder/")


# human (hs) CD8
CD8 <- read.csv(file = "temp_folder/human_v2.0/ref_CD8_v2.0.txt", sep = "\t", header = FALSE)
colnames(CD8) <- c("CDR3b", "TRBV", "TRBJ")
CD8_CDR3b_length <- read.csv(file = "temp_folder/human_v2.0/ref_L_CD8_v2.0.txt", sep = "\t", header = FALSE)
colnames(CD8_CDR3b_length) <- c("CDR3_length_AA", "p")
CD8_TRBV_usage <- read.csv(file = "temp_folder/human_v2.0/ref_V_CD8_v2.0.txt", sep = "\t", header = FALSE)
colnames(CD8_TRBV_usage) <- c("TRBV", "p")
hs_CD8_ref <- list(data = CD8, 
                   CDR3b_length = CD8_CDR3b_length, 
                   TRBV_usage = CD8_TRBV_usage)
rm(CD8, CD8_CDR3b_length, CD8_TRBV_usage)
save(hs_CD8_ref, file = "hs_CD8_ref.RData")


# human (hs) CD4
CD4 <- read.csv(file = "temp_folder/human_v2.0/ref_CD4_v2.0.txt", sep = "\t", header = FALSE)
colnames(CD4) <- c("CDR3b", "TRBV", "TRBJ")
CD4_CDR3b_length <- read.csv(file = "temp_folder/human_v2.0/ref_L_CD4_v2.0.txt", sep = "\t", header = FALSE)
colnames(CD4_CDR3b_length) <- c("CDR3_length_AA", "p")
CD4_TRBV_usage <- read.csv(file = "temp_folder/human_v2.0/ref_V_CD4_v2.0.txt", sep = "\t", header = FALSE)
colnames(CD4_TRBV_usage) <- c("TRBV", "p")
hs_CD4_ref <- list(data = CD4, 
                   CDR3b_length = CD4_CDR3b_length, 
                   TRBV_usage = CD4_TRBV_usage)
rm(CD4, CD4_CDR3b_length, CD4_TRBV_usage)
save(hs_CD4_ref, file = "hs_CD4_ref.RData")



# mouse (mm) CD8
CD8 <- read.csv(file = "temp_folder/mouse_v1.0/ref_CD8_ms.txt", sep = "\t", header = FALSE)
colnames(CD8) <- c("CDR3b", "TRBV", "TRBJ")
CD8_CDR3b_length <- read.csv(file = "temp_folder/mouse_v1.0/ref_L_CD8_ms.txt", sep = "\t", header = FALSE)
colnames(CD8_CDR3b_length) <- c("CDR3_length_AA", "p")
CD8_TRBV_usage <- read.csv(file = "temp_folder/mouse_v1.0/ref_V_CD8_ms.txt", sep = "\t", header = FALSE)
colnames(CD8_TRBV_usage) <- c("TRBV", "p")
mm_CD8_ref <- list(data = CD8, 
                   CDR3b_length = CD8_CDR3b_length, 
                   TRBV_usage = CD8_TRBV_usage)
rm(CD8, CD8_CDR3b_length, CD8_TRBV_usage)
save(mm_CD8_ref, file = "mm_CD8_ref.RData")


# mouse (mm) CD4
CD4 <- read.csv(file = "temp_folder/mouse_v1.0/ref_CD4_ms.txt", sep = "\t", header = FALSE)
colnames(CD4) <- c("CDR3b", "TRBV", "TRBJ")
CD4_CDR3b_length <- read.csv(file = "temp_folder/mouse_v1.0/ref_L_CD4_ms.txt", sep = "\t", header = FALSE)
colnames(CD4_CDR3b_length) <- c("CDR3_length_AA", "p")
CD4_TRBV_usage <- read.csv(file = "temp_folder/mouse_v1.0/ref_V_CD4_ms.txt", sep = "\t", header = FALSE)
colnames(CD4_TRBV_usage) <- c("TRBV", "p")
mm_CD4_ref <- list(data = CD4, 
                   CDR3b_length = CD4_CDR3b_length, 
                   TRBV_usage = CD4_TRBV_usage)
rm(CD4, CD4_CDR3b_length, CD4_TRBV_usage)
save(mm_CD4_ref, file = "mm_CD4_ref.RData")



# get cluster sizes n = 100, 200, 300, ... , 1000, 2000, .... 