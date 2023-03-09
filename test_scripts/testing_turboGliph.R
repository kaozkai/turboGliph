require(turboGliph)
data("hs_CD8_ref")

data_sample <- hs_CD8_ref[sample(x = 1:nrow(hs_CD8_ref), size = 2000, replace = F), 1:3]
data_sample$sample_id <- "s"

data_ref <- hs_CD8_ref[, 1:3]
ks <- c(2, 3, 4)
B <- 1000
cores <- 1

source("R/util_v1.R")
source("R/util_v2.R")
source("R/util_v1_v2.R")
source("R/gliph_v1.R")
source("R/gliph_v2.R")


#### Gliph 1-2 + Jan's versions #####


t_gliph_v1_a <- Sys.time()
o1 <- gliph_v1(
    data_sample = data_sample,
    data_ref = data_ref,
    ks = ks,
    B = B,
    cores = cores,
    control = NULL)
t_gliph_v1_b <- Sys.time()


t_gliph_v2_a <- Sys.time()
o2 <- gliph_v2(
    data_sample = data_sample,
    data_ref = data_ref,
    ks = ks,
    cores = cores,
    control = NULL)
t_gliph_v2_b <- Sys.time()


t_jan_gliph_v1_a <- Sys.time()
o3 <- turboGliph::turbo_gliph(
    cdr3_sequences = data_sample,
    result_folder = "/home/sktron/Desktop/tmp/",
    refdb_beta = data_ref,
    lcminp = 0.05,
    gccutoff = 1,
    structboundaries = F,
    boundary_size = 0,
    cluster_min_size = 1,
    accept_sequences_with_C_F_start_end = FALSE)
t_jan_gliph_v1_b <- Sys.time()


t_jan_gliph_v2_a <- Sys.time()
o4 <- turboGliph::gliph2(
    cdr3_sequences = data_sample,
    result_folder = "/home/sktron/Desktop/tmp/",
    refdb_beta = data_ref,
    lcminp = 0.05,
    # gccutoff = 1,
    # structboundaries = F,
    # boundary_size = 0,
    # cluster_min_size = 1,
    accept_sequences_with_C_F_start_end = FALSE)
t_jan_gliph_v2_b <- Sys.time()


t_gliph_v

o1 <- o1[[1]]
o2 <- o2[[1]]

e1 <- o1$motif_enrichment
e2 <- o2$motif_enrichment
e3 <- o3$motif_enrichment$all_motifs
e4 <- o4$motif_enrichment$all_motifs


e1 <- e1[e1$filter == TRUE,]
e2 <- e2[e2$filter == TRUE,]
e3 <- o3$motif_enrichment$selected_motifs
e4 <- o4$motif_enrichment$selected_motifs


intersect(e1$motif, e3$Motif)
setdiff(e1$motif, e3$Motif)
setdiff(e3$Motif, e1$motif)


t <- o1$global_pairs
length(unique(c(unique(cdr3)[t[, 1]],
                unique(cdr3)[t[, 2]])))
unique_t <- unique(c(unique(cdr3)[t[, 1]],
                     unique(cdr3)[t[, 2]]))


nrow(o1$global_pairs)
w <- o2$connections
w <- w[w[, 3]=="global", ]
length(unique(c(w[, 1], w[, 2])))
unique_w <- unique(c(w[, 1], w[, 2]))


setdiff(unique_t, unique_w)
setdiff(unique_w, unique_t)




w <- cbind(w, stringdist::stringdist(a = w[,1], b = w[,2], 
                                     method = "hamming"))
# w <- w[w[, 3] == "global", ]
w <- w[, 4] == 1
table(w[,4])


d <- stringdist::stringdistmatrix(a = cdr3, b = cdr3, method = "hamming")
q <- which(d<=1, arr.ind = T)








#### Gliph 2 #####
require(turboGliph)
data("hs_CD8_ref")

cdr3 <- sample(x = hs_CD8_ref$data$CDR3b, size = 5000, replace = F)
cdr3_ref <- hs_CD8_ref$data$CDR3b
ks <- c(2, 3, 4)
B <- 1000
cores <- 1
control <- NULL
source("R/gliph_v2.R")
control <- get_control(control_in = control)
















