require(turboGliph)
data("hs_CD8_ref")

cdr3 <- sample(x = hs_CD8_ref$data$CDR3b, size = 20000, replace = F)
cdr3_ref <- hs_CD8_ref$data$CDR3b
ks <- c(2, 3, 4)
B <- 1000
cores <- 1
# control <- NULL
# control <- get_control(control_in = control)

source("R/gliph_v1.R")


t1 <- Sys.time()
o1 <- gliph_v1(
    cdr3 = unique(cdr3),
    cdr3_ref = cdr3_ref,
    ks = ks,
    B = B,
    cores = cores,
    control = NULL)
t2 <- Sys.time()

t3 <- Sys.time()
o2 <- turboGliph::turbo_gliph(
    cdr3_sequences = cdr3,
    result_folder = "/home/sktron/Desktop/tmp/",
    refdb_beta = hs_CD8_ref$data,
    lcminp = 0.05,
    gccutoff = 1,
    structboundaries = F,
    boundary_size = 0,
    cluster_min_size = 1,
    accept_sequences_with_C_F_start_end = FALSE)
t4 <- Sys.time()

t4-t3
t2-t1


e1 <- o1$motif_enrichment[o1$motif_enrichment$filter == TRUE,]
e2 <- o2$motif_enrichment$selected_motifs
intersect(e1$motif, e2$Motif)


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
