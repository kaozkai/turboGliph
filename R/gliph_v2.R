
# Description:
# This is rewritten code for GLIPH v2 (details are missing, TODO)
gliph_v2 <- function(cdr3,
                     cdr3_ref,
                     ks = c(2, 3, 4),
                     B = 1000,
                     cores = 1,
                     control = NULL) {
    
    # 0. parameter check
    
    
    # a. control check
    control <- get_control(control_in = control)
    
    
    
    # b. filter sample data based on inputs (e.g. remove CDR3 with small size, 
    # remove CDR3s without C and F, ...,)
    
    
    
    # 1. load reference
    
    
    # 2. local clustering
    # a. get local motifs
    motifs <- lapply(X = ks, 
           FUN = get_motifs_v2,
           cdr3 = cdr3, 
           cdr3_ref = cdr3_ref,
           min_o = control$filtering$local_min_o) 
    motifs <- do.call(rbind, motifs)
    
    
    # b. compute enrichment with fisher's exact test
    motif_stats <- t(apply(
        X = motifs[, c("f_sample", "f_ref", "n_sample", "n_ref")],
        MARGIN = 1, FUN = get_motif_enrichment_fet_v2))
    motifs$ove <- motif_stats[,1]
    motifs$p_value <- motif_stats[,2]
    rm(motif_stats)
    
    
    # c. add filter flag
    motif_enrichment <- get_motif_filter_v2(
        m = motifs,
        ks = ks,
        min_p = control$filtering$local_min_p, 
        min_ove = control$filtering$local_min_ove, 
        min_o = control$filtering$local_min_o)
    
    
    # d. find motifs in CDR3
    local_pairs <- get_local_pair(
        cdr3 = cdr3, 
        motif = motif_enrichment$motif[motif_enrichment$filter==TRUE])
    
    
    # 3. global TODO SK
    global_pairs <- get_global_pairs(
        cdr3 = cdr3, 
        global_max_dist = control$filtering$global_max_dist)
    
    
    # 4. scoring -> new function
    
    
    
    # 5. return TODO: format properly
    return(list(local_pairs = local_pairs,
                global_pairs = global_pairs,
                motif_enrichment = motif_enrichment,
                control = control))
}



# Description:
# Computes motif frequencies for a sample and reference
get_motifs_v2 <- function(x, cdr3, cdr3_ref, min_o) {
    
    # find kmers in sample
    kmers_s <- stringdist::qgrams(cdr3, q = x)
    if(ncol(kmers_s)==0) {
        # in real applications this should not happen -> input checks should
        # catch such errors
        stop("no kmers found in sample")
    }
    kmers_s <- kmers_s[1, kmers_s[1,]>=min_o]
    if(length(kmers_s)==0) {
        # in real applications this should not happen -> input checks should
        # catch such errors
        stop("no kmers found in sample")
    }
    
    
    # find kmers in reference
    kmers_r <- stringdist::qgrams(cdr3_ref, q = x)
    if(ncol(kmers_r)==0) {
        # in real applications this should not happen -> input checks should
        # catch such errors
        stop("no kmers found in reference")
    }
    kmers_r <- kmers_r[1,]
    
    # we are only interested in enrichment of motifs in sample relative to
    # reference. Remove all motifs from reference not found in sample.
    kmers_r <- kmers_r[names(kmers_r) %in% names(kmers_s)]
    
    
    # convert table to data.frame
    kmers_s <- data.frame(motif = names(kmers_s), 
                          f_sample = as.numeric(kmers_s))
    kmers_s$n_sample <- sum(kmers_s$f_sample)
    kmers_r <- data.frame(motif = names(kmers_r), 
                          f_ref = as.numeric(kmers_r))
    kmers_r$n_ref <- sum(kmers_r$f_ref)
    
    m <- merge(x = kmers_s, y = kmers_r, by = "motif", all = TRUE)
    m[is.na(m[,"f_sample"]), "f_sample"] <- 0
    m[is.na(m[,"f_ref"]), "f_ref"] <- 0
    m$k <- x
    
    # return for statistical test
    return(motifs = m)
}



# Description:
# Compute motif enrichment with Fisher's exact test (fet) based on data 
# collected by function get_motifs_v2
get_motif_enrichment_fet_v2 <- function(x) {
    # Description of parameters used in hypergeometric test (below)
    # f_sample = x[1]
    # f_ref = x[2]
    # n_sample = x[3]
    # n_ref = x[4]
    #
    # q = f_sample, 
    # m = f_ref+f_sample, 
    # n = n_ref+n_sample-(f_ref+f_sample),
    # k = n_ref+n_sample

    # ove TODO: check how is this done in gliph2
    ove <- (x[1]/x[3])/((x[2]/x[4]))
    
    q <- x[1]
    m <- x[1]+x[2]
    k <- x[3]+x[4]
    n <- k-m
    
    # TODO: do some testing on demo data 
    p <- stats::phyper(q = q, m = m, n = n, k = k, lower.tail = FALSE)
    
    return(c(ove, p))
}



# Description:
# Given data from function get_motif_enrichment, this function adds a 
# filter column with filter = T if motif is enriched (given the input 
# standards) and filter = F if not-enriched
get_motif_filter_v2 <- function(m, 
                             ks,
                             min_p, 
                             min_ove, 
                             min_o) {
    m$filter <- FALSE
    for(i in 1:length(ks)) {
        j <- which(m$p<=min_p&
                       m$ove>=min_ove[i]&
                       m$f_sample >= min_o&
                       m$k == ks[i])
        if(length(j)!=0) {
            m$filter[j] <- TRUE
        }
    }
    return(m)
}

