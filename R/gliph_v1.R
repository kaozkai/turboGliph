
# Description:
# This is rewritten code for GLIPH v1 (details are missing, TODO)
gliph_v1 <- function(cdr3,
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
           FUN = get_motifs_v1,
           cdr3 = cdr3, 
           cdr3_ref = cdr3_ref, 
           B = B,
           min_o = control$filtering$local_min_o,
           cores = cores)
    names(motifs) <- as.character(ks)
    
    
    # b. compute local enrichment scores
    future::plan(future::multisession, workers = cores)
    motif_enrichment <- do.call(rbind, future.apply::future_lapply(
        X = ks, 
        FUN = get_motif_enrichment_v1, 
        m = motifs, 
        B = B, 
        cores = cores,
        future.seed = TRUE))
    future::plan(future::sequential())
    
    
    # c. add filter flag
    motif_enrichment <- get_motif_filter_v1(
        ks = ks,
        m = motif_enrichment, 
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
# Setup control list.
get_control <- function(control_in) {
    
    control <- list(
        refdb = list(cdr3b = NULL, # what is the default?
                     v_usage_freq = NULL,
                     cdr3_length_freq = NULL),
        input = list(min_seq_length = 8,
                     cdr3b_start_end_with_C_F = TRUE, #accept_sequences_with_C_F_start_end = TRUE,
                     structboundaries = TRUE,
                     boundary_size = 3),
        clustering = list(global_vgene = FALSE,
                          positional_motifs = FALSE,
                          public_tcrs = TRUE),
        filtering = list(global_max_dist = 1,
                         local_min_p = 0.05, 
                         local_min_ove = c(10^3, 10^2, 10^1), 
                         local_min_o = 3),
        sampling = list(cdr3_len_stratify = FALSE,
                        vgene_stratify = FALSE),
        extra = list(make_depth_fig = FALSE),
        scoring = list(cluster_min_size = 2,
                       hla_cutoff = 0.1,
                       ref_cluster_size = "original"))
    
    if(missing(control_in)|is.null(control_in)) {
        return(control)
    }
    
    # here run check on individual parameters in control
    check_control(control_in)
    
    # if failed check -> stop
    # else return control
    
    
    # edit control by user-defined control_in
    # There are packages to update list given another list, but here we can 
    # live with the following "inefficiency" as the list is generally small 
    # (~10 elements)
    for(i in 1:length(control_in)) {
        n_i <- names(control_in)[i]
        for(j in 1:length(control_in[[i]])) {
            n_j <- names(control_in[[i]])[j]
            control[[n_i]][[n_j]] <- control_in[[n_i]][[n_j]]
        }
    }
    
    return(control)
}



# Description:
# Check parameters in list control.
check_control <- function(control) {
    
    
}



# Description:
# Computes motif frequencies for a sample and reference
# uses functions:
# * get_kmers_freq_ref
# * get_kmers_freq_sample
get_motifs_v1 <- function(cdr3, 
                       cdr3_ref, 
                       B, 
                       ks, 
                       cores, 
                       min_o) {
    
    # Description:
    # x = k in k-mer
    get_kmers_freq_ref <- function(x, 
                                   cdr3, 
                                   N, 
                                   B, 
                                   relevant_motifs, 
                                   cores){
        
        get_qgrams <- function(x, q, cdr3, N, relevant_motifs) {
            draw_cdr3 <- sample(x = cdr3, size = N, replace = TRUE)
            o <- stringdist::qgrams(draw_cdr3, q = q)
            if(ncol(o)==0) {
                return(NA)
            }
            o <- o[1,]
            o <- o[names(o) %in% relevant_motifs]
            if(length(o)==0) {
                return(NA)
            }
            return(o)
        }
        
        future::plan(future::multisession, workers = cores)
        o <- future.apply::future_lapply(
            X = 1:B, 
            q = x,
            N = N, 
            relevant_motifs = relevant_motifs, 
            cdr3 = cdr3,
            FUN = get_qgrams,
            future.seed = TRUE)
        future::plan(future::sequential())
        return(o)
    }
    
    
    # Description:
    # x = k in k-mer
    get_kmers_freq_sample <- function(x, 
                                      cdr3, 
                                      min_o){
        
        o <- stringdist::qgrams(cdr3, q = x)
        if(ncol(o)==0) {
            return(NA)
        }
        o <- o[1, o[1,]>=min_o]
        if(length(o)==0) {
            return(NA)
        }
        return(o)
        # return(data.frame(motif = colnames(o), 
        #                   f = o[1,]))
    }
    
    # find motifs in sample
    motif_sample <- lapply(
        X = ks, # edit here
        FUN = get_kmers_freq_sample, 
        cdr3 = cdr3, # edit here 
        min_o = min_o)
    names(motif_sample) <- ks
    
    
    # do sampling & find motifs
    found_kmers <- as.vector(unlist(
        lapply(X = motif_sample, FUN = names)))
    motif_ref <- lapply(
        X = ks,
        FUN = get_kmers_freq_ref, 
        cdr3 = cdr3_ref,
        B = B,
        N = length(cdr3),
        relevant_motifs = found_kmers,
        cores = cores)
    names(motif_ref) <- ks
    
    return(list(motif_sample = motif_sample,
                motif_ref = motif_ref))
}


# Description:
# Computes motif enrichment with data collected by function:
# * get_motifs
get_motif_enrichment_v1 <- function(x, 
                                 m, 
                                 B, 
                                 cores) {
    
    get_e <- function(x, B) {
        # o
        o <- x[1]
        
        # e
        e <- x[-1]
        
        # mean, max
        e_mean <- base::mean(e)
        e_max <- base::max(e)
        
        # OvE = /e
        ove <- o/e_mean
        
        # prob e>=o
        p <- sum(e>=o)/B
        
        # return
        return(c(e_mean, e_max, o, ove, p))
    }
    
    motif_sample <- m[[as.character(x)]]$motif_sample[[1]]
    motif_ref <- m[[as.character(x)]]$motif_ref[[1]]
    
    # matrix of k-mer counts
    f_m <- matrix(data = 0, ncol = B+1, nrow = length(motif_sample))
    rownames(f_m) <- names(motif_sample)
    f_m[names(motif_sample), 1] <- motif_sample
    for(i in 1:length(motif_ref)) {
        s <- motif_ref[[i]]
        f_m[names(s),i+1] <- s 
    }
    
    e <- t(apply(X = f_m, MARGIN = 1, FUN = get_e, B = B))
    
    # format output
    e <- data.frame(e)
    colnames(e) <- c("sim_mean", "sim_max", 
                     "obs", "ove", "p")
    e$motif <- rownames(e)
    e$k <- x
    return(e)
}


# Description:
# Given data from function get_motif_enrichment, this function adds a 
# filter column with filter = T if motif is enriched (given the input 
# standards) and filter = F if not-enriched
get_motif_filter_v1 <- function(m, 
                             ks, 
                             min_p, 
                             min_ove, 
                             min_o) {
    m$filter <- FALSE
    for(i in 1:length(ks)) {
        j <- which(m$p<=min_p & m$ove >= min_ove[i] & 
                       m$obs >= min_o & m$k == ks[i])
        if(length(j)!=0) {
            m$filter[j] <- TRUE
        }
    }
    return(m)
}




# Description:
# Look for enriched motifs in CDR3s. CDR3 sequence pairs of local 
# clusters are returned. Used by gliph_v1 and gliph_v2.
get_local_pair <- function(cdr3, 
                           motif) {
    # if no enriched motifs
    if(length(motif)==0) {
        return(NULL)
    }
    
    get_motif_in_cdr <- function(x, motif, cdr3) {
        j <- which(regexpr(pattern = motif[x], text = cdr3)!=-1)
        if(length(j)==1) {
            return(data.frame(from = j, to = j, motif = motif[x]))
        }
        u <- t(utils::combn(x = j, m = 2))
        u <- rbind(u, cbind(j,j))
        return(data.frame(from = u[,1], to = u[,2], motif = motif[x]))
    }
    
    return(do.call(rbind, lapply(X = 1:length(motif), 
                                 motif = motif, 
                                 FUN = get_motif_in_cdr, 
                                 cdr3 = cdr3)))
}


# Description:
# Look for pairs of global connections. Used by gliph_v1 and gliph_v2.
get_global_pairs <- function(cdr3, 
                             global_max_dist) {
    # Jan is looping over all sequences and computing distances with the rest.
    # two problems:
    # a) slower than passing vector to stringdist (but lower max. memory footprint)
    # b) hamming distance between sequences with unequal lengths computed nontheless
    
    # at the very least loop over cdr3.lengths -> faster, but somewhat more 
    # memory is data passed together to stringdist
    
    cdr3_len <- base::nchar(cdr3)
    cdr3_lens <- unique(cdr3_len)
    
    get_hamming_dist <- function(x, cdr3, cdr3_len, global_max_dist) {
        is <- which(cdr3_len == x)
        if(length(is)==1) {
            return(NULL)
        }
        if(length(is)==2) {
            d <- stringdist::stringdist(
                a = cdr3[is[1]],
                b = cdr3[is[2]],
                method = "hamming")
            if(d>global_max_dist) {
                return(NULL)
            } 
            return(c(cdr3_is[is[1]], cdr3_is[is[2]]))
        }
        
        d <- stringdist::stringdistmatrix(
            a = cdr3[is],
            b = cdr3[is],
            method = "hamming")
        d[upper.tri(x = d, diag = TRUE)] <- NA
        # d[1:nrow(d), 1:nrow(d)] <- NA
        js <- which(d<=global_max_dist, arr.ind = TRUE)
        if(nrow(js)==0) {
            return(NULL)
        }
        return(cbind(is[js[,1]], is[js[,2]]))
    }
    
    # Description:
    # same as get_hamming_dist but slower (x3), however has much 
    # smaller memory footprint -> appropriate for large input
    # Similar but faster implementation than that in Jan's code
    get_hamming_dist_memory <- function(x, cdr3, cdr3_len, global_max_dist) {
        is <- which(cdr3_len == x)
        if(length(is)==1) {
            return(NULL)
        }
        
        get_pairdist  <- function(x, a, len_a, global_max_dist) {
            d <- stringdist::stringdist(a = a[x], 
                                        b = a[(x+1):len_a], 
                                        method = "hamming")
            js <- which(d<=global_max_dist)
            if(length(js)==0) {
                return(NULL)
            }
            js <- x+js
            return(cbind(rep(x = x, times = length(js)), js))
        }
        
        hd <- lapply(X = 1:(length(is)-1),
                     FUN = get_pairdist, 
                     a = cdr3[is],
                     len_a <- length(is),
                     global_max_dist = global_max_dist)
        hd <- do.call(rbind, hd)
        if(is.null(hd)) {
            return(hd)
        }
        # map to original indices
        return(cbind(is[hd[,1]], is[hd[,2]]))
    }
    
    
    hd <- lapply(X = cdr3_lens, 
                 FUN = get_hamming_dist, 
                 cdr3 = cdr3, 
                 cdr3_len = cdr3_len, 
                 global_max_dist = global_max_dist)
    hd <- do.call(rbind, hd)
    return(hd)
    # microbenchmark(
    #     times = 10,
    #     "a" = {
    #         hd <- lapply(X = cdr3_lens,
    #                      FUN = get_hamming_dist,
    #                      cdr3 = cdr3,
    #                      cdr3_len = cdr3_len,
    #                      global_max_dist = global_max_dist)
    #     },
    #     "b" = {
    #         lapply(X = cdr3_lens,
    #                FUN = get_hamming_dist_memory,
    #                cdr3 = cdr3,
    #                cdr3_len = cdr3_len,
    #                global_max_dist = global_max_dist)
    #     })
    # 
    # membench <- bench::mark(
    #     lapply(X = cdr3_lens,
    #            FUN = get_hamming_dist,
    #            cdr3 = cdr3,
    #            cdr3_len = cdr3_len,
    #            global_max_dist = global_max_dist),
    #     lapply(X = cdr3_lens,
    #            FUN = get_hamming_dist_memory,
    #            cdr3 = cdr3,
    #            cdr3_len = cdr3_len,
    #            global_max_dist = global_max_dist),
    #     future_lapply(X = cdr3_lens,
    #                   FUN = get_hamming_dist,
    #                   cdr3 = cdr3,
    #                   cdr3_len = cdr3_len,
    #                   global_max_dist = global_max_dist),
    #     future_lapply(X = cdr3_lens,
    #                   FUN = get_hamming_dist_memory,
    #                   cdr3 = cdr3,
    #                   cdr3_len = cdr3_len,
    #                   global_max_dist = global_max_dist),
    #     check = F)
}

