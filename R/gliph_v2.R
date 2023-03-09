# Description:
# This is rewritten code for GLIPH v2 (details are missing, TODO)
gliph_v2 <- function(data_sample,
                     data_ref,
                     ks = c(2, 3, 4),
                     cores = 1,
                     control = NULL) {
    
    # 0. parameter check
    
    
    # a. control check
    control <- get_control(control_in = control)
    
    # b. filter sample data based on inputs (
    # e.g. remove CDR3 with small size, 
    #       remove CDR3s without C and F, ...,)
    
    # get chains to be analyzed
    chains <- get_chains(colnames(data_sample))
    
    
    # gliph 1 and 2 would only consider non-redundant sequences
    # TODO: implement this here -> this will affect present edge
    # computation
    
    # run analysis for each chain (if available)
    clust <- vector(mode = "list", length = length(chains))
    names(clust) <- chains
    edges <- vector(mode = "list", length = length(chains))
    names(edges) <- chains
    
    for(chain in chains) {
        
        if(control$trim_flanks) {
            # NAs ignored by qqgram, how about global dist? TODO
            data_sample[, chain] <- get_trimmed_flanks(
                x = data_sample[, chain],
                flank_size = control$flank_size)
            data_ref[, chain] <- get_trimmed_flanks(
                x = data_ref[, chain],
                flank_size = control$flank_size)
        }
        
        # run local + global clustering
        clust[[chain]] <- get_chain_run_v2(
            cdr3 = data_sample[, chain],
            cdr3_ref = data_ref[, chain],
            ks = ks, 
            cores = cores, 
            control = control)
        
        # create edges of a graph
        edges[[chain]] <- get_edges(
            local_pairs = runs[[chain]]$local_pairs, 
            global_pairs = runs[[chain]]$global_pairs, 
            data_sample = data_sample,
            chain = chain)
    }
    
    return(list(clust = clust, 
                edges = edges))
}

