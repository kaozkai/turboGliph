#' Check turboGliph Input parameters
#'
#' @param ... turboGliph parameters
#'
#' @return Returns checked / corrected parameters
#' @export
#'
#' @examples
parameter_check <- function(...){
  
  p <- list(...)
  
  ### result_folder
  if(!base::is.character(p$result_folder)) base::stop("result_folder has to be a character object")
  if(base::length(p$result_folder) > 1) base::stop("result_folder has to be a single path")
  p$save_results <- FALSE
  if(p$result_folder != ""){
    if(base::substr(p$result_folder,base::nchar(p$result_folder),base::nchar(p$result_folder)) != "/") p$result_folder <- base::paste0(p$result_folder,"/")
    if(!base::dir.exists(p$result_folder)) base::dir.create(p$result_folder)
    p$save_results <- TRUE
    if(base::file.exists(base::paste0(p$result_folder,"kmer_resample_",p$sim_depth,"_log.txt"))){
      p$save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(p$result_folder,"kmer_resample_",p$sim_depth,"_log.txt"),"\n"))
    }
    if(base::file.exists(base::paste0(p$result_folder,"kmer_resample_",p$sim_depth,"_minp",p$lcminp,"_ove",base::paste(p$lcminove, collapse = "_"),".txt"))){
      p$save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(p$result_folder,"kmer_resample_",p$sim_depth,"_minp",p$lcminp,"_ove",base::paste(p$lcminove, collapse = "_"),".txt"),"\n"))
    }
    if(base::file.exists(base::paste0(p$result_folder,"kmer_resample_",p$sim_depth,"_all_motifs.txt"))){
      p$save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(p$result_folder,"kmer_resample_",p$sim_depth,"_all_motifs.txt"),"\n"))
    }
    if(base::file.exists(base::paste0(p$result_folder, "clone_network.txt"))){
      p$save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(p$result_folder, "clone_network.txt"),"\n"))
    }
    if(base::file.exists(base::paste0(p$result_folder,"convergence_groups.txt"))){
      p$save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(p$result_folder,"convergence_groups.txt"),"\n"))
    }
    if(base::file.exists(base::paste0(p$result_folder,"parameter.txt"))){
      p$save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(p$result_folder,"parameter.txt"),"\n"))
    }
  }
  
  ### refdb_beta
  # if(!(refdb_beta %in% base::c("gliph_reference", "human_v1.0_CD4", "human_v1.0_CD8", "human_v1.0_CD48", "human_v2.0_CD4",
  #                        "human_v2.0_CD8", "human_v2.0_CD48", "mouse_v1.0_CD4", "mouse_v1.0_CD8", "mouse_v1.0_CD48")) &&
  #    !base::is.data.frame(refdb_beta)){
  if(!base::is.data.frame(p$refdb_beta)){
    if(base::length(p$refdb_beta) != 1 || !is.character(p$refdb_beta)){
      base::stop("refdb_beta has to be a data frame (containing CDR3b sequences in the first column and optional V-gene information in the second column) or the value 'gliph_reference'")
    } else if(!(p$refdb_beta %in% base::c("gliph_reference"))){
      base::stop("refdb_beta has to be a data frame (containing CDR3b sequences in the first column and optional V-gene information in the second column) or the value 'gliph_reference'")
    }
  }
  
  ### v_usage_freq
  if(!base::is.null(p$v_usage_freq)){
    if(base::is.data.frame(p$v_usage_freq)){
      if(base::ncol(p$v_usage_freq) < 2) base::stop("v_usage_freq has to be a data frame containing V-gene information in the first column and the corresponding frequency in a naive  T-cell repertoire in the second column.")
      if(base::nrow(p$v_usage_freq) < 1) base::stop("v_usage_freq has to contain at least one row.")
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(p$v_usage_freq[,2])))) == TRUE){
        base::stop("The second column of v_usage_freq must contain the frequency of the corresponding V-gene in the first column in a naive T-cell repertoire.")
      } else p$v_usage_freq[,2] <- as.numeric(p$v_usage_freq[,2])
      
    } else {base::stop("v_usage_freq has to be a data frame containing V-gene information in the first column and the corresponding frequency in a naive T-cell repertoire in the second column.")}
  }
  
  ### cdr3_length_freq
  if(!base::is.null(p$cdr3_length_freq)){
    if(base::is.data.frame(p$cdr3_length_freq)){
      if(base::ncol(p$cdr3_length_freq) < 2) base::stop("cdr3_length_freq has to be a data frame containing CDR3 lengths in the first column and the corresponding frequency in a naive  T-cell repertoire in the second column.")
      if(base::nrow(p$cdr3_length_freq) < 1) base::stop("cdr3_length_freq has to contain at least one row.")
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(p$cdr3_length_freq[,2])))) == TRUE){
        base::stop("The second column of cdr3_length_freq must contain the frequency of the corresponding CDR3 length in the first column in a naive T-cell repertoire.")
      } else p$cdr3_length_freq[,2] <- as.numeric(p$cdr3_length_freq[,2])
      
    } else {base::stop("cdr3_length_freq has to be a data frame containing CDR3 lengths in the first column and the corresponding frequency in a naive T-cell repertoire in the second column.")}
  }
  
  ### ref_cluster_size
  if(!(p$ref_cluster_size %in% base::c("original", "simulated") || !base::is.character(p$ref_cluster_size) || base::length(p$ref_cluster_size) > 1)){
    base::stop("ref_cluster_size has to be either 'original' or 'simulated'.")
  }
  
  ### sim_depth
  if(!base::is.numeric(p$sim_depth)) base::stop("sim_depth has to be numeric")
  if(base::length(p$sim_depth) > 1) base::stop("sim_depth has to be a single number")
  if(p$sim_depth < 1) base::stop("sim_depth must be at least 1")
  p$sim_depth <- base::round(p$sim_depth)
  
  ### lcminp
  if(!base::is.numeric(p$lcminp)) base::stop("lcminp has to be numeric")
  if(base::length(p$lcminp) > 1) base::stop("lcminp has to be a single number")
  if(p$lcminp <= 0) base::stop("lcminp must be greater than 0")
  
  ### lcminove
  if(!base::is.numeric(p$lcminove)) base::stop("lcminove has to be numeric")
  if(base::length(p$lcminove) > 1 && base::length(p$lcminove) != base::length(p$motif_length)) base::stop("lcminove has to be a single number or of same length as motif_length")
  if(base::any(p$lcminove < 1)) base::stop("lcminove must be at least 1")
  
  ### kmer_mindepth
  if(!base::is.numeric(p$kmer_mindepth)) base::stop("kmer_mindepth has to be numeric")
  if(base::length(p$kmer_mindepth) > 1) base::stop("kmer_mindepth has to be a single number")
  if(p$kmer_mindepth < 1) base::stop("kmer_mindepth must be at least 1")
  p$kmer_mindepth <- base::round(p$kmer_mindepth)
  
  ### accept_sequences_with_C_F_start_end
  if(!base::is.logical(p$accept_sequences_with_C_F_start_end)) base::stop("accept_sequences_with_C_F_start_end has to be logical")
  
  ### min_seq_length
  if(!base::is.numeric(p$min_seq_length)) base::stop("min_seq_length has to be numeric")
  if(base::length(p$min_seq_length) > 1) base::stop("min_seq_length has to be a single number")
  if(p$min_seq_length < 0) base::stop("min_seq_length must be at least 0")
  p$min_seq_length <- base::round(p$min_seq_length)
  
  ### gccutoff
  if(!base::is.null(p$gccutoff) && !base::is.numeric(p$gccutoff)) base::stop("gccutoff has to be NULL or numeric")
  if(!base::is.null(p$gccutoff) && base::length(p$gccutoff)>1) base::stop("gccutoff has to be NULL or a single number")
  if(!base::is.null(p$gccutoff) && p$gccutoff < 0) base::stop("gccutoff must be at least 0")
  
  ### structboundaries
  if(!base::is.logical(p$structboundaries)) base::stop("structboundaries has to be logical")
  
  ### boundary_size
  if(!base::is.numeric(p$boundary_size)) base::stop("boundary_size has to be numeric")
  if(base::length(p$boundary_size) > 1) base::stop("boundary_size has to be a single number")
  if(p$boundary_size < 0) base::stop("boundary_size must be at least 0")
  p$boundary_size <- base::round(p$boundary_size)
  if(p$structboundaries == TRUE) p$min_seq_length <- base::max(p$min_seq_length, ((p$boundary_size * 2)+1))
  
  ### motif_length
  if(!base::is.numeric(p$motif_length)) base::stop("motif_length has to be numeric")
  if(base::any(p$motif_length < 1)) base::stop("values of motif_length must be at least 1")
  p$motif_length <- base::round(p$motif_length)
  
  ### discontinuous
  if(!base::is.logical(p$discontinuous)) base::stop("discontinuous has to be logical")
  
  ### make_depth_fig
  if(!base::is.logical(p$make_depth_fig)) base::stop("make_depth_fig has to be logical")
  
  ### local_similarities
  if(!base::is.logical(p$local_similarities)) base::stop("local_similarities has to be logical")
  
  ### global_similarities
  if(!base::is.logical(p$global_similarities)) base::stop("global_similarities has to be logical")
  if(p$local_similarities == FALSE && p$global_similarities == FALSE) base::stop("Either local_similarities or global_similarities have to be TRUE")
  
  ### global_vgene
  if(!base::is.logical(p$global_vgene)) base::stop("global_vgene has to be logical")
  
  ### positional_motifs
  if(!base::is.logical(p$positional_motifs)) base::stop("positional_motifs has to be logical")
  
  ### cdr3_len_stratify
  if(!base::is.logical(p$cdr3_len_stratify)) base::stop("cdr3_len_stratify has to be logical")
  
  ### vgene_stratify
  if(!base::is.logical(p$vgene_stratify)) base::stop("vgene_stratify has to be logical")
  
  ### public_tcrs
  if(!base::is.logical(p$public_tcrs)) base::stop("public_tcrs has to be logical")
  
  ### cluster_min_size
  if(!base::is.numeric(p$cluster_min_size)) base::stop("cluster_min_size has to be numeric")
  if(base::length(p$cluster_min_size) > 1) base::stop("cluster_min_size has to be a single number")
  if(p$cluster_min_size < 1) base::stop("cluster_min_size must be at least 1")
  p$cluster_min_size <- base::round(p$cluster_min_size)
  
  ### hla_cutoff
  if(!base::is.numeric(p$hla_cutoff)) base::stop("hla_cutoff has to be numeric")
  if(base::length(p$hla_cutoff) > 1) base::stop("hla_cutoff has to be a single number")
  if(p$hla_cutoff > 1 || p$hla_cutoff < 0) base::stop("hla_cutoff must be between 0 and 1")
  
  ### n_cores
  if(!base::is.null(p$n_cores))
  {
    if(!base::is.numeric(p$n_cores)) base::stop("n_cores has to be numeric")
    if(base::length(p$n_cores) > 1) base::stop("n_cores has to be a single number")
    if(p$n_cores < 1) base::stop("n_cores must be at least 1")
    p$n_cores <- base::round(p$n_cores)
  }
  
  ### min_seq_length second round
  if(p$structboundaries == TRUE) p$min_seq_length <- base::max(p$min_seq_length, 2*p$boundary_size)
  
  return(p)
  
}