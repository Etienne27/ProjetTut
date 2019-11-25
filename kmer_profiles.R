


####################################################################################################################################################
####################################################### k-mer profiles #############################################################################



################################### Load embed dataset #####################################

load('D:/Documents/M2MIASHS/ProjetTut/format-tataru-dataset/output/dataset_embed.Rdata')



############################# Function for computing kmer profiles ##########################

compute_kmer_profiles <- function(ASV_seq, k=3) {
  
  start_time <- Sys.time()
  
  if (!require(pbapply)) {
    writeLines('\nInstalling required package "pbapply". \n')
    install.packages("pbapply")
  }
  if (!require(Matrix)) {
    writeLines('\nInstalling required package "Matrix". \n')
    install.packages("Matrix")
  }
  if (!require(Rcpp)) {
    writeLines('\nInstalling required package "Rcpp". \n')
    install.packages("Rcpp")
  }
  
  cppFunction('StringVector kmerC(std::string DNAseq, int k) {
    int n = DNAseq.length();
    StringVector kmer(n-k+1);
    for( int i=0; i < n-k+1; i++ ){
      kmer[i] = DNAseq.substr(i, k);
    }
    return kmer;
  }')
  
  writeLines('\nComputing k-mer sub-sequences: \n')
  
  pboptions(char='=')
  kmers <- t(pbsapply(asv.seq, kmerC, k))
  kmers_unique <- unname(unique(c(kmers)))
  
  writeLines(paste0('\nFound ', length(kmers_unique), ' unique sub-sequences of length ', k, ' in all the ASV sequences. \n'))
  writeLines('Computing sub-sequences frequency across all ASV sequences: \n')
  
  kmer_profiles <- matrix(nrow=nrow(kmers), ncol=length(kmers_unique))
  pb <- txtProgressBar(0, nrow(kmers), style = 3)
  for (i in 1:nrow(kmers)) {
    kmer_profiles[i,] <- as.numeric(table(factor(kmers[i,],lev=kmers_unique)))
    setTxtProgressBar(pb, i)
  }
  
  colnames(kmer_profiles) <- kmers_unique
  rownames(kmer_profiles) <- rownames(kmers)
  kmer_profiles <- kmer_profiles[, sort(colnames(kmer_profiles))]
  
  end_time <- Sys.time()
  writeLines('\n')
  writeLines(paste0('Total time taken: ', round(end_time - start_time, 4), 'secs.'))
  
  return(Matrix(kmer_profiles))
}


########################################### Usage ##################################################

kmer_profiles <- compute_kmer_profiles(asv.seq, k=5) ### Return sparse matrix
kmer_profiles_mat <- as.matrix(kmer_profiles)

#write.csv(kmer_profiles_mat, file='D:/Documents/M2MIASHS/ProjetTut/format-tataru-dataset/outputCSV/kmers.csv')