


####################################################################################################################################################
####################################################### Glove embedding ############################################################################



################################### Load embed dataset #####################################

load('D:/Documents/M2MIASHS/ProjetTut/format-tataru-dataset/output/dataset_embed.Rdata')



################################### Compute co-occurence matrix ####################################

### Compute presence/occurence matrix from abundance matrix
presence_matrix <- Matrix(1*(as.matrix(X)!=0))

### Compute co-occurence matrix as t(presence_matrix) %*% presence_matrix
cooc <- crossprod(presence_matrix)
diag(cooc) <- 0 



################################ Glove algorithm #####################################

require(text2vec)

### Create new glove object
gloveModel <- GloVe$new(113, vocabulary = asv.seq, x_max = 100, learning_rate = 0.1)

### Get vectors
vectors  <- gloveModel$fit_transform(cooc, n_iter = 10)

### Get context
context <- t(gloveModel$components)

### Get taxa embedding as vectors + context as advised in http://nlp.stanford.edu/pubs/glove.pdf
taxa_vectors <- vectors + context

#write.csv(taxa_vectors, file='D:/Documents/M2MIASHS/ProjetTut/format-tataru-dataset/outputCSV/embedded.csv')
