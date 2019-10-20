#' Get chunked sets of row-wise or column-wise indices of a matrix 
#' 
#' @name getMatrixBlocks
#'
#' @param mat Input matrix
#' @param chunk.size The size of the chunks to use for coercion
#' @param by.row Whether to chunk in a row-wise fashion
#' @param by.col Whether to chunk in a column-wise fashion
#'
#' @return A set of chunked indices
#' @export
#'
#' @examples
#' #make a sparse binary matrix
#' library(Matrix)
#' m <- 100
#' n <- 1000
#' mat <- round(matrix(runif(m*n), m, n))
#' mat.sparse <- Matrix(mat, sparse = TRUE)
#' 
#' #get row-wise chunks of 10
#' chunks <- getMatrixBlocks(mat.sparse, chunk.size = 10)

getMatrixBlocks <- function(mat, chunk.size = 1e5,
                            by.row = TRUE, by.col = FALSE) {
  message("Using chunk size: ", chunk.size)
  if (by.row) {
    message("Breaking into row chunks.")
    return(split(1:nrow(mat), ceiling(seq_along(1:nrow(mat))/chunk.size)))
  }
  
  #assumes column-wise chunking
  message("Breaking into column chunks.")
  return(split(1:ncol(mat), ceiling(seq_along(1:ncol(mat))/chunk.size)))
}

#' Convert a sparse matrix to a dense matrix in a block-wise fashion 
#' 
#' @name sparseToDenseMatrix
#'
#' @param mat Input sparse matrix
#' @param blockwise Whether to do the coercion in a block-wise manner
#' @param by.row Whether to chunk in a row-wise fashion
#' @param by.col Whether to chunk in a column-wise fashion
#' @param chunk.size The size of the chunks to use for coercion
#' @param parallel Whether to perform the coercion in parallel
#' @param cores The number of cores to use in the parallel coercion
#'
#' @return A dense matrix of the same dimensions as the input
#' 
#' @import Matrix
#' @import parallel
#' 
#' @export 
#' 
#' @examples
#' #make a sparse binary matrix
#' library(Matrix)
#' m <- 100
#' n <- 1000
#' mat <- round(matrix(runif(m*n), m, n))
#' mat.sparse <- Matrix(mat, sparse = TRUE)
#' 
#' #coerce back
#' mat.dense <- sparseToDenseMatrix(mat.sparse, chunk.size = 10)
#' 
#' #make sure they are the same dimensions
#' dim(mat) == dim(mat.dense)
#' 
#' #make sure they are the same numerically
#' all(mat == mat.dense)

sparseToDenseMatrix <- function(mat, blockwise = TRUE,
                                by.row = TRUE, by.col = FALSE,
                                chunk.size = 1e5, parallel = FALSE,
                                cores = 2) {
  if (isFALSE(blockwise)) return(as(mat, "matrix"))
  
  #do block-wise reconstruction of matrix
  chunks <- getMatrixBlocks(mat, chunk.size = chunk.size,
                            by.row = by.row, by.col = by.col)
  
  if (by.row & parallel) {
    return(do.call("rbind", mclapply(chunks, function(r) {
      return(as(mat[r,], "matrix"))
    }, mc.cores = cores)))
  }
  
  if (by.row & !parallel) {
    return(do.call("rbind", lapply(chunks, function(r) {
      return(as(mat[r,], "matrix"))
    })))
  }
  
  #assumes column-wise conversion
  if (by.col & parallel) {
    return(do.call("cbind", mclapply(chunks, function(r) {
      return(as(mat[,r], "matrix"))
    }, mc.cores = cores)))
  }
  
  return(do.call("cbind", lapply(chunks, function(r) {
    return(as(mat[,r], "matrix"))
  })))
  
}