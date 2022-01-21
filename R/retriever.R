#' @export retriever
#' @title Generate robust disease-specific drug-response profiles
#' @description Generate disease-specific drug-response profiles that are independent of time, concentration, and cell-line. Based on the cell lines used as surrogates, the returned profiles represent the unique transcriptional changes induced by a compound in a given disease.
#' @details In the first step, we take the response profile of a given cell line to the same compound under the same concentration at different time points and averaged them. Then, the descriptive power of the generated profile to represent the drug response at a given concentration in the cell line, independently of time, is evaluated using Spearman's correlation coefficient. If the correlation with the generated profile is larger than 0.6, then the averaged profile is returned. The original drug response profiles that do not reach the threshold are removed, and the averaged profile is recomputed using the ones above the threshold, this procedure ensures the removal of aberrant or insufficient cellular responses. Only averaged signatures of at least two profiles are used in the second step.
#'
#' In the second step, we take the stable time-independent signature profiles of the compounds at a particular concentration in the same cell line. To remove the concentration dependency, we applied the same procedure described in the first step over the averaged profiles. The profiles returned by the second step are the stable ones independent of the time and concentration under the same cell line.
#'
#' In the last step, to generate disease-specific drug response profiles we again applied the procedure described in the step 1 to the stable response profiles to the same compounds in the three cell lines used as surrogate of the triple-negative breast cancer in the LINCS-L1000 project. The profiles returned by the third step are robust disease-specific transcriptional signatures representing the changes that are unique for a compound in triple-negative breast cancer.
#' @param cellLines A vector of cell lines from where the disease-specific drug response profiles will be generated.
#' @param corThreshold A numeric value between 0 and 1 that represent the minimum Spearman correlation coefficient value required by a profile to remain in the analysis.
#' @return Robust disease-specific drug-response profiles that represent the unique transcriptional changes induced by a compound in a given disease.
#' @author Daniel Osorio <daniecos@uio.no>
#' @examples
#' \donttest{
#'
#' # Generate a robust profiles across different breast cancer cell lines.
#' BRCA <- retriever(cellLines = c('MDAMB231', 'MCF7', 'SKBR3', 'HS578T', 'BT20'))
#' }

retriever <- function(cellLines, corThreshold = 0.6){

  # Loading data and getting associated metadata
  if(!requireNamespace("ccdata", quietly = TRUE)){
    message("Bioconductor package \"ccdata\" needed for this function to work. Please install it.")
  } else {
    utils::data('l1000_es', package = 'ccdata', envir = environment())
    l1000_es <- l1000_es
    colnames(l1000_es) <- gsub('\\[|\\]', '', colnames(l1000_es))
    profilesMetadata <- getProfilesMetadata(colnames(l1000_es))
    genesMetadata <- rownames(l1000_es)

    # Filtering cell lines
    mCellLines <- match.arg(cellLines, unique(profilesMetadata$cellLine), several.ok = TRUE)
    if(!all(cellLines %in% mCellLines)){
      stop(message = paste0('\n', cellLines[!cellLines %in% mCellLines], ' cell line not found'))
    }
    profilesMetadata <- profilesMetadata[profilesMetadata$cellLine %in% mCellLines,]

    # Step 1
    S1 <- profilesMetadata[,1:3]
    S1 <- apply(S1, 1, function(X){paste0(X, collapse = '_')})
    S1 <- unique(S1)
    message('Step 1: ')
    S1Profiles <- pbapply::pbsapply(S1, function(N){
      X <- l1000_es[,profilesMetadata$name[grepl(N, profilesMetadata$name)], drop = FALSE]
      if(ncol(X)>1){
        X <- preprocessCore::normalize.quantiles(as.matrix(X))
      }
      corValue <- stats::cor(data.frame(rowMeans(X),X),method = 'sp')[,1]
      corValue <- (corValue[-1] > corThreshold)
      if(sum(corValue) > 1){
        return(rowMeans(X[,corValue]))
      } else {
        return(rep(NA, 1001))
      }
    })
    colnames(S1Profiles) <- S1
    S1Profiles <- S1Profiles[,stats::complete.cases(t(S1Profiles))]
    S1_names <- colnames(S1Profiles)
    rownames(S1Profiles) <- genesMetadata

    # Step 2
    S2 <- profilesMetadata[,1:2]
    S2 <- apply(S2, 1, function(X){paste0(X, collapse = '_')})
    S2 <- unique(S2)
    message('Step 2: ')
    S2Profiles <- pbapply::pbsapply(S2, function(N){
      X <- S1Profiles[,grepl(N, S1_names, fixed = TRUE), drop = FALSE]
      if(ncol(X)>1){
        X <- preprocessCore::normalize.quantiles(as.matrix(X))
      }
      corValue <- stats::cor(data.frame(rowMeans(X),X),method = 'sp')[,1]
      corValue <- (corValue[-1] > corThreshold)
      if(sum(corValue) > 1){
        return(rowMeans(X[,corValue]))
      } else {
        return(rep(NA, 1001))
      }
    })
    colnames(S2Profiles) <- S2
    S2Profiles <- S2Profiles[,stats::complete.cases(t(S2Profiles))]
    S2_names <- colnames(S2Profiles)
    rownames(S2Profiles) <- genesMetadata

    # Step 3
    S3 <- unique(profilesMetadata[,1])
    message('Step 3:')
    S3Profiles <- pbapply::pbsapply(S3, function(N){
      X <- S2Profiles[,grepl(N, S2_names, fixed = TRUE), drop = FALSE]
      if(ncol(X)>1){
        X <- preprocessCore::normalize.quantiles(as.matrix(X))
      }
      corValue <- stats::cor(data.frame(rowMeans(X),X),method = 'sp')[,1]
      corValue <- (corValue[-1] > corThreshold)
      if(sum(corValue) > 1){
        return(rowMeans(X[,corValue]))
      } else {
        return(rep(NA, 1001))
      }
    })
    colnames(S3Profiles) <- S3
    S3Profiles <- S3Profiles[,stats::complete.cases(t(S3Profiles))]
    S3_names <- colnames(S3Profiles)
    rownames(S3Profiles) <- rownames(l1000_es)

    return(S3Profiles)
  }
}
