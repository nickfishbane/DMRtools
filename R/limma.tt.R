#' Generate limma top-table for differential methylation
#'
#' Creates a custom output for testing differential methylation of M-values
#' according to a variable within a linear model for each CpG site in a 450K array.
#'
#' \code{limma.tt()} Runs limma on M-values based on a linear model including a variable of interest (for inference) and covariates.
#' \code{limma.tt.paired()} Is a custom function which takes in an array of M-values of paired samples. In this analysis, the variable of interest is a two-level factor and differences in M-values are calculated per subject according to this factor. Covariates are included in the linear model as either differences according to the variable of interest or fixed covariates. The inference is then done on the intercept in the linear model to test for differences according to the variable of interest per subject, adjusting for appropriate covariates.
#'
#' @param mval.mat Matrix of M-values with CpG sites in rows and samples in columns
#' @param betas.mat Matrix of \eqn{\beta}-values corresponding to \code{mval.mat}
#' @param pheno Data frame of phenotypes describing samples in same order as columns of \code{mval.mat} and \code{betas.mat}
#' @param subj.fld Character value for the ID of subjects in the experiment. Must correspond to a column name of \code{pheno}. For \code{limma.tt.paired}, each different subject must have exactly two samples in the array, one for each level of \code{var.interest}.
#' @param var.interest Character value for one column name in \code{pheno}. For \code{limma.tt}, this is simply the variable in the linear model on which limma will conduct inference. For \code{limma.tt.paired}, this variable must be a two-level factor for which each subject has two samples in the array, one for each factor level.
#' @param covars Character vector of covariates to include in the linear model. These should all be names of columns in \code{pheno}.
#' @param covars.diff Character vector of covariates which vary according to the levels of \code{var.interest} in the paired analysis. The differences of these covariates will be added to the linear model. Must correspond to columns of \code{pheno}.
#' @param covars.fixed Character vector of covariates which do not vary according to the levels of \code{var.interest}. These will be added as is to the linear model in the paired analysis. Must correspond to columns of \code{pheno}.
#' @param fDat Data.frame of CpG site annotations. Must contain columns named "GeneSymbol" and "cpgSite".
#'
#' @return A data.frame giving the standard output from a limma analysis for the M-values, as well as a column denoting effect sizes for corresponding \eqn{\beta}-values, called \code{betaFC}, and a column of corresponding genes for the CpG sites.
#'
#' @export
#' @examples
#'

limma.tt <- function(mval.mat, betas.mat, pheno, var.interest, covars = NULL, fDat){
  library(limma)
  stopifnot(all(c(covars, var.interest) %in% names(pheno)),
            all(colnames(betas.mat) == colnames(mval.mat)),
            all(c('GeneSymbol', 'cpgSite') %in% names(fDat)))
  form <- paste('~', var.interest)
  if(!is.null(covars))
    form <- paste(form, '+', paste(covars, collapse = ' + '))
  des <- model.matrix(as.formula(form), data = pheno)
  cat('Limma on m-values... \n')
  fit <- lmFit(mval.mat, des)
  eb.fit <- eBayes(fit)
  tt <- topTable(eb.fit, coef = 2, number = nrow(mval.mat), adjust.method = 'BH')
  tt$Gene <- fDat$GeneSymbol[match(rownames(tt), fDat$cpgSite)]
  cat('Limma on beta values... \n')
  fit.beta <- lmFit(betas.mat, des)
  eb.fit <- eBayes(fit.beta)
  tt.beta <- topTable(eb.fit, coef = 2, number = nrow(betas.mat), adjust.method = 'BH')
  tt$betaFC <- tt.beta[rownames(tt), 'logFC']
  tt
}

#' @rdname limma.tt
limma.tt.paired <- function(mval.mat, betas.mat, pheno, subj.fld, var.interest,
                            covars.diff = NULL, covars.fixed = NULL, fDat){
  library(limma)
  library(plyr)
  library(dplyr)
  stopifnot(all(c(covars.diff, covars.fixed, var.interest, subj.fld) %in% names(pheno)),
            all(colnames(betas.mat) == colnames(mval.mat)),
            all(c('GeneSymbol', 'cpgSite') %in% names(fDat)),
            ncol(mval.mat) == nrow(pheno))
  id.fld <- names(pheno)[which(laply(pheno, function(x) all(x == colnames(betas.mat))))[1]]
  phen <- pheno[c(id.fld, subj.fld, var.interest, covars.diff, covars.fixed)]
  names(phen)[1:3] <- c('id', 'subj', 'var.int')
  if(!is.null(covars.diff)){
    dups <- phen %>% group_by(subj) %>% dplyr::select(-one_of(covars.diff)) %>%
      summarise(n.samp = n_distinct(var.int))
  } else dups <- phen %>% group_by(subj) %>% summarise(n.samp = n_distinct(var.int))
  stopifnot(all(dups$n.samp == 2))
  levs <- sort(unique(phen$var.int))
  mval.diff <- t(daply(phen, .(subj), function(x){
    samps <- x$id[match(levs, x$var.int)]
    mval.mat[,samps[2]] - mval.mat[,samps[1]]}))
  beta.diff <- t(daply(phen, .(subj), function(x){
    samps <- x$id[match(levs, x$var.int)]
    betas.mat[,samps[2]] - betas.mat[,samps[1]]}))

  ## create difference for each covars.diff:
  for(v in covars.diff){
    tmp <- phen[c('subj', 'var.int', v)]; names(tmp)[3] <- 'v'
    tmp <- tmp %>% group_by(subj) %>% arrange(var.int) %>% summarise(v = diff(v))
    phen <- phen %>% dplyr::select(-one_of(v)) %>% merge(tmp, all.x = T)
    names(phen)[ncol(phen)] <- paste0(v, '_Diff')
  }
  covars <- setdiff(names(phen), c('subj', 'id', 'var.int'))
  if(length(covars) == 0) covars <- NULL

  phen <- phen %>% tbl_df %>% dplyr::select(-id, -var.int) %>% group_by(subj) %>%
    distinct() %>% arrange(subj)
  form <- '~ 1'
  if(!is.null(covars))
    form <- paste(form, '+', paste(covars, collapse = ' + '))
  des <- model.matrix(as.formula(form), data = phen)
  stopifnot(all(phen$subj == colnames(beta.diff)),
            all(phen$subj == colnames(mval.diff)))
  cat('\n', 'Limma on m-values...')
  fit <- lmFit(mval.diff, des)
  eb.fit <- eBayes(fit)
  tt <- topTable(eb.fit, coef = 1, number = nrow(mval.diff), adjust.method = 'BH')
  cat('\n Re-annotating gene names...')
  tt$Gene <- fDat$GeneSymbol[match(rownames(tt), fDat$cpgSite)]
  cat('\n Limma on beta values... \n')
  fit.beta <- lmFit(beta.diff, des)
  eb.fit <- eBayes(fit.beta)
  tt.beta <- topTable(eb.fit, coef = 1, number = nrow(beta.diff), adjust.method = 'BH')
  tt$betaFC <- tt.beta[rownames(tt), 'logFC']
  tt
}
