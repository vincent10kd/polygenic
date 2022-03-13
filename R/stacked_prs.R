#' @title Stacked P+T polygenic risk scores by regularized regression
#'
#' @description Uses individual-level data to combine ('stack') polygenic risk scores based on
#' different tuning parameter combinations to optimize predictive performance. Based on
#' Privé et al. AJHG 2019.
#'
#' @param prs_grid The element 'result_df' from the list output by prs_grid().
#' @param y The observed outcome (trait).
#' @param nfolds Number of folds to use in k-fold cross-validation within cv.glmnet().
#' @param family The distribution family of the outcome, e.g. binomial, gaussian.
#' @param plot_error If TRUE, plots the error of the stacked PRS by the tuning parameter lambda.
#' @param seed If not NULL, sets a seed prior to running cv.glmnet().
#' @param ... Additional arguments to pass to cv.glmnet().
#'
#' @return A list containing the model, the stacked PRS, and the non-zero model coefficients.
#' @export
#' @importFrom glmnet "cv.glmnet"

stacked_prs <- function(prs_grid, y, nfolds=10,
                        family='binomial',
                        plot_error=TRUE,
                        seed=NULL,
                        ...){
  X <- as.matrix(prs_grid)
  if(!is.null(seed)) set.seed(seed)
  gmod <- glmnet::cv.glmnet(X, y, nfolds=nfolds,
                            family = family,
                            ...)
  if(plot_error==TRUE){plot(gmod)}
  lasso_coefs <- coef(gmod$glmnet.fit)[,which(gmod$lambda==gmod$lambda.min)]
  print(lasso_coefs[lasso_coefs!=0])
  reslist <- list(model=gmod,
                  stacked_prs=predict(gmod,newx=X,s='lambda.min'),
                  model_coefs=lasso_coefs[lasso_coefs!=0])
  return(reslist)
}
