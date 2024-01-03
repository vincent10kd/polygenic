#' @title Obtain conditional estimates from marginal estimates and a correlation matrix
#'
#' @description Returns conditional estimates based on marginal estimates (such as public GWAS summary statistics) and a correlation matrix.
#' Allows post-hoc adjustments for covariates without access to the original data.
#' Among other things, this is useful for fine-mapping and for decorrelating variants to prevent double-counting in polygenic risk scores.
#'
#' @param marginal_coefs A numeric vector of k marginal coefficient estimates from linear or logistic regression models.
#' @param X An n*k matrix of predictors, with k equal to the length of marginal_coefs. Is not necessary when covX or corX and MAF are provided.
#' @param covX A k*k covariance matrix.
#' @param corX A k*k correlation matrix.
#' @param MAF A numeric vector of k minor allele frequencies.
#' @param N A single constant representing the sample size of the GWAS from which the marginal estimates were derived.
#' @param marginal_se A numeric vector of k marginal standard error estimates from linear or logistic regression models.
#' @param estimate_se If TRUE and if marginal_se is not NULL, returns conditional standard error estimates.
#' @param y A continuous or binary response vector. If provided together with X, standard OLS is fit to obtain the estimates.
#' @param varY The variance of the response vector. Is estimated if not provided.
#' @param ridge If TRUE, applies a ridge penalty.
#' @param lambda The parameter controlling the degree of (ridge) regularization.
#' @param binary If TRUE, assumes that marginal_coefs are log odds ratios and that the standard errors should be estimated based on a different variance approximation for the response vector.
#' @param caseprop The proportion of the response vector that represents 'cases', i.e. the proportion of the binary-valued vector that is equal to 1.
#'
#' @return A data.frame of conditional regression estimates, including standard errors and p-values depending on the input arguments.
#' @examples
#' marg_coefs <- rbind(summary(lm(Sepal.Length~Petal.Width, iris))$coef[2,1:2],
#' summary(lm(Sepal.Length~Petal.Length, iris))$coef[2,1:2],
#' summary(lm(Sepal.Length~Sepal.Width, iris))$coef[2,1:2])
#'
#' covX <- cov(iris[c('Petal.Width','Petal.Length','Sepal.Width')])

#' N <- nrow(iris)

#' marg2con(marginal_coefs = marg_coefs[,1],
#'         covX = covX,
#'         N = N,
#'         marginal_se = marg_coefs[,2],
#'         estimate_se = TRUE)

#' @export

marg2con <- function(marginal_coefs,
                     X=NULL,
                     covX=NULL,
                     corX=NULL,
                     MAF=NULL,
                     N=NULL,
                     marginal_se=NULL,
                     estimate_se=FALSE,
                     y=NULL,
                     varY=NULL,
                     ridge=FALSE,
                     lambda=0,
                     binary=FALSE,
                     caseprop=0.25){

  beta <- as.matrix(marginal_coefs)
  if(!is.null(X) & is.null(covX)){
    Xm <- scale(X, scale=FALSE)
    covX <- (t(Xm) %*% Xm) /(N-1)
  }
  if(is.null(X) & is.null(covX)){
    if(!is.null(corX) & !is.null(MAF)){
      sdX <- sqrt((2*MAF*(1-MAF)))
      covX <- diag(sdX) %*% corX %*% diag(sdX)
    }
    else stop('One of X or covX should be provided. Alternatively, corX and MAF are required.')
  }

  if(ridge){
    lambdaI <- lambda * diag(ncol(X2))
    cond <- solve(covX + lambdaI) %*% diag(diag(covX)) %*% beta
    cat('\n Ridge-regularized conditional estimates from marginal estimates and covariance matrix:\n')
    return(data.frame(beta=cond))
  }
  cond <- solve(covX, diag(diag(covX)) %*% beta)

  if(estimate_se){
    if(!is.null(y)){

      X2 <- cbind(1,X)
      beta_adj <- solve(t(X2)%*%X2, t(X2) %*% y)
      resids <- y - (X2 %*% beta_adj)
      df <- N - length(beta_adj)
      mse <- sum(resids^2) / df
      resid_se <- sqrt(mse)
      se <- sqrt(diag(mse * solve(t(X2)%*%X2)))
      t <- beta_adj/se
      p <- 2 * pt(abs(t), df, lower=F)

      res <- data.frame(beta = round(beta_adj,3),
                        se = round(se,3),
                        t = round(t,3),
                        p = p)
      rownames(res)[1] <- '(Intercept)'

      ypred <- X2 %*% beta_adj
      Rsq <- cor(y, ypred)^2
      Rsq_adj <- 1-(((1-Rsq)*(N-1))/df)

      cat('\n Conditional estimates from OLS:\n')
      print(res)
      cat('\nR-squared: ', round(Rsq,3),
          '; Adjusted R-squared: ', round(Rsq_adj,3),
          '\nMSE: ',mse,
          '; Residual SE: ', round(resid_se,3),
          '; Degrees of freedom: ',df,'\n', sep='')
      invisible(res)

    }
    else{
      if(is.null(N)){
        stop('A sample size (N) should be provided to approximate conditional SEs.')
      }
      if(!is.null(X)){
        X2 <- scale(X, scale=FALSE)
        covX <- (t(X2) %*% X2) /(N-1)
      }
      if(is.null(marginal_se)){
        stop('Marginal standard errors are required to approximate conditional SEs.')
      }
      beta_adj <- solve(covX, diag(diag(covX)) %*% marginal_coefs)

      B <- t(beta_adj) %*% diag(diag(covX)) %*% marginal_coefs
      df <- N - length(beta_adj)
      D <- diag(covX)
      est_varY <- median(D*marginal_se^2*(N-1)+D*marginal_coefs^2)
      if(binary){
        est_varY <- 1
      }
      varY <- ifelse(is.null(varY), est_varY, varY)
      mse <- varY - B
      resid_se <- sqrt(mse)
      if(binary){
        N_eff <- N*caseprop*(1-caseprop)
        N <- N_eff
        mse <- 1
      }
      se <- sqrt(mse[1]/N * diag(solve(covX)))
      t <- beta_adj/se
      p <- 2 * pt(abs(t), df, lower=F)

      res <- data.frame(beta = beta_adj,
                        se = se,
                        t = t,
                        p = p)

      if(!is.null(X)){
        ypred <- X2 %*% beta_adj
        Rsq <- B/varY
        Rsq_adj <- 1-(((1-Rsq)*(N-1))/df)
      }

      if(!binary){
        cat("\n Conditional estimates with SEs per Yang et al.'s method:\n")
      }
      else{
        cat("\n Conditional estimates with SEs, adapted for binary trait (Pirinen 2019):\n")
      }
      print(res)
      if(!binary & !is.null(X)){
        cat('\nR-squared: ', round(Rsq,3),
            '; Adjusted R-squared: ', round(Rsq_adj,3),
            '\nMSE: ',mse,
            '; Residual SE: ', round(resid_se,3),
            '; Degrees of freedom: ',df,'\n', sep='')
      }
      invisible(res)
    }
  }
  else{
    cat('\n Conditional estimates from marginal estimates and covariance matrix:\n')
    res <- data.frame(beta=cond)
    print(res)
  }
}
