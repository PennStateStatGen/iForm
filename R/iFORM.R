#' Interaction Screening for Ultra-High Dimensional Data
#'
#' Extended variable selection approaches to jointly model main and interaction effects from high-dimensional data orignally proposed by Hao and Zhang (2014) and extended by Gosik and Wu (2016).
#' Based on a greedy forward approach, their model can identify all possible interaction effects through two algorithms, iFORT and iFORM, which have been proved to possess sure screening property in an ultrahigh-dimensional setting.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param data data.frame of your data with the response and all p predictors
#' @param heredity  a string specifying the heredity to be considered. NULL, weak, strong
#' @param higher_order  logical TRUE indicating to include order-3 interactions in the search (default FALSE)
#' @return a summary of the linear model returned after the selection procedure
#' @author Kirk Gosik
#' @details
#' Runs the iFORM selection procedure on the dataset and returns a linear model
#' of the final selected model.  The model is of an R object of class "lm"
#' @seealso \code{lm}
#' @seealso \code{model.frame}
#' @export
#' @importFrom stats lm
#' @importFrom stats model.frame


iForm <- function(formula,
                  data,
                  heredity = "strong",
                  higher_order = FALSE) {

  dat <- model.frame(formula, data)

  y <- dat[ , 1]
  x <- dat[ , -1]
  p <- ncol(x)
  n <- nrow(x)
  C <- names(x)
  S <- NULL
  M <- NULL
  bic <- NULL

  fit <- iformselect(x, y, p, n, C, S, bic, heredity, higher_order)

  y <- fit$y
  S <- fit$S
  bic <- fit$bic

  model_formula <- as.formula(paste("y ~0+", paste(S[1:which.min(bic)], collapse = "+")))
  lm(model_formula, data = x)

}



#' Selection for iForm procedure
#'
#' This is a helper function to run the selection procedure under differen heredity principles
#' and different levels of interactions included in the selection.
#'
#'
#' @param x data.frame of your data with all p predictors
#' @param y vector of observed responses
#' @param p number of predictors in the dataset
#' @param n number of observations in the dataset
#' @param C vector of candidate predictors to consider in this step of the procedure
#' @param S vector of solution predictor selected from previous steps of the procedure
#' @param bic vector of bic values calculated for each step of the procedure
#' @param heredity  a string specifying the heredity to be considered. NULL, weak, strong
#' @param higher_order  logical TRUE indicating to include order-3 interactions in the search (default FALSE)
#' @return the response vector, the solution set of predictors and the calculated bic values
#' @author Kirk Gosik
#' @details
#' Runs the iFORM selection procedure for specifed heredity and level of interactions.
#' It returns the solution to be fit from \code{iForm}
#' @export



iformselect <- function( x, y, p, n, C, S, bic, heredity, higher_order ) {

  repeat{

    RSS <- rss_map_func(C = C, S = S, y = y, data = x)

    S <- c(S, C[which.min(unlist(RSS))])
    C <- C[-which.min(unlist(RSS))]

    order2 <- switch( heredity,
                      `none` = NULL,
                      `strong` = strong_order2(S = S, data = x),
                      `weak` = weak_order2(S = S, C = C, data = x)
    )

    C <- union(C, order2)

    if( higher_order ) {

      order3 <- switch( heredity,
                        `strong` = strong_order3( S = S, data = x ),
                        `weak` = weak_order3(S = S, C = C, data = x)
      )

      C <- union(C, order3)

    }

    bic_val <- log(min(unlist(RSS))/n) + length(S) * (log(n) + 2 * log(p))/n
    bic <- append(bic, bic_val)
    if(length(bic) > n/log(n)) break

  }

  list(y = y, S = S, bic = bic)
}




#' Finding minimum RSS
#'
#' Extended variable selection approaches to jointly model main and interaction effects from high-dimensional data orignally proposed by Hao and Zhang (2014) and extended by Gosik and Wu (2016).
#' Based on a greedy forward approach, their model can identify all possible interaction effects through two algorithms, iFORT and iFORM, which have been proved to possess sure screening property in an ultrahigh-dimensional setting.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param data data.frame of your data with the response and all p predictors
#' @param heredity  a string specifying the heredity to be considered. NULL, weak, strong
#' @param higher_order  logical TRUE indicating to include order-3 interactions in the search (default FALSE)
#' @return a summary of the linear model returned after the selection procedure
#' @author Kirk Gosik
#' @details
#' Runs the iFORM selection procedure on the dataset and returns a linear model
#' of the final selected model.  The model is of an R object of class "lm"
#' @seealso \code{lm}
#' @seealso \code{model.frame}
#' @export
#' @importFrom stats lm
#' @importFrom stats model.frame


rss_map_func <- function( C, S, y, data ) {

  sapply(C, function(candidates) {
    var_names <- c(S, candidates)

    X <- model.matrix(as.formula(paste("~0+", paste(var_names, collapse = "+"))), data = data)

    tryCatch({
      sum((y - X %*% (solve(t(X) %*% X)) %*% (t(X) %*% y)) ^ 2)
    }, error = function(e) Inf)


  })

}


## Heredity Selection

# strong order 2
strong_order2 <- function(S, data) {

  tryCatch({

    main_effects <- sort(S[S %in% names(data)])
    combn(main_effects, 2, paste0, collapse = ":")

  }, error = function(e) NULL)

}


# weak order 2
weak_order2 <- function( S, C, data ) {

  tryCatch({

    main_effects <- sort(S[S %in% names(data)])
    as.vector(outer(main_effects, C, paste, sep = ":"))

  }, error = function(e) NULL)

}


# strong order 3
strong_order3 <- function( S, data ) {

  tryCatch({

    main_effects <- sort(S[S %in% names(data)])
    combn(main_effects, 3, paste0, collapse = ":")

  }, error = function(e) NULL)

}


# weak order 3
weak_order3 <- function( S, C, data ) {

  tryCatch({

    interaction_effects <- unlist(
      Map(function(int_term) paste0(int_term, collapse = ":"),
          Filter(function(vec) {length(vec) == 2}, strsplit(S, "[.]|[:]"))
      )
    )

    weak_three <- as.vector(
      outer(interaction_effects, C, paste, sep = ":")
    )

    as.vector(
      unlist(
        Map(paste0, collapse = ":",
            Map(sort, strsplit(weak_three, ":"))
        )
      )
    )

  }, error = function(e) NULL)
}




