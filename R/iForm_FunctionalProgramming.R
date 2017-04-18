library(foreach)
library(doParallel)
library(iterators)
library(Rcpp)

iForm_funcprog <- function(formula, data, heredity = "strong", higher_order = FALSE) {

  dat <- model.frame(formula, data)
  y <- dat[ , 1]
  x <- dat[ , -1]
  p <- ncol(x)
  n <- nrow(x)
  C <- names(x)
  S <- NULL
  M <- NULL
  bic <- NULL

  fit <- iformselect_funcprog(x, y, p, n, C, S, bic, heredity, higher_order)

  y <- fit$y
  S <- fit$S
  bic <- fit$bic

  model_formula <- as.formula(paste("y ~0+", paste(S[1:which.min(bic)], collapse = "+")))
  lm(model_formula, data = x)

}






iformselect_funcprog <- function( x, y, p, n, C, S, bic, heredity, higher_order ) {

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
    if(length(bic) > 20) break

  }

  list(y = y, S = S, bic = bic)
  stopCluster(cl); gc(reset = TRUE)
}




iformselect_iter <- function( x, y, p, n, C, S, bic, heredity, higher_order ) {

    RSS <- rss_map_func(C = C, S = S, y = y, data = x)

    S <- c(S, C[Position(min, RSS)])
    C <- C[-Position(min, RSS)]

    order2 <- switch( heredity,
                      `none` = NULL,
                      `strong` = strong_order2(S = S, data = x),
                      `weak` = weak_order2(S = S, C = C, data = x)
    )

    S <- union(S, order2)

    if( higher_order ) {

      order3 <- switch( heredity,
                        `strong` = strong_order3( S = S, data = x ),
                        `weak` = weak_order3(S = S, C = C, data = x)
      )

      S <- union(S, order3)

    }

    bic_val <- log(Reduce(min, RSS)/n) + length(S) * (log(n) + 2 * log(p))/n
    bic <- append(bic, bic_val)

  list(x = x,
       y = y,
       p = p,
       n = n,
       C = C,
       S = S,
       bic = bic,
       heredity = heredity,
       higher_order = higher_order)
}



## Help Functions ############################
{
# ## Iterative function application:
# Funcall <- function(f, ...) f(...)
#
# Iterate <- function(f, n = 1)
#   function(x) Reduce(Funcall, rep.int(list(f), n), x, right = TRUE)
#
# Iterate(iformselect_iter(x = x,
#                              y = y,
#                              p = p,
#                              n = n,
#                              C = C,
#                              S = S,
#                              bic = bic,
#                              heredity = heredity,
#                              higher_order = higher_order), n = 5)
#
# Iterate(iformselect_iter, 5)
}




rss_map_func <- function( C, S, y, data ) {

  sapply(C, function(candidates) {
    var_names <- c(S, candidates)

    X <- model.matrix(as.formula(paste("~0+", paste(var_names, collapse = "+"))), data = data)

    tryCatch({
      sum((y - X %*% (solve(t(X) %*% X)) %*% (t(X) %*% y)) ^ 2)
    }, error = function(e) Inf)


  })

}

## fastRSS in Rcpp

# adapted for just RSS calculation
cppFunction(depends = 'RcppArmadillo',
            'List fastRSS(NumericVector yr, NumericMatrix Xr) {
            int n = Xr.nrow(), k = Xr.ncol();
            arma::mat X(Xr.begin(), n, k, false);
            arma::colvec y(yr.begin(), yr.size(), false);
            arma::colvec coef = arma::solve(X, y);
            arma::colvec resid = y - X*coef;
            double rss = pow(arma::norm(resid), 2);

            return List::create(Named("RSS") = rss);
            }')

rss_map_func_cpp <- function( C, S, y, data ) {

  sapply(C, function(candidates) {
    var_names <- c(S, candidates)

    X <- model.matrix(as.formula(paste("~0+", paste(var_names, collapse = "+"))), data = data)

    tryCatch({
      fastRSS(y, X)
    }, error = function(e) Inf)


  })

}


rss_map_func_parallel <- function( C, S, y, data, no_cores = 2 ) {

  cl <- makeCluster(no_cores)
  registerDoParallel(cl)

  foreach(candidate = C, .combine = "c") %dopar%{
    var_names <- c(S, candidate)

    X <- model.matrix(as.formula(paste("~0+", paste(var_names, collapse = "+"))), data = data)

    tryCatch({
      sum((y - X %*% (solve(t(X) %*% X)) %*% (t(X) %*% y)) ^ 2)
    }, error = function(e) Inf)


  }

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


