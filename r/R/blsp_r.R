#' Run blsp-c MCMC sampler
#'
#' @param data           data vector
#' @param doy            date of year vector
#' @param theta_mean     vector[7] of initial parameter means for sampler
#' @param theta_sd       vector[7] of initial parameter stddev for sampler
#' @param iterations     number of iterations to sample
#' @param burn           burn in period for mcmc sampler
#' @param store.samples  NOT USED TODO
#'
#' @return BLSPModel object
blsp_test <- function(data,
                      doy,
                      theta_mean,
                      theta_sd,
                      iterations = 7000,
                      burn = 2000,
                      store.samples = FALSE
                      ) {
    ## condition y_idx data for internal time series struct
    yidx <- c(0, lengths(data))
    for (i in 2:length(yidx)) yidx[i] <- yidx[i] + yidx[i - 1]

    ## FFI for MCMC sampler
    samples <- .Call(
        C_run_blsp, unlist(data), unlist(doy), yidx, theta_mean, theta_sd,
        as.integer(iterations), as.integer(burn)
    ))

    ## output matrix
    sample_mat <- matrix(samples, (length(yidx) - 1) * iterations, 9, byrow = TRUE)

    ## construct result
    result <- list(
        nyrs = (length(yidx) - 1),
        iterations = iterations,
        burn = burn,
        samples = sample_mat
    )
    class(result) <- "BLSPModel"

    ## TODO: return quantiles instead of every sample
    return(result)
}
