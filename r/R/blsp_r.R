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
                      theta_mean = c(0.05, 1, 120, 8, 290, 8, 0.001),
                      theta_sd = rep(0.1, 7),
                      iterations = 7000,
                      burn = 2000,
                      store.samples = TRUE) {
    ## condition y_idx data for internal time series struct
    yidx <- c(0, lengths(data))
    for (i in 2:length(yidx)) yidx[i] <- yidx[i] + yidx[i - 1]

    ## FFI for MCMC sampler
    samples <- .Call(
        C_run_blsp, unlist(data), unlist(doy), yidx, theta_mean, theta_sd,
        as.integer(iterations), as.integer(burn)
    )

    ## output matrix
    sample_mat <- matrix(samples, (length(yidx) - 1) * iterations, 9,
        byrow = TRUE
    )

    mydt <- data.table::as.data.table(sample_mat)
    data.table::setnames(mydt, names(mydt), c(
        "year", "iter", "m1", "m2", "m3", "m4", "m5", "m6", "m6"
    ))

    mydt <- mydt[mydt$iter >= (burn - 1), ] # iter count in C starts at 0

    dt <- list()
    for (i in 0:(length(yidx) - 2)) {
        year_dt <- mydt[mydt$year == i, ]
        theta3 <- quantile(year_dt[, "m3"], c(0.025, 0.5, 0.975))
        theta5 <- quantile(year_dt[, "m5"], c(0.025, 0.5, 0.975))
        dt[[i + 1]] <- list(
            theta3[1], theta3[2], theta3[3],
            theta5[1], theta5[2], theta5[3]
        )
    }

    dtp <- data.table::rbindlist(dt)
    colnames(dtp) <- c(
        "midgup_lower", "midgup", "midgup_upper",
        "midgdown_lower", "midgdown", "midgdown_upper"
    )


    ## construct result
    samples_ret <- NA
    if (store.samples) samples_ret <- sample_mat
    result <- list(
        nyrs = (length(yidx) - 1),
        iterations = iterations,
        burn = burn,
        pheno = dtp,
        samples = samples_ret
    )

    class(result) <- "BLSPModel"

    ## TODO: return quantiles instead of every sample
    return(result)
}
