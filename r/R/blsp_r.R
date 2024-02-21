blsp_test <- function(data,
                      doy,
                      theta_mean,
                      theta_sd,
                      iterations = 7000,
                      burn = 2000,
                      store.samples = FALSE) {
    yidx <- c(0, lengths(data))
    # condition y_idx data for internal time series struct
    for (i in 2:length(yidx)) yidx[i] <- yidx[i] + yidx[i - 1]

    samples <- .Call(
        C_run_blsp, unlist(data), unlist(doy), yidx, theta_mean, theta_sd,
        as.integer(iterations), as.integer(burn)
    ) # , PACKAGE = "blspR")

    sample_mat <- matrix(samples, (length(yidx) - 1) * iterations, 9, byrow = TRUE)

    result <- list(
        nyrs = (length(yidx) - 1),
        iterations = iterations,
        burn = burn,
        samples = sample_mat
    )

    class(result) <- "BLSPModel"
    return(result)
}
