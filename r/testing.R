data_path <- file.path("../test_data.Rdata")
load(data_path)

devtools::clean_dll()
devtools::load_all(reset=T)


theta_mean <- c(
    -1.998392, 0.960355, 120.702350, 9.263498,
    288.853856, 9.166685, -6.592421
)
theta_sd <- c(
    0.07057906, 0.05609551, 1.08944966, 0.88183154,
    1.55979462, 1.20727157, 0.19881890
)


tt <- blspR::blsp_test(y, t, theta_mean, theta_sd)

mydt <- data.table::as.data.table(tt$samples)
dt <- list()
for (i in 0:19) {
    year_dt = mydt[V2 == i]
    theta3 <- quantile(year_dt[, V5], c(0.025, 0.5, 0.975))
    theta5 <- quantile(year_dt[, V7], c(0.025, 0.5, 0.975))
    dt[[i+1]] <- list(theta3[1], theta3[2], theta3[3],
        theta5[1], theta5[2], theta5[3])
}

dtp <- data.table::rbindlist(dt)
colnames(dtp) <- c("midgup_lower", "midgup", "midgup_upper",
    "midgdown_lower", "midgdown", "midgdown_upper")

dtp
