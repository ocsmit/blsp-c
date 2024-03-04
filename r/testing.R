data_path <- file.path("test_data.Rdata")
load(data_path)

devtools::clean_dll()
devtools::load_all(reset = T)

theta_mean <- c(
    -1.998392, 0.960355, 120.702350, 9.263498,
    288.853856, 9.166685, -6.592421
)
theta_sd <- c(
    0.07057906, 0.05609551, 1.08944966, 0.88183154,
    1.55979462, 1.20727157, 0.19881890
)

tt <- blspR::blsp_test(y, t, theta_mean, theta_sd)
print(tt$pheno)
