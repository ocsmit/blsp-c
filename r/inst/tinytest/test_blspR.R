data(SimulatedPhenology)

tt <- with(
    SimulatedPhenology,
    blspR::blsp_test(y, t, theta_mean, theta_sd, iterations = 50000)
)
tinytest::expect_equivalent(tt$pheno, SimulatedPhenology$expected_result$pheno)

names(tt)
class(tt)

head(tt$samples)
