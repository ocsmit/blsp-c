data_path <- file.path("../test_data.Rdata")

load(data_path)
ls()

devtools::load_all()

blspR::blsp_test(unlist(y), unlist(t), yr_idx)
