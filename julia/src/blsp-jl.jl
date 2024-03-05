module blsp-jl


# This will currently segfault because no inputs
@ccall "/usr/local/lib/libblsp".BLSP_sampler()::Int32


end
