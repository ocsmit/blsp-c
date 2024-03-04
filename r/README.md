# BLSP-C R bindings

BLSP-C will need to be installed on your system, if it is not please see [here](../README.md).
The compiler on macOS might not have `/usr/local` in its path for the linker to search.
If this is the case add the following to your `~/.R/Makevars` file (create it if it does not exist)
```sh
CPPFLAGS += -I/usr/local/include
LDFLAGS += -L/usr/local/lib
```

## R Devtools
[devtools](https://devtools.r-lib.org/) makes it easy to install directly from github or from source. To install from github use the following
```R
# from github
devtools::install_github("ocsmit/blsp-c/r")
```
otherwise clone the repo and run
```R
devtools::clean_dll()
devtools::load_all(reset = T)
```

## R CMD
To install using the standard R tooling run the following from your shell
```sh
R CMD INSTALL .
```
