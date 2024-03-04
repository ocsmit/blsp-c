# MCMC sampler for the Bayesian Land Surface Phenology model.

important: *This is under active development and subject to change at any point.*

This is a monorepo containing the `C` source code for a Metroplois Hastings + Gibbs MCMC sampler, and interfaces for several higher level languages.

For more info please see the paper
> Gao, X., Gray, J. M., & Reich, B. J. (2021). Long-term, medium spatial resolution annual land surface phenology with a Bayesian hierarchical model. Remote Sensing of Environment, 261, 112484. https://doi.org/10.1016/j.rse.2021.112484

and the [original algorithms git repo](https://github.com/ncsuSEAL/Bayesian_LSP).


Depends on [GNU Scientific Library](https://www.gnu.org/software/gsl/)


# Building for Linux & macOS
## Dependencies
### Autotools
Ensure that you have autotools installed.
They are required to generate the build scripts and compile the software.

On macOS this can be accomplished with homebrew.
```sh
brew install autoconf automake libtool
```
Linux systems can use their respective package managers, e.g.
```sh
# ubuntu
apt install autoconf libtool automake
```

### GNU Scientific Library

On mac
```sh
brew install gsl
```
Linux systems can use their respective package managers, e.g.
```sh
# ubuntu
apt install libgsl-dev
```

### Installing from source

Once dependencies are installed, you can build and install.

Generate the require build scripts with
```sh
./autogen.sh
```

To build and install on your system run the following
```sh
./compile
make install # May need sudo
```
This assumes that GSL is discoverable on your system (installed to a common location).

# Language bindings
  - [Building R package](r/README.md)
  - Python (todo)
  - Julia (todo)
