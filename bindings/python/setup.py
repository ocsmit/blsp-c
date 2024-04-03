#!/usr/bin/env python3

from distutils.core import Extension, setup

package = "blsp"
version = "0.1"

module = Extension(
    package,
    sources=["blsp_py.c"],
    libraries=["blsp", "gsl"],
    include_dirs=["/usr/local/include", "/opt/homebrew/include"],
    library_dirs=["/usr/local/lib"],
)

setup(
    name=package,
    version=version,
    ext_modules=[module],
)
