#! /bin/sh
# aclocal \
# && automake --add-missing \
# && autoconf
autoreconf --install || exit 1
