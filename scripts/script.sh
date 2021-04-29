#!/bin/bash

# This file was adapted from work by Keith James (keithj) and Jaime Tovar Corona
# (jmtc). The original source can be found as part of the wtsi-npg/data_handling
# and wtsi-npg/qc projects here:
#
#   https://github.com/wtsi-npg/data_handling
#   https://github.com/wtsi-npg/npg_qc


set -e -x

unset PERL5LIB

cpanm --notest --installdeps . 
perl Build.PL
./Build

./Build test --verbose

