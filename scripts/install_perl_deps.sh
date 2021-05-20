#!/bin/bash

set -e -x

WTSI_NPG_GITHUB_URL=$1
WTSI_NPG_BUILD_BRANCH=$2

eval $(perl -I ~/perl5ext/lib/perl5/ -Mlocal::lib=~/perl5ext)
cpanm --quiet --notest LWP::Protocol::https
cpanm --quiet --notest Alien::Tidyp # For npg_tracking
cpanm --quiet --notest Module::Build
cpanm --quiet --notest Net::SSLeay
cpanm --quiet --notest https://github.com/chapmanb/vcftools-cpan/archive/v0.953.tar.gz # for npg_qc

# Git branch to merge to or custom branch
# WTSI NPG Perl repo dependencies

repos="" 
for repo in perl-dnap-utilities ml_warehouse npg_tracking; do
  cd /tmp
  # Always clone master when using depth 1 to get current tag
  git clone --branch master --depth 1 ${WTSI_NPG_GITHUB_URL}/${repo}.git ${repo}.git
  cd /tmp/${repo}.git
  # Shift off master to appropriate branch (if possible)
  git ls-remote --heads --exit-code origin ${WTSI_NPG_BUILD_BRANCH} && git pull origin ${WTSI_NPG_BUILD_BRANCH} && echo "Switched to branch ${WTSI_NPG_BUILD_BRANCH}"
  repos=$repos" /tmp/${repo}.git"
done

for repo in $repos
do
    export PERL5LIB=$repo/blib/lib:$PERL5LIB:$repo/lib
done

for repo in $repos
do
    cd $repo
    cpanm  --quiet --notest --installdeps .
    perl Build.PL
    ./Build
done

eval $(perl -I ~/perl5ext/lib/perl5/ -Mlocal::lib=~/perl5npg)
for repo in $repos
do
    cd $repo
    cpanm  --quiet --notest --installdeps .
    ./Build install
done
cd
