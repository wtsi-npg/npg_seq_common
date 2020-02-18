#!/bin/bash

# This file was adapted from work by Keith James (keithj) and Jaime Tovar Corona
# (jmtc). The original source can be found as part of the wtsi-npg/data_handling
# and wtsi-npg/qc projects here:
#
#   https://github.com/wtsi-npg/data_handling
#   https://github.com/wtsi-npg/npg_qc


set -e -x

### Install third party tools ###

pushd /tmp

# bwa0_6

git clone --branch ${BWA0_6_VERSION} --depth 1 https://github.com/lh3/bwa.git bwa0_6
pushd bwa0_6
make
ln -s /tmp/bwa0_6/bwa /tmp/bin/bwa0_6
popd

# blat

wget https://users.soe.ucsc.edu/~kent/src/blatSrc${BLAT_VERSION}.zip
unzip -q blatSrc${BLAT_VERSION}
pushd blatSrc
MACHTYPE="`uname -m`"
mkdir -p $HOME/bin/$MACHTYPE
mkdir -p $HOME/lib/$MACHTYPE
make
popd

# bowtie

git clone --branch ${BOWTIE_VERSION} --depth 1 https://github.com/dkj/bowtie.git bowtie
pushd bowtie
make
popd

ln -s /tmp/bowtie/bowtie /tmp/bin/bowtie
ln -s /tmp/bowtie/bowtie-build /tmp/bin/bowtie-build
ln -s /tmp/bowtie/bowtie-inspect /tmp/bin/bowtie-inspect


# bowtie2

git clone --branch ${BOWTIE2_VERSION} --depth 1 https://github.com/BenLangmead/bowtie2.git bowtie2
pushd bowtie2
make
popd

ln -s /tmp/bowtie2/bowtie2 /tmp/bin/bowtie2
ln -s /tmp/bowtie2/bowtie2-align-l /tmp/bin/bowtie2-align-l
ln -s /tmp/bowtie2/bowtie2-align-s /tmp/bin/bowtie2-align-s

ln -s /tmp/bowtie2/bowtie2-build /tmp/bin/bowtie2-build
ln -s /tmp/bowtie2/bowtie2-build-l /tmp/bin/bowtie2-build-l
ln -s /tmp/bowtie2/bowtie2-build-s /tmp/bin/bowtie2-build-s

ln -s /tmp/bowtie2/bowtie2-inspect /tmp/bin/bowtie2-inspect
ln -s /tmp/bowtie2/bowtie2-inspect-l /tmp/bin/bowtie2-inspect-l
ln -s /tmp/bowtie2/bowtie2-inspect-s /tmp/bin/bowtie2-inspect-s

# staden_io_lib

wget http://sourceforge.net/projects/staden/files/io_lib/${STADEN_IO_LIB_VERSION}/io_lib-${STADEN_IO_LIB_VERSION}.tar.gz/download -O io_lib.tar.gz
tar xzf io_lib.tar.gz
pushd io_lib-${STADEN_IO_LIB_VERSION}
./configure --prefix=/tmp
make
make install
popd

# htslib/samtools

git clone --branch ${HTSLIB_VERSION} --depth 1 https://github.com/samtools/htslib htslib
pushd htslib
autoreconf -fi
./configure --prefix=/tmp --enable-plugins
make
make install
popd

git clone --branch ${SAMTOOLS_VERSION} --depth 1 https://github.com/samtools/samtools samtools
pushd samtools
mkdir -p acinclude.m4
pushd acinclude.m4
curl -L https://github.com/samtools/samtools/files/62424/ax_with_htslib.m4.txt > ax_with_htslib.m4
curl -L 'http://git.savannah.gnu.org/gitweb/?p=autoconf-archive.git;a=blob_plain;f=m4/ax_with_curses.m4;hb=0351b066631215b4fdc3c672a8ef90b233687655' > ax_with_curses.m4
popd
aclocal -I acinclude.m4
autoreconf -i
./configure --prefix=/tmp --with-htslib=/tmp/htslib --enable-plugins --without-curses
make
# for other tools to find samtools and its alias samtools_irods
ln -s /tmp/samtools/samtools /tmp/bin/samtools_irods
ln -s /tmp/samtools/samtools /tmp/bin/samtools
popd

# picard
wget https://sourceforge.net/projects/picard/files/picard-tools/${PICARD_VERSION}/picard-tools-${PICARD_VERSION}.zip/download -O picard-tools-${PICARD_VERSION}.zip
unzip picard-tools-${PICARD_VERSION}.zip

# biobambam

wget https://github.com/gt1/biobambam2/releases/download/${BIOBAMBAM_VERSION}/biobambam2-${BIOBAMBAM_VERSION}-x86_64-etch-linux-gnu.tar.gz -O biobambam2.tar.gz
mkdir biobambam2
tar xzf biobambam2.tar.gz -C biobambam2 --strip-components 1

# star

wget https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.zip
unzip ${STAR_VERSION}.zip
ln -s /tmp/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/STAR /tmp/bin/star

# minimap2

curl -L https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 | tar -jxvf -
ln -s /tmp/minimap2-${MINIMAP2_VERSION}_x64-linux/minimap2 /tmp/bin/minimap2

# hisat2

wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-${HISAT2_VERSION}-Linux_x86_64.zip
unzip hisat2-${HISAT2_VERSION}-Linux_x86_64.zip
for file in hisat2  hisat2-align-s  hisat2-align-l \
            hisat2-build  hisat2-build-s  hisat2-build-l \
            hisat2-inspect  hisat2-inspect-s  hisat2-inspect-l
do
    ln -s /tmp/hisat2-${HISAT2_VERSION}/${file} /tmp/bin/${file}
done

popd

# Third party tools install done
cpanm --quiet --notest LWP::Protocol::https

# Git branch to merge to or custom branch
WTSI_NPG_BUILD_BRANCH=${WTSI_NPG_BUILD_BRANCH:=$TRAVIS_BRANCH}
# WTSI NPG Perl repo dependencies
repos="perl-dnap-utilities ml_warehouse npg_tracking"

for repo in $repos
do
  cd /tmp
  # Always clone master when using depth 1 to get current tag
  git clone --branch master --depth 1 ${WTSI_NPG_GITHUB_URL}/${repo}.git ${repo}.git
  cd /tmp/${repo}.git
  # Shift off master to appropriate branch (if possible)
  git ls-remote --heads --exit-code origin ${WTSI_NPG_BUILD_BRANCH} && git pull origin ${WTSI_NPG_BUILD_BRANCH} && echo "Switched to branch ${WTSI_NPG_BUILD_BRANCH}"
  cpanm --quiet --notest --installdeps . || find /home/travis/.cpanm/work -cmin -1 -name '*.log' -exec tail -n20  {} \;
  perl Build.PL
  ./Build
  ./Build install

done

cd "$TRAVIS_BUILD_DIR"

