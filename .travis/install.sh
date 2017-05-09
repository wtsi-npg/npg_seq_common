#!/bin/bash

# This file was adapted from work by Keith James (keithj) and Jaime Tovar Corona
# (jmtc). The original source can be found as part of the wtsi-npg/data_handling
# and wtsi-npg/qc projects here:
#
#   https://github.com/wtsi-npg/data_handling
#   https://github.com/wtsi-npg/npg_qc


set -e -x

sudo apt-get install libgd2-xpm-dev # For npg_tracking
sudo apt-get install liblzma-dev # For npg_qc

### Install third party tools ###

pushd /tmp

# bwa

git clone --branch ${BWA_VERSION} --depth 1 https://github.com/wtsi-npg/bwa.git bwa
pushd bwa
make
ln -s /tmp/bwa/bwa /tmp/bin/bwa
popd

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


# samtools 0.1.19

wget http://sourceforge.net/projects/samtools/files/samtools/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2/download -O samtools-${SAMTOOLS_VERSION}.tar.bz2
tar jxf samtools-${SAMTOOLS_VERSION}.tar.bz2
pushd samtools-${SAMTOOLS_VERSION}
make
ln -s /tmp/samtools-${SAMTOOLS_VERSION}/samtools /tmp/bin/samtools
popd

# staden_io_lib

wget http://sourceforge.net/projects/staden/files/io_lib/${STADEN_IO_LIB_VERSION}/io_lib-${STADEN_IO_LIB_VERSION}.tar.gz/download -O io_lib.tar.gz
tar xzf io_lib.tar.gz
pushd io_lib-${STADEN_IO_LIB_VERSION}
./configure --prefix=/tmp
make
make install
popd


# pb_calibration # for calibration_pu

git clone --branch ${PB_CALIBRATION_VERSION} --depth 1 https://github.com/wtsi-npg/pb_calibration.git
pushd pb_calibration/src
autoreconf --force --install
./configure --with-samtools=/tmp/samtools-0.1.19 --with-io_lib=/tmp --prefix=/tmp
make
make install
popd

# htslib/samtools

git clone --branch ${HTSLIB_VERSION} --depth 1 https://github.com/wtsi-npg/htslib.git htslib
pushd htslib
autoreconf -fi
./configure --prefix=/tmp --enable-plugins
make
make install
popd


git clone --branch ${SAMTOOLS1_VERSION} --depth 1 https://github.com/wtsi-npg/samtools.git samtools-irods
pushd samtools-irods
mkdir -p acinclude.m4
pushd acinclude.m4
curl -L https://github.com/samtools/samtools/files/62424/ax_with_htslib.m4.txt > ax_with_htslib.m4
curl -L 'http://git.savannah.gnu.org/gitweb/?p=autoconf-archive.git;a=blob_plain;f=m4/ax_with_curses.m4;hb=0351b066631215b4fdc3c672a8ef90b233687655' > ax_with_curses.m4
popd
aclocal -I acinclude.m4
autoreconf -i
./configure --prefix=/tmp --with-htslib=/tmp/htslib --enable-plugins
make
ln -s /tmp/samtools-irods/samtools /tmp/bin/samtools_irods
popd


# picard
wget https://sourceforge.net/projects/picard/files/picard-tools/${PICARD_VERSION}/picard-tools-${PICARD_VERSION}.zip/download -O picard-tools-${PICARD_VERSION}.zip
unzip picard-tools-${PICARD_VERSION}.zip

# biobambam

wget https://github.com/gt1/biobambam2/releases/download/${BIOBAMBAM_VERSION}/biobambam2-${BIOBAMBAM_VERSION}-x86_64-etch-linux-gnu.tar.gz -O biobambam2.tar.gz
mkdir biobambam2
tar xzf biobambam2.tar.gz -C biobambam2 --strip-components 1


popd

# Third party tools install done

# CPAN as in npg_npg_deploy
cpanm --notest --reinstall App::cpanminus
cpanm --no-lwp --notest https://github.com/wtsi-npg/perl-dnap-utilities/releases/download/${DNAP_UTILITIES_VERSION}/WTSI-DNAP-Utilities-${DNAP_UTILITIES_VERSION}.tar.gz

# WTSI NPG Perl repo dependencies
cd /tmp
git clone --branch devel --depth 1 https://github.com/wtsi-npg/ml_warehouse.git ml_warehouse.git
git clone --branch devel --depth 1 https://github.com/wtsi-npg/npg_tracking.git npg_tracking.git
git clone --branch devel --depth 1 https://github.com/wtsi-npg/npg_qc.git npg_qc.git


repos="/tmp/ml_warehouse.git /tmp/npg_tracking.git /tmp/npg_qc.git"

for repo in $repos
do
      cd "$repo"
        cpanm --quiet --notest --installdeps . || find /home/travis/.cpanm/work -cmin -1 -name '*.log' -exec tail -n20  {} \;
          perl Build.PL
            ./Build
              ./Build install
          done

          cd "$TRAVIS_BUILD_DIR"

