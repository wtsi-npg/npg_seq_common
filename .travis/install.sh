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
# for compiling tools in npg_qc since they expect to find samtools headers in /include
# relative to which samtools in PATH
ln -s /tmp/samtools /tmp/samtools/lib
ln -s /tmp/samtools /tmp/samtools/include
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

# salmon pre-0.10.x release
wget https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON0_8_VERSION}/Salmon-${SALMON0_8_VERSION}_linux_x86_64.tar.gz
tar xzf Salmon-${SALMON0_8_VERSION}_linux_x86_64.tar.gz
mkdir -p /tmp/salmon/${SALMON0_8_VERSION}
cp -r /tmp/Salmon-${SALMON0_8_VERSION}_linux_x86_64/* /tmp/salmon/${SALMON0_8_VERSION}
mv /tmp/salmon/${SALMON0_8_VERSION}/bin/salmon /tmp/salmon/${SALMON0_8_VERSION}/bin/salmon0_8


# salmon post-0.10.x release
wget https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON0_10_VERSION}/salmon-${SALMON0_10_VERSION}-linux_x86_64.tar.gz
tar xzf salmon-${SALMON0_10_VERSION}-linux_x86_64.tar.gz
mkdir -p /tmp/salmon/${SALMON0_10_VERSION}
cp -r /tmp/salmon-${SALMON0_10_VERSION}-linux_x86_64/* /tmp/salmon/${SALMON0_10_VERSION}
mv /tmp/salmon/${SALMON0_10_VERSION}/bin/salmon /tmp/salmon/${SALMON0_10_VERSION}/bin/salmon0_10

# tophat2
wget https://ccb.jhu.edu/software/tophat/downloads/tophat-${TOPHAT2_VERSION}.Linux_x86_64.tar.gz
tar xzf tophat-${TOPHAT2_VERSION}.Linux_x86_64.tar.gz
cp -r /tmp/tophat-${TOPHAT2_VERSION}.Linux_x86_64/* /tmp/bin/

# gffread
wget http://ccb.jhu.edu/software/stringtie/dl/gffread-${GFFREAD_VERSION}.Linux_x86_64.tar.gz
tar xzf gffread-${GFFREAD_VERSION}.Linux_x86_64.tar.gz
ln -s /tmp/gffread-${GFFREAD_VERSION}.Linux_x86_64/gffread /tmp/bin/gffread

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

