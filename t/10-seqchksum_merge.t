use strict;
use warnings;
use Test::More tests => 22;
use Cwd qw/abs_path getcwd/;
use File::Temp qw/tempdir/;
use File::Slurp;
use Digest::MD5;
use JSON;

use npg_tracking::util::build qw/git_tag/;

# test the seqchksum_merge script

my $tmpdir = tempdir('seqchksum_merge_test_XXXXXX', CLEANUP => 1, DIR => '/tmp' );
print "Created temporary directory: ".abs_path($tmpdir)."\n";
my $startDir = getcwd();
my $seqchksum_data = abs_path('t/data/seqchksum_merge');
unless (-d $seqchksum_data) {
    die "Cannot find data directory $seqchksum_data\n";
}

is(system("$startDir/bin/seqchksum_merge.pl"), 0, 'seqchksum_merge.pl exit status');

# first, a normal successful run
is(system("$startDir/bin/seqchksum_merge.pl $seqchksum_data/target.seqchksum $seqchksum_data/phix.seqchksum > $tmpdir/merged.matching.seqchksum"), 0, 'seqchksum_merge.pl create matching test output exit status');
is(system("cmp -s $tmpdir/merged.matching.seqchksum $seqchksum_data/merged.seqchksum"), 0, 'seqchksum_merge.pl matching cmp exit status');

# merge sam, bam, cram and seqchksum files - N.B. sort order isn't guaranteed so sort results before comparing
is(system("export LC_ALL=C; $startDir/bin/seqchksum_merge.pl $seqchksum_data/format1.sam $seqchksum_data/format2.bam $seqchksum_data/format3.cram $seqchksum_data/format4.seqchksum | sort > $tmpdir/merged.formats.matching.seqchksum"), 0, 'seqchksum_merge.pl create different formats test output exit status');
is(system("cmp -s $tmpdir/merged.formats.matching.seqchksum $seqchksum_data/merged.formats.seqchksum"), 0, 'seqchksum_merge.pl different formats cmp exit status');

# the files used for merge contain different read groups
is(system("$startDir/bin/seqchksum_merge.pl $seqchksum_data/group1.seqchksum $seqchksum_data/group2.seqchksum > $tmpdir/merged.different.seqchksum"), 0, 'seqchksum_merge.pl create matching test output exit status');
is(system("cmp -s $tmpdir/merged.different.seqchksum $seqchksum_data/merged.different.seqchksum"), 0, 'seqchksum_merge.pl different read groups exit status');

# the files used for the merge do not "partition" the reads in the combined file (extra occurrence of target.seqchksum), so comparison fails
is(system("$startDir/bin/seqchksum_merge.pl $seqchksum_data/target.seqchksum $seqchksum_data/phix.seqchksum $seqchksum_data/target.seqchksum > $tmpdir/merged.extra_input.seqchksum"), 0, 'seqchksum_merge.pl create test output 2 exit status');
isnt(system("cmp -s $tmpdir/merged.extra_input.seqchksum $seqchksum_data/merged.seqchksum"), 0, 'seqchksum_merge.pl extra component file cmp exit status');

# the files used for the merge do not "partition" the reads in the combined file (missing target.seqchksum), so comparison fails
is(system("$startDir/bin/seqchksum_merge.pl $seqchksum_data/phix.seqchksum > $tmpdir/merged.missing_input.seqchksum"), 0, 'seqchksum_merge.pl create test output 2 exit status');
isnt(system("cmp -s $tmpdir/merged.missing_input.seqchksum $seqchksum_data/merged.seqchksum"), 0, 'seqchksum_merge.pl missing component file cmp exit status');

# errors because format of input files is not consistent
isnt(system("$startDir/bin/seqchksum_merge.pl $seqchksum_data/target.missing_col.seqchksum $seqchksum_data/phix.seqchksum > /dev/null"), 0, 'seqchksum_merge.pl input missing col status');
isnt(system("$startDir/bin/seqchksum_merge.pl $seqchksum_data/target.extra_col.seqchksum $seqchksum_data/phix.seqchksum > /dev/null"), 0, 'seqchksum_merge.pl input extra col status');
# column 4 flagged (by default) as constant does not match across input files
isnt(system("$startDir/bin/seqchksum_merge.pl $seqchksum_data/target.seqchksum $seqchksum_data/phix.inconstant.seqchksum > /dev/null"), 0, 'seqchksum_merge.pl inconstant constant exit status');

# Initial input file contains only comments. No useful output, but not an error.
is(system("$startDir/bin/seqchksum_merge.pl $seqchksum_data/target.all_comments.seqchksum $seqchksum_data/phix.seqchksum > /dev/null"), 0, 'seqchksum_merge.pl with all-comments input exit status');

# incorrectly identify column 6 as both a count (accumulate) and a partition field. Expect fatal error from script.
isnt(system("$startDir/bin/seqchksum_merge.pl -a5 -p5 $seqchksum_data/target.seqchksum $seqchksum_data/phix.seqchksum > /dev/null"), 0, 'seqchksum_merge.pl chksum, conflicting column type definition exit status');

# incorrectly identify column 5 (containing check sum values) as a count (accumulate) field. Merged file successfully generated, comparison fails
is(system("$startDir/bin/seqchksum_merge.pl -a5 $seqchksum_data/target.seqchksum $seqchksum_data/phix.seqchksum > $tmpdir/merged.chksum_as_acc.seqchksum"), 0, 'seqchksum_merge.pl chksum as accumulate test output exit status');
isnt(system("cmp -s $tmpdir/merged.chksum_as_acc.seqchksum $seqchksum_data/merged.seqchksum"), 0, 'seqchksum_merge.pl chksum fld as acc cmp exit status');

# incorrectly identify column 3 (containing count values) as a chksum field. Merged file successfully generated, comparison fails
is(system("$startDir/bin/seqchksum_merge.pl -c3 $seqchksum_data/target.seqchksum $seqchksum_data/phix.seqchksum > $tmpdir/merged.acc_as_chksum.seqchksum"), 0, 'seqchksum_merge.pl accumulate as chksum test output exit status');
isnt(system("cmp -s $tmpdir/merged.acc_as_chksum.seqchksum $seqchksum_data/merged.seqchksum"), 0, 'seqchksum_merge.pl acc fld as chksum cmp exit status');

# incorrectly identify column 6 (containing check sum values) as a partition field. Merged file successfully generated, comparison fails
is(system("$startDir/bin/seqchksum_merge.pl -p6 $seqchksum_data/target.seqchksum $seqchksum_data/phix.seqchksum > $tmpdir/merged.chksum_as_acc.seqchksum"), 0, 'seqchksum_merge.pl chksum as partition test output exit status');
isnt(system("cmp -s $tmpdir/merged.chksum_as_partition.seqchksum $seqchksum_data/merged.seqchksum"), 0, 'seqchksum_merge.pl chksum fld as partition cmp exit status');

1;
