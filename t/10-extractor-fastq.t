use strict;
use warnings;
use Test::More tests => 28;
use Test::Exception;
use File::Temp qw/tempdir/;
use Cwd qw/cwd/;
use File::Compare;

my $source2 = q[t/data/s2.fastq];
my $dir = tempdir(CLEANUP => 1);

my $empty = $dir . q[/empty_fastq.fastq];
`touch $empty`;

my @subs = qw(read_count split_reads );
use_ok( 'npg_common::extractor::fastq', @subs);

foreach my $sub (@subs) {
   can_ok(__PACKAGE__, $sub);
}

{
   is(read_count($source2), '26', 'read count when fastqcheck file is not present');  
}

{
   throws_ok {split_reads()} qr/Input file name should be given/, 'error when no args are given';
   throws_ok {split_reads($empty)} qr/Read lengths for the target files should be given as an array reference/, 'error when no target read length is given';
   throws_ok {split_reads($empty, [0])} qr/Target read length should be positive/, 'error when zero target read length is given';
   throws_ok {split_reads($empty, [55, 0])} qr/Target read length should be positive/,'error when zero target read length is given';
   throws_ok {split_reads($empty, [55, -55])} qr/Target read length should be positive/,'error when negative target read length is given';  
}

{
   my $current = cwd;
   chdir $dir;
   lives_ok {split_reads($empty, [37, 37])} 'splitting an empty file leaves';
   ok(-e q[37_1.fastq] && -e q[37_2.fastq], 'output for splitting an empty file exists');
   chdir $current;

   my $out1 = $dir . q[/dodo1];
   my $out2 = $dir . q[/dodo2];
   lives_ok {split_reads($empty, [37, 37], [$out1, $out2])} 'splitting an empty file leaves';
   ok(-e $out1 && -e $out2, 'output for splitting an empty file exists');
}

{
   my $out1 = $dir . q[/dodo3];
   my $out2 = $dir . q[/dodo4];
   lives_ok {split_reads($empty, [], [$out1, $out2])} 'splitting an empty file in halves leaves';
   ok(-e $out1 && -e $out2, 'output for splitting an empty file exists');
}

{
   my $current = cwd;
   chdir $dir;
   throws_ok {split_reads($current . q[/t/data/s2.fastq], [37])} qr/is too short/, 'error when the read is too short';
   chdir $current;
}

{  
   my $current = cwd;
   chdir $dir;
   lives_ok {split_reads($current . q[/t/data/ss1.fastq], [3])} 'splitting reads lives';
   ok(-e q[3_1.fastq], 'output exists');
   ok(!-e q[3_2.fastq], 'second output does not exist');
   chdir $current;
   ok(compare(q[t/data/ss1_split3.fastq], $dir . q[/3_1.fastq])==0, 'correct output for a single file'); 
}

{  
   my $current = cwd;
   chdir $dir;
   lives_ok {split_reads($current . q[/t/data/1008_s_1.fastq], [37,37])} 'splitting reads in two with explicit lengths lives';
   ok(-e q[37_1.fastq] && -e q[37_2.fastq], 'output exists');
   chdir $current;
   ok(compare(q[t/data/1008_1_1.fastq], $dir . q[/37_1.fastq])==0, 'correct output for a forward file');
   ok(compare(q[t/data/1008_1_2.fastq], $dir . q[/37_2.fastq])==0, 'correct output for a reverse file'); 
}

{  
   my $current = cwd;
   chdir $dir;
   my $f1 = q[half1.fastq];
   my $f2 = q[half2.fastq];
   lives_ok {split_reads($current . q[/t/data/1008_s_1.fastq], [], [$f1,$f2])} 'splitting reads in halves';
   ok(-e $f1 && -e $f2, 'output exists');
   chdir $current;
   ok(compare(q[t/data/1008_1_1.fastq], $dir . q[/] . $f1)==0, 'correct output for a forward file');
   ok(compare(q[t/data/1008_1_2.fastq], $dir . q[/] . $f2)==0, 'correct output for a reverse file'); 
}

1;
