use strict;
use warnings;
use Test::More tests => 6;
use Test::Exception;
use Test::Deep;
use Cwd;
use File::Temp qw(tempdir);

use npg_qc::autoqc::results::bam_flagstats;
use Test::MockModule;
my $mock = Test::MockModule->new('npg_qc::autoqc::results::bam_flagstats');
$mock->mock(execute  => sub { return 1; });
$mock->mock(set_info => sub { return 1; });
$mock->mock(store    => sub { return 1; });

local $ENV{'http_proxy'}='http://wibble.com';

subtest 'subtest 1' => sub {
  use_ok('npg_common::sequence::BAM_MarkDuplicate');
};

subtest 'subtest 2' => sub {
  my $num_tests = 6;
  plan tests => $num_tests;

  SKIP: {
    skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
          $num_tests unless ($ENV{'TOOLS_INSTALLED'});
    my $temp_dir = tempdir( CLEANUP => 1);
    my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
                 input_bam         => 'input.bam',
                 output_bam        => 'output.bam',
                 metrics_json_dir  => 'qc',
                 id_run            => 35,
                 position          => 1);
    isa_ok($bam, 'npg_common::sequence::BAM_MarkDuplicate');

    $bam = npg_common::sequence::BAM_MarkDuplicate->new(
                 input_bam        => 'input.bam',
                 output_bam       => 'output.bam',
                 metrics_json_dir => 'qc',
                 no_alignment     => 0,
                 id_run           => 35,
                 position         => 1
               );
    lives_ok {$bam->temp_dir} 'temp dir generated';
    lives_ok {$bam->metrics_file()} 'temp metrics file';
    $bam->metrics_file('metrics.txt');
    $bam->temp_dir($temp_dir);
    like($bam->mark_duplicate_cmd(), qr/bammarkduplicates2 I=$temp_dir\/sorted.bam O=\/dev\/stdout tmpfile=$temp_dir\/ M=metrics\.txt/,
      'correct picard command with absolute path to jar');
    like($bam->bamseqchksum_cmd(q{bam}), qr{\Qbamseqchksum verbose=0 inputformat=bam\E},
      'correct bamseqchksum command for a bam file');
    like($bam->bamseqchksum_cmd(q{cram}), qr{\Qbamseqchksum verbose=0 inputformat=cram\E},
      'correct bamseqchksum command for a cram file with no reference');
  };
};

my $bammarkduplicates = `which bammarkduplicates2`;
$bammarkduplicates = `readlink -f $bammarkduplicates`;
chomp $bammarkduplicates;
my $bamseqchksum = `which bamseqchksum`;
$bamseqchksum = `readlink -f $bamseqchksum`;
chomp $bamseqchksum;
my $scramble = `which scramble`;
$scramble = `readlink -f $scramble`;
chomp $scramble;
my $cram_index = `which cram_index`;
$cram_index = `readlink -f $cram_index`;
chomp $cram_index;

subtest 'subtest 3' => sub {
  my $num_tests = 31;
  plan tests => $num_tests;

  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
         $num_tests unless ($ENV{'TOOLS_INSTALLED'});

      my $temp_dir = tempdir( CLEANUP => 1);
      my $input = join q[/], $temp_dir, '4392_1.bam';
      my $output_root     = join q[/], $temp_dir, 'output_mk';
      my $output_bam      = $output_root . q[.bam];

      my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
                 input_bam        => $input,
                 output_bam       => $output_bam,
                 metrics_json_dir => $temp_dir,
                 temp_dir         => $temp_dir,
                 reference        => 't/data/references/Plasmodium_falciparum/default/all/bwa0_6/Pf3D7_v3.fasta',
                 id_run           => 35,
                 position         => 1
               );
      my $md_metrics_file = $bam->metrics_file;
      is($md_metrics_file, $output_root.q[.markdups_metrics.txt], 'metrics file path generated relative to mk root');
      my $expected_mark_duplicate_cmd =
        qq{bammarkduplicates2 I=$input O=/dev/stdout tmpfile=$temp_dir/ M=$md_metrics_file};
      like($bam->mark_duplicate_cmd(), qr/$expected_mark_duplicate_cmd/, 'correct biobambam command');
      ok($bam->no_alignment(), 'input bam with alignment');
      $bam->no_alignment(1);
      $bam->clear_no_alignment();
      $bam->no_alignment(0); ## has alignments
      
      $bam->input_bam('t/data/sequence/5551_3#6.bam'); 
      
      my $current_dir = getcwd();
      system "cp -pv $current_dir/t/data/sequence/plasmodium.bam $temp_dir";
      $bam->input_bam("$temp_dir/plasmodium.bam");

      my $bam_pb_cal_cmd = $bam->pb_cal_cmd();

      my $expected_bamseqchk_cmd = $bamseqchksum .
        q{ verbose=0 inputformat=cram reference=t/data/references/Plasmodium_falciparum/default/all/fasta/Pf3D7_v3.fasta};
      is($bam->bamseqchksum_cmd(q{cram}), $expected_bamseqchk_cmd,
        'correct bamseqchksum command for a cram file with reference');

      my $expected_tee_cmd = qq{set -o pipefail;$bammarkduplicates } .
        qq{I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$md_metrics_file | tee};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.md5.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.flagstat.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.stats1.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.stats2.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.bschk.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.alt.bschk.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.index.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.pb_cal.fifo};
      $expected_tee_cmd .= qq{ $temp_dir/output_mk.bam.scramble.fifo};
      $expected_tee_cmd .= qq{ > $temp_dir/output_mk.bam};

      is($bam->_tee_cmd, $expected_tee_cmd, 'entire tee command generated correctly');

      my $cram_seqchksum_file_name_mk = qq{$temp_dir/output_mk.cram.seqchksum};
      my $cram_seqchksum_fifo_name_mk = qq{$temp_dir/output_mk.cram.seqchksum.fifo};
      my $cram_crai_file_name_mk = qq{$temp_dir/output_mk.cram.crai};
      my $cram_crai_fifo_name_mk = qq{$temp_dir/output_mk.cram.crai.fifo};
      my $cram_md5_file_name_mk = qq{$temp_dir/output_mk.cram.md5};
      my $cram_md5_fifo_name_mk = qq{$temp_dir/output_mk.cram.md5.fifo};
      my $cram_file_name_mk = qq{$temp_dir/output_mk.cram};
      my $cram_fifo_name_mk = qq{$temp_dir/output_mk.cram.fifo};
      my $bam_file_name_mk = $output_bam;
      my $bam_seqchksum_file_name_mk = qq{$temp_dir/output_mk.bam.seqchksum};
      my $bam_seqchksum_fifo_name_mk = qq{$temp_dir/output_mk.bam.seqchksum.fifo};

      my @expected_fork_cmds = ();
      my $samtools_cmd = $bam->samtools_cmd;
      my $samtools_irods_cmd = $bam->samtools_irods_cmd();

      my $expected_md5_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.md5.fifo | };
      $expected_md5_cmd .= qq{md5sum -b | tr -d }.q{"\n *-" }. qq{ > $temp_dir/output_mk.bam.md5};

      my $expected_flagstat_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.flagstat.fifo | };
      $expected_flagstat_cmd .= qq{$samtools_cmd flagstat -  > $temp_dir/output_mk.flagstat};

      my $expected_stats1_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.stats1.fifo | };
      $expected_stats1_cmd .= qq{$samtools_irods_cmd stats -F 0x900 > $temp_dir/output_mk_F0x900.stats};
      my $expected_stats2_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.stats2.fifo | };
      $expected_stats2_cmd .= qq{$samtools_irods_cmd stats -F 0xB00 > $temp_dir/output_mk_F0xB00.stats};

      my $expected_index_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.index.fifo | } .
        qq{$samtools_cmd index /dev/stdin /dev/stdout > $temp_dir/output_mk.bai};

      my $expected_pb_cal_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.pb_cal.fifo | };
      $expected_pb_cal_cmd .= qq{$bam_pb_cal_cmd -p $temp_dir/output_mk -filter-bad-tiles 2 -};

      my $expected_bamchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.bschk.fifo | };
      $expected_bamchksum_cmd .= qq{$bamseqchksum verbose=0 inputformat=bam};
      $expected_bamchksum_cmd .= qq{ | tee $temp_dir/output_mk.bam.seqchksum.fifo > $temp_dir/output_mk.bam.seqchksum};
  
      my $expected_altchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_mk.bam.alt.bschk.fifo | };
      $expected_altchksum_cmd .= qq{$bamseqchksum verbose=0 inputformat=bam hash=sha512primesums512};
      $expected_altchksum_cmd .= qq{ > $temp_dir/output_mk.bam.sha512primesums512.seqchksum};
  
      my $expected_scramble_cmd = qq{$scramble -I bam -O cram -r } .
        qq{t/data/references/Plasmodium_falciparum/default/all/fasta/Pf3D7_v3.fasta < $temp_dir/output_mk.bam.scramble.fifo };
      $expected_scramble_cmd .= qq{| tee $temp_dir/output_mk.cram.fifo $cram_crai_fifo_name_mk } .
        qq{$cram_md5_fifo_name_mk > $temp_dir/output_mk.cram};

      my $expected_cramchksum_cmd =  qq{set -o pipefail; cat $cram_fifo_name_mk | $bamseqchksum verbose=0 inputformat=cram } .
        q{reference=t/data/references/Plasmodium_falciparum/default/all/fasta/Pf3D7_v3.fasta };
      $expected_cramchksum_cmd .= qq{| tee $cram_seqchksum_fifo_name_mk > $cram_seqchksum_file_name_mk};

      my $expected_cramindex_cmd = qq{set -o pipefail; cat $cram_crai_fifo_name_mk | $cram_index - $cram_crai_file_name_mk};
      my $expected_crammd5_cmd = qq{set -o pipefail; cat $cram_md5_fifo_name_mk | md5sum -b | tr -d }.q{"\n *-" }. 
        qq{ > $cram_md5_file_name_mk};

      my $expected_diff_cmd = qq{diff $bam_seqchksum_fifo_name_mk $cram_seqchksum_fifo_name_mk};

      push  @expected_fork_cmds, $expected_tee_cmd;
      push  @expected_fork_cmds, $expected_bamchksum_cmd;
      push  @expected_fork_cmds, $expected_altchksum_cmd;
      push  @expected_fork_cmds, $expected_scramble_cmd;
      push  @expected_fork_cmds, $expected_cramchksum_cmd;
      push  @expected_fork_cmds, $expected_cramindex_cmd;
      push  @expected_fork_cmds, $expected_crammd5_cmd;
      push  @expected_fork_cmds, $expected_diff_cmd;
      push  @expected_fork_cmds, $expected_md5_cmd;
      push  @expected_fork_cmds, $expected_flagstat_cmd;
      push  @expected_fork_cmds, $expected_stats1_cmd;
      push  @expected_fork_cmds, $expected_stats2_cmd;
      push  @expected_fork_cmds, $expected_index_cmd;
      push  @expected_fork_cmds, $expected_pb_cal_cmd;

      my $expected_fork_cmds = \@expected_fork_cmds;
      cmp_deeply($bam->fork_cmds(), $expected_fork_cmds, 'commands for fork generated correctly') or
        diag explain [$bam->fork_cmds(),$expected_fork_cmds];

      lives_ok{$bam->process()} q{Processed OK};

      is (!-e "$temp_dir/output_mk.bam.md5.fifo", 1, 'md5 FIFO removed');
      is (!-e "$temp_dir/output_mk.bam.flagstat.fifo", 1, 'flagstat FIFO removed');
      is (!-e "$temp_dir/output_mk.bam.stats1.fifo", 1, 'stats1 FIFO removed');
      is (!-e "$temp_dir/output_mk.bam.stats2.fifo", 1, 'stats2 FIFO removed');
      is (!-e "$temp_dir/output_mk.bam.index.fifo", 1, 'index FIFO removed');
      is (!-e "$temp_dir/output_mk.bam.pb_cal.fifo", 1, 'pb_cal FIFO removed');
      is (!-e "$temp_dir/output_mk.bam.scramble.fifo", 1, 'scramble FIFO removed');
      is (!-e "$temp_dir/output_mk.cram.fifo", 1, 'cram FIFO removed');
      is (!-e "$temp_dir/output_mk.bam.bschk.fifo", 1, 'bamseqchksum input FIFO removed');
      is (!-e "$temp_dir/output_mk.bam.seqchksum.fifo", 1, 'bamseqchksum output FIFO removed');
      is (!-e "$temp_dir/output_mk.cram.seqchksum.fifo", 1, 'bamseqchksum output FIFO removed');

      is (!-z "$temp_dir/output.bam", 1, 'BAM file created with contents');      
      is (!-z "$temp_dir/output.bai", 1, 'BAM index created with contents');      
      is (!-z "$temp_dir/output.bam.md5", 1, 'BAM md5 created with contents');      
      is (!-z "$temp_dir/output.flagstat", 1, 'BAM flagstat created with contents');      
      is (!-z "$temp_dir/output_quality_cycle_caltable.txt", 1, 'Quality caltable created with contents');
      is (!-z "$temp_dir/output_quality_cycle_surv.txt", 1, 'Quality surv created with contents');
      is (!-z "$temp_dir/output_quality_error.txt", 1, 'Quality error table created with contents');
      is (!-z "$temp_dir/output.cram", 1, 'CRAM file created with contents');
      is (!-z "$temp_dir/output.cram.crai", 1, 'CRAM index file created with contents');
      is (!-z "$temp_dir/output.cram.md5", 1, 'CRAM md5 file created with contents');
      is (!-z "$temp_dir/output.seqchksum", 1, 'BAM seqchksum file created with contents');
      is (!-z "$temp_dir/output.sha512primesums512.seqchksum", 1, 'sha512primesums512 seqchksum file created with contents');
      is (!-e "$temp_dir/output_mk.cram.seqchksum", 1, 'CRAM seqchksum file created with contents has been removed');
      #ok(-e "$temp_dir/plasmodium.bam_flagstats.json", 'file with serialized bam_flagstats object exists');
  }    
};

subtest 'subtest 4' => sub {
  my $num_tests = 7;
  plan tests => $num_tests;

  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
            $num_tests unless ($ENV{'TOOLS_INSTALLED'});

      my $temp_dir = tempdir( CLEANUP => 1);
      my $input = join q[/], $temp_dir, 'non_aligned.bam';
      system("cp t/data/sequence/15156_1#54.bam $input");
      my $output_root     = join q[/], $temp_dir, 'non_aligned_output';
      my $output_bam      = $output_root . q[.bam];

      my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
                  input_bam        => $input,
                  output_bam       => $output_bam,
                  metrics_json_dir => $temp_dir,
                  temp_dir         => $temp_dir,
                  id_run           => 35,
                  position         => 1
                );
      my $md_metrics_file = $bam->metrics_file;
      my $expected_mark_duplicate_cmd = qq{$bammarkduplicates I=$input } .
        qq{O=/dev/stdout tmpfile=$temp_dir/ M=$md_metrics_file};
      is($bam->mark_duplicate_cmd(), $expected_mark_duplicate_cmd, 'correct biobambam command');
      is ( $bam->no_alignment(), 1, 'input bam with no alignment');
      lives_ok{$bam->process()} q{Processed OK};
      
      ok (-e "$temp_dir/non_aligned.cram", 'non-aligned CRAM file created');
      ok (-e "$temp_dir/non_aligned.cram.md5", 'non-aligned CRAM md5 file created');
      ok (-e "$temp_dir/non_aligned.seqchksum", 'non-aligned BAM seqchksum file created');
      ok (-e "$temp_dir/non_aligned.sha512primesums512.seqchksum",
       'non-aligned sha512primesums512 seqchksum file created');
      #ok (-e "$temp_dir/non_aligned.bam_flagstats.json", 'file with serialized bam_flagstats object exists');
  }
};

subtest 'subtest 5' => sub {
  my $num_tests = 16;
  plan tests => $num_tests;

  SKIP: {
      skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
            $num_tests unless ($ENV{'TOOLS_INSTALLED'});

      my $temp_dir = tempdir( CLEANUP => 1);
      my $input = join q[/], $temp_dir, 'no_align.bam';
      system("cp t/data/sequence/unaligned.bam $input");
      my $output_root     = join q[/], $temp_dir, 'output_no_align';
      my $output_bam      = $output_root . q[.bam];
      my $md_metrics_file = $output_root . '.markdups_metrics.txt';
        
      my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
                 input_bam        => $input,
                 output_bam       => $output_bam,
                 metrics_json_dir => $temp_dir,
                 temp_dir         => $temp_dir,
                 metrics_file     => $md_metrics_file,
                 no_alignment     => 1,
                 id_run           => 1234,
                 position         => 1,
               );
      my $expected_mark_duplicate_cmd = qq{$bammarkduplicates I=$input O=/dev/stdout tmpfile=$temp_dir/ M=$md_metrics_file};
      is($bam->mark_duplicate_cmd(), $expected_mark_duplicate_cmd, 'correct biobambam command');
      ok( $bam->no_alignment(), 'input bam without alignment');
      like($bam->bamseqchksum_cmd(q{bam}), qr{\Qbamseqchksum verbose=0 inputformat=bam\E}, 'correct bamseqchksum command for a bam file with reference but no alignment');
      like($bam->bamseqchksum_cmd(q{cram}), qr{\Qbamseqchksum verbose=0 inputformat=cram\E}, 'correct bamseqchksum command for a cram file with reference but no alignment');

      my $mdup_cmd = $bam->mark_duplicate_cmd;
      my $samtools_cmd = $bam->samtools_cmd;
      my $samtools_irods_cmd = $bam->samtools_irods_cmd();

      my $expected_tee_cmd = qq{set -o pipefail;$mdup_cmd | tee};
      $expected_tee_cmd .= qq{ $temp_dir/output_no_align.bam.md5.fifo $temp_dir/output_no_align.bam.flagstat.fifo $temp_dir/output_no_align.bam.stats1.fifo $temp_dir/output_no_align.bam.stats2.fifo $temp_dir/output_no_align.bam.bschk.fifo $temp_dir/output_no_align.bam.alt.bschk.fifo $temp_dir/output_no_align.bam.scramble.fifo > $temp_dir/output_no_align.bam};
      is($bam->_tee_cmd, $expected_tee_cmd, 'entire tee command generated correctly if no_alignment flag used');

      my @expected_fork_cmds = ();

      my $expected_md5_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.md5.fifo | };
      $expected_md5_cmd .= qq{md5sum -b | tr -d }.q{"\n *-" }. qq{ > $temp_dir/output_no_align.bam.md5};

      my $expected_flagstat_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.flagstat.fifo | };
      $expected_flagstat_cmd .= qq{$samtools_cmd flagstat -  > $temp_dir/output_no_align.flagstat};

      my $expected_stats1_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.stats1.fifo | };
      $expected_stats1_cmd .= qq{$samtools_irods_cmd stats -F 0x900 > $temp_dir/output_no_align_F0x900.stats};
      my $expected_stats2_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.stats2.fifo | };
      $expected_stats2_cmd .= qq{$samtools_irods_cmd stats -F 0xB00 > $temp_dir/output_no_align_F0xB00.stats};

      my $expected_altchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.alt.bschk.fifo | };
      $expected_altchksum_cmd .= qq{$bamseqchksum verbose=0 inputformat=bam hash=sha512primesums512};
      $expected_altchksum_cmd .= qq{ > $temp_dir/output_no_align.bam.sha512primesums512.seqchksum};

      my $cram_md5_file_name_mk = qq{$temp_dir/output_no_align.cram.md5};
      my $cram_md5_fifo_name_mk = qq{$temp_dir/output_no_align.cram.md5.fifo};
      my $cram_seqchksum_file_name_mk = qq{$temp_dir/output_no_align.cram.seqchksum};
      my $cram_seqchksum_fifo_name_mk = qq{$temp_dir/output_no_align.cram.seqchksum.fifo};
  
      my $expected_scramble_cmd = qq{$scramble -I bam -O cram < $temp_dir/output_no_align.bam.scramble.fifo };
      $expected_scramble_cmd .= qq{| tee $temp_dir/output_no_align.cram.fifo $cram_md5_fifo_name_mk > $temp_dir/output_no_align.cram};

      my $expected_cramchksum_cmd =  qq{set -o pipefail; cat $temp_dir/output_no_align.cram.fifo | $bamseqchksum verbose=0 inputformat=cram };
      $expected_cramchksum_cmd .= qq{| tee $cram_seqchksum_fifo_name_mk > $cram_seqchksum_file_name_mk};

      my $expected_crammd5_cmd = qq{set -o pipefail; cat $cram_md5_fifo_name_mk | md5sum -b | tr -d }.q{"\n *-" }. qq{ > $cram_md5_file_name_mk};

      my $expected_bamchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_no_align.bam.bschk.fifo | };
      $expected_bamchksum_cmd .= qq{$bamseqchksum verbose=0 inputformat=bam};
      $expected_bamchksum_cmd .= qq{ | tee $temp_dir/output_no_align.bam.seqchksum.fifo};
      $expected_bamchksum_cmd .= qq{ > $temp_dir/output_no_align.bam.seqchksum};
  
      my $expected_diff_cmd = qq{diff $temp_dir/output_no_align.bam.seqchksum.fifo $cram_seqchksum_fifo_name_mk};

      push  @expected_fork_cmds, $expected_tee_cmd;
      push  @expected_fork_cmds, $expected_bamchksum_cmd;
      push  @expected_fork_cmds, $expected_altchksum_cmd;
      push  @expected_fork_cmds, $expected_scramble_cmd;
      push  @expected_fork_cmds, $expected_cramchksum_cmd;
      push  @expected_fork_cmds, $expected_crammd5_cmd;
      push  @expected_fork_cmds, $expected_diff_cmd;
      push  @expected_fork_cmds, $expected_md5_cmd;
      push  @expected_fork_cmds, $expected_flagstat_cmd;
      push  @expected_fork_cmds, $expected_stats1_cmd;
      push  @expected_fork_cmds, $expected_stats2_cmd;

      my $expected_fork_cmds = \@expected_fork_cmds;
      cmp_deeply($bam->fork_cmds(), $expected_fork_cmds, 'commands for ForkManager generated correctly') or
        diag explain [$bam->fork_cmds(),$expected_fork_cmds];

      lives_ok{$bam->process()} q{Processed OK};

      is (!-z "$temp_dir/no_align.bam", 1, 'BAM file created with contents');      
      is (!-e "$temp_dir/no_align.bai", 1, 'BAM index NOT created if no_alignment flag used');      
      is (!-z "$temp_dir/no_align.bam.md5", 1, 'BAM md5 created with contents');      
      is (!-z "$temp_dir/no_align.flagstat", 1, 'BAM flagstat created with contents');      
      is (!-z "$temp_dir/no_align_quality_cycle_caltable.txt", 1, 'Quality caltable created with contents');
      is (!-z "$temp_dir/no_align_quality_cycle_surv.txt", 1, 'Quality surv created with contents');
      is (!-z "$temp_dir/no_align_quality_error.txt", 1, 'Quality error table created with contents');
      is (-e "$temp_dir/no_align.cram", 1, 'CRAM file created if no_alignment flag used');
      is (!-z "$temp_dir/no_align.bam.seqchksum", 1, 'BAM seqchksum file created with contents if no_alignment flag used');
      #ok (-e "$temp_dir/no_align.bam_flagstats.json", 'file with serialized bam_flagstats object exists');
  }
};

subtest 'subtest 6' => sub {
  my $num_tests = 19;
  plan tests => $num_tests;

  SKIP: {
    skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
       $num_tests unless ($ENV{'TOOLS_INSTALLED'});

    my $temp_dir = tempdir( CLEANUP => 1);
    my $input = join q[/], $temp_dir, 'phix.bam';
    system("cp t/data/sequence/phix.bam $input");
    my $output_root     = join q[/], $temp_dir, 'output_phix';
    my $output_bam      = $output_root . q[.bam];
    my $md_metrics_file = $output_root . '.markdups_metrics.txt';

    my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
                 input_bam        => $input,
                 output_bam       => $output_bam,
                 metrics_json_dir => $temp_dir,
                 temp_dir         => $temp_dir,
                 metrics_file     => $md_metrics_file,
                 subset           => 'phix',
                 id_run           => 1234,
                 position         => 2,
               );
      my $expected_mark_duplicate_cmd = qq{$bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$md_metrics_file};
      is($bam->mark_duplicate_cmd(), $expected_mark_duplicate_cmd, 'correct biobambam command');
      ok(!$bam->no_alignment(), 'input PhiX bam with alignment');      
      like($bam->bamseqchksum_cmd(q{bam}), qr{\Qbamseqchksum verbose=0 inputformat=bam\E}, 'correct bamseqchksum command for a bam file with reference but no alignment');
      like($bam->bamseqchksum_cmd(q{cram}), qr{\Qbamseqchksum verbose=0 inputformat=cram\E}, 'correct bamseqchksum command for a cram file with reference but no alignment');

      my $samtools_cmd = $bam->samtools_cmd;
      my $samtools_irods_cmd = $bam->samtools_irods_cmd();

      my $expected_tee_cmd = qq{set -o pipefail;$bammarkduplicates I=$temp_dir/sorted.bam O=/dev/stdout tmpfile=$temp_dir/ M=$md_metrics_file | tee};
      $expected_tee_cmd .= qq{ $temp_dir/output_phix.bam.md5.fifo $temp_dir/output_phix.bam.flagstat.fifo $temp_dir/output_phix.bam.stats1.fifo $temp_dir/output_phix.bam.stats2.fifo $temp_dir/output_phix.bam.bschk.fifo $temp_dir/output_phix.bam.alt.bschk.fifo $temp_dir/output_phix.bam.index.fifo $temp_dir/output_phix.bam.pb_cal.fifo $temp_dir/output_phix.bam.scramble.fifo > $temp_dir/output_phix.bam};
      is($bam->_tee_cmd, $expected_tee_cmd, 'entire tee command generated correctly for PhiX');

      my @expected_fork_cmds = ();

      my $expected_md5_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.md5.fifo | };
      $expected_md5_cmd .= qq{md5sum -b | tr -d }.q{"\n *-" }. qq{ > $temp_dir/output_phix.bam.md5};

      my $expected_flagstat_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.flagstat.fifo | };
      $expected_flagstat_cmd .= qq{$samtools_cmd flagstat -  > $temp_dir/output_phix.flagstat};

      my $expected_stats1_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.stats1.fifo | };
      $expected_stats1_cmd .= qq{$samtools_irods_cmd stats -F 0x900 > $temp_dir/output_phix_F0x900.stats};
      my $expected_stats2_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.stats2.fifo | };
      $expected_stats2_cmd .= qq{$samtools_irods_cmd stats -F 0xB00 > $temp_dir/output_phix_F0xB00.stats};

      my $expected_index_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.index.fifo | $samtools_cmd index /dev/stdin /dev/stdout > $temp_dir/output_phix.bai};

      my $bam_pb_cal_cmd = $bam->pb_cal_cmd();
      my $expected_pb_cal_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.pb_cal.fifo | };
      $expected_pb_cal_cmd .= qq{$bam_pb_cal_cmd -p $temp_dir/output_phix -filter-bad-tiles 2 -};

      my $expected_bamchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.bschk.fifo | };
      $expected_bamchksum_cmd .= qq{$bamseqchksum verbose=0 inputformat=bam};
      $expected_bamchksum_cmd .= qq{ | tee $temp_dir/output_phix.bam.seqchksum.fifo};
      $expected_bamchksum_cmd .= qq{ > $temp_dir/output_phix.bam.seqchksum};
 
      my $cram_seqchksum_file_name_phix = qq{$temp_dir/output_phix.cram.seqchksum};
      my $cram_seqchksum_fifo_name_phix = qq{$temp_dir/output_phix.cram.seqchksum.fifo};
      my $cram_crai_file_name_phix = qq{$temp_dir/output_phix.cram.crai};
      my $cram_crai_fifo_name_phix = qq{$temp_dir/output_phix.cram.crai.fifo};
      my $cram_md5_file_name_phix = qq{$temp_dir/output_phix.cram.md5};
      my $cram_md5_fifo_name_phix = qq{$temp_dir/output_phix.cram.md5.fifo};
      my $cram_file_name_phix = qq{$temp_dir/output_phix.cram};
      my $cram_fifo_name_phix = qq{$temp_dir/output_phix.cram.fifo};

      my $expected_altchksum_cmd = qq{set -o pipefail; cat $temp_dir/output_phix.bam.alt.bschk.fifo | };
      $expected_altchksum_cmd .= qq{$bamseqchksum verbose=0 inputformat=bam hash=sha512primesums512};
      $expected_altchksum_cmd .= qq{ > $temp_dir/output_phix.bam.sha512primesums512.seqchksum};
  
      my $expected_scramble_cmd = qq{$scramble -I bam -O cram < $temp_dir/output_phix.bam.scramble.fifo };
      $expected_scramble_cmd .= qq{| tee $temp_dir/output_phix.cram.fifo $cram_crai_fifo_name_phix $cram_md5_fifo_name_phix > $temp_dir/output_phix.cram};

      my $expected_cramchksum_cmd =  qq{set -o pipefail; cat $cram_fifo_name_phix | $bamseqchksum verbose=0 inputformat=cram };
      $expected_cramchksum_cmd .= qq{| tee $cram_seqchksum_fifo_name_phix > $cram_seqchksum_file_name_phix};

      my $expected_cramindex_cmd = qq{set -o pipefail; cat $cram_crai_fifo_name_phix | $cram_index - $cram_crai_file_name_phix};
      my $expected_crammd5_cmd = qq{set -o pipefail; cat $cram_md5_fifo_name_phix | md5sum -b | tr -d }.q{"\n *-" }. qq{ > $cram_md5_file_name_phix};

      my $expected_diff_cmd = qq{diff $temp_dir/output_phix.bam.seqchksum.fifo $cram_seqchksum_fifo_name_phix};

      push  @expected_fork_cmds, $expected_tee_cmd;
      push  @expected_fork_cmds, $expected_bamchksum_cmd;
      push  @expected_fork_cmds, $expected_altchksum_cmd;
      push  @expected_fork_cmds, $expected_scramble_cmd;
      push  @expected_fork_cmds, $expected_cramchksum_cmd;
      push  @expected_fork_cmds, $expected_cramindex_cmd;
      push  @expected_fork_cmds, $expected_crammd5_cmd;
      push  @expected_fork_cmds, $expected_diff_cmd;
      push  @expected_fork_cmds, $expected_md5_cmd;
      push  @expected_fork_cmds, $expected_flagstat_cmd;
      push  @expected_fork_cmds, $expected_stats1_cmd;
      push  @expected_fork_cmds, $expected_stats2_cmd;
      push  @expected_fork_cmds, $expected_index_cmd;
      push  @expected_fork_cmds, $expected_pb_cal_cmd;

      my $expected_fork_cmds = \@expected_fork_cmds;
      cmp_deeply($bam->fork_cmds(), $expected_fork_cmds, 'commands for ForkManager generated correctly') or diag explain [$bam->fork_cmds(),$expected_fork_cmds];

      lives_ok{$bam->process()} q{Processed OK};

      is (!-z "$temp_dir/phix.bam", 1, 'BAM file created with contents for PhiX');      
      is (!-z "$temp_dir/phix.bai", 1, 'BAM index created with contents for PhiX');      
      is (!-z "$temp_dir/metrics_phix.bam.json", 1, 'metrics json created with contents for PhiX');      
      is (!-z "$temp_dir/phix.markdups_metrics.txt", 1, 'metrics txt created with contents for PhiX');      
      is (!-z "$temp_dir/phix.bam.md5", 1, 'BAM md5 created with contents for PhiX');      
      is (!-z "$temp_dir/phix.flagstat", 1, 'BAM flagstat created with contents for PhiX');      
      is (!-z "$temp_dir/phix_quality_cycle_caltable.txt", 1, 'Quality caltable created with contents for PhiX');
      is (!-z "$temp_dir/phix_quality_cycle_surv.txt", 1, 'Quality surv created with contents for PhiX');
      is (!-z "$temp_dir/phix_quality_error.txt", 1, 'Quality error table created with contents for PhiX');

      is (!-z "$temp_dir/phix.cram", 1, 'CRAM file created with contents for PhiX');
      is (!-z "$temp_dir/phix.cram.seqchksum", 1, 'CRAM seqchksum file created with contents for PhiX');
      is (!-z "$temp_dir/phix.bam.seqchksum", 1, 'BAM seqchksum file created with contents for PhiX');
      #ok(-e "$temp_dir/phix_phix.bam_flagstats.json", 'file with serialized bam_flagstats object exists');
  }
};

1;
