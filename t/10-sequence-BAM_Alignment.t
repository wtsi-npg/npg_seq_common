use strict;
use warnings;
use English qw{-no_match_vars};
use Test::More tests => 107;
use Test::Deep;
use Test::Exception;
use File::Temp qw(tempfile tempdir);
use File::Which qw[which];
use Perl6::Slurp;
use Cwd qw(getcwd abs_path);
use npg_tracking::data::reference::list;

use npg_qc::autoqc::results::bam_flagstats;
use Test::MockModule;
my $mock = Test::MockModule->new('npg_qc::autoqc::results::bam_flagstats');
$mock->mock(execute  => sub { return 1; });
$mock->mock(set_info => sub { return 1; });
$mock->mock(store    => sub { return 1; });

local $ENV{'http_proxy'}='http://wibble.com';

use_ok('npg_common::sequence::BAM_Alignment');


my $REP_ROOT = $npg_tracking::data::reference::list::REP_ROOT;

my $temp_dir = tempdir( CLEANUP => 1 );
my $current_dir = getcwd();
my $java_memory = q[-Xmx100m];
my $JAVA_CMD = q[java];
my $human = q[/my/ref/repos/Homo_sapiens/strain];
my $test_classpath = q[t/bin/aligners/picard/current:t/bin/aligners/illumina2bam/current];

{
  local $ENV{CLASSPATH} = $test_classpath;
  my $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 'input.bam',
                 is_paired_read    => 1,
                 non_consent_split => 1,
                 spiked_phix_split => 1,       
                 read_length       => 37,
                 output_prefix     => 'output',
                 bwa_cmd           => 't/bin/aligners/bwa/bwa-0.5.8c/bwa',
                 samtools_cmd      => 't/bin/aligners/samtools/current/samtools',
                 temp_dir          => $temp_dir,
                 reference         => 'reference.fasta',
                 java_xmx_flag     => $java_memory,
                 id_run            => 35,
                 position          => 1
               );
               
  isa_ok( $bam, 'npg_common::sequence::BAM_Alignment', 'object test' );
  
  ok( ! $bam->no_alignment(), 'output bam will be created with alignment' );
  
  is($bam->output(), 'output.bam', 'correct output file name');
  is($bam->non_consented_output(),'output_human.bam', 'correct output file name for human');
  is($bam->spiked_phix_output(),'output_phix.bam', 'correct output file name for spiked phix');
  is($bam->_bam_output_fifo(), $temp_dir.q{/output_fifo.bam}, 'correct bam output fifo file name');
  is($bam->_human_bam_output_fifo(), $temp_dir.q{/human_output_fifo.bam}, 'correct human bam output fifo file name');
  
  my $rel_aligner_path = q[t/bin/bwa];
  my $expected_aligner_cmd = join q[/], getcwd(), $rel_aligner_path;
  is ( $bam->bwa_cmd(), $expected_aligner_cmd, 'correct aligner command with absolute path' );

  my $rel_samtools_path = q[t/bin/aligners/samtools/current/samtools];
  my $abs_samtools_path = abs_path($rel_samtools_path);

  is ($bam->samtools_cmd(), $abs_samtools_path , 'correct samtools command with absolute path' );
  is( $bam->bwa_aln_options(), '-q 15 -t 4', 'default bwa aln options');
  
  like( $bam->_generate_bwa_aln_command('reference.fasta', 'input.bam', 1),
        qr{t/bin/bwa aln -q 15 -t 4 reference.fasta -b1 input.bam > $temp_dir/1.sai$},
        'correct bwa aln command for read 1'
      );
  like( $bam->_generate_bwa_aln_command('human_reference.fasta', 'input.bam', 0, 1),
        qr{t/bin/bwa aln -q 15 -t 4 human_reference.fasta -b0 input.bam > $temp_dir/human0.sai$},
        'correct bwa aln human command for single read'
      );
  like ( $bam->_generate_bwa_sam_command('reference.fasta', ['input.bam','input.bam']),
        qr{t/bin/bwa sampe -t 1 reference.fasta $temp_dir/1.sai $temp_dir/2.sai input.bam input.bam$},
         'correct bwa sampe command for paired read'
       );
  like ( $bam->_generate_bwa_sam_command('reference.fasta', ['input.bam'], 1),
         qr{t/bin/bwa samse -t 1 reference.fasta $temp_dir/human0.sai input.bam$},
         'correct bwa samse command for single read and human split alignment'
       );
  $bam->bwa_sampe_options('-s');
  like ( $bam->_generate_bwa_sam_command('reference.fasta', ['input.bam','input.bam'], 1),
        qr{t/bin/bwa sampe -t 1 -s reference.fasta $temp_dir/human1.sai $temp_dir/human2.sai input.bam input.bam$},
         'correct bwa sampe command for paired read and human split with -s option'
       );    
  is( scalar @{$bam->_bwa_aln_commands()}, 2, '2 bwa aln commands for paired-read input file');
  $bam->human_reference('human_reference.fasta');
  is( scalar @{$bam->_human_bwa_aln_commands()}, 2, '2 bwa aln commands for paired-read input file against human');
  like ( $bam->bwa_sam_command(),
         qr{t/bin/bwa sampe -t 1 -s reference.fasta $temp_dir/1.sai $temp_dir/2.sai input.bam input.bam$},
         'correct bwa sampe command for paired read with -s option'
       );
  like ( $bam->_human_bwa_sam_command(),
         qr{t/bin/bwa sampe -t 1 -s human_reference.fasta $temp_dir/human1.sai $temp_dir/human2.sai input.bam input.bam$},
         'correct human bwa sampe command for paired read with -s option'
       );
  
  my @bwa_sam_commands = sort keys %{$bam->_bwa_sam_commands()};
  is(scalar @bwa_sam_commands, 4, '4 bwa sam commands returned inclusing cat - tee commands');
  is($bwa_sam_commands[2], "cat input.bam | tee $temp_dir/fifo1.bam $temp_dir/human_fifo1.bam > /dev/null", 'first cat tee command');
  is($bwa_sam_commands[3], "cat input.bam | tee $temp_dir/fifo2.bam $temp_dir/human_fifo2.bam > /dev/null", 'second cat tee command');
  like ( $bam->bwa_sam_command(),
         qr{t/bin/bwa sampe -t 1 -s reference.fasta $temp_dir/1.sai $temp_dir/2.sai $temp_dir/fifo1.bam $temp_dir/fifo2.bam$},
         'correct bwa sampe command for paired read with -s option using fifo'
       );
  like ( $bam->_human_bwa_sam_command(),
         qr{t/bin/bwa sampe -t 1 -s human_reference.fasta $temp_dir/human1.sai $temp_dir/human2.sai $temp_dir/human_fifo1.bam $temp_dir/human_fifo2.bam$},
         'correct human bwa sampe command for paired read with -s option using fifo'
       );
  ok( -e "$temp_dir/fifo1.bam", 'first fifo created');
  ok( -e "$temp_dir/fifo2.bam", 'second fifo created');
  ok( -e "$temp_dir/human_fifo1.bam", 'first human fifo created');
  ok( -e "$temp_dir/human_fifo2.bam", 'second human fifo created');

  my $illumina2bam_path = abs_path( q[t/bin/aligners/illumina2bam/current] );
  my $expected_alignment_filter_jar = qq[$illumina2bam_path/AlignmentFilter.jar];
  my $expected_filter_cmd = qq{$JAVA_CMD $java_memory -jar $expected_alignment_filter_jar VALIDATION_STRINGENCY=SILENT CREATE_MD5_FILE=true METRICS_FILE=output.bam_alignment_filter_metrics.json IN=input.bam IN=${temp_dir}/human_output_fifo.bam IN=${temp_dir}/output_fifo.bam OUT=output_phix.bam OUT=output_human.bam OUT=output.bam};
  like( $bam->_bam_alignment_filter_command(), qr[$expected_filter_cmd], 'correct filter command with absolute jar path');
  ok( -e "$temp_dir/output_fifo.bam", 'bam output fifo created before alignment filter');
  ok( -e "$temp_dir/human_output_fifo.bam", 'human bam output fifo created before alignment filter');

  my $expected_bam_merger_jar = qq[$illumina2bam_path/BamMerger.jar];
  my $expected_bam_merger_cmd = qq{$JAVA_CMD $java_memory -jar $expected_bam_merger_jar};
  like( $bam->_bam_merger_cmd(), qr[$expected_bam_merger_cmd], "correct bam merger command with absolute jar path");

  my $picard_path = abs_path( q[t/bin/aligners/picard/current] ); 
  my $expected_sam_converter_jar = $picard_path . q[/SamFormatConverter.jar];
  my $expected_picard_converter_cmd = qq[$JAVA_CMD $java_memory -jar $expected_sam_converter_jar];
  like ($bam->_picard_converter_cmd(), qr[$expected_picard_converter_cmd], "correct Sam format converter command with absolute jar path");
}

{
  local $ENV{CLASSPATH} = $test_classpath;
  local $ENV{LSB_MCPU_HOSTS} = q(sf-3-3-12 6);
  my $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 'input.bam',
                 is_paired_read    => 1,
                 non_consent_split => 1,
                 spiked_phix_split => 1,
                 read_length       => 37,
                 output_prefix     => 'output',
                 bwa_cmd           => 't/bin/bwa',
                 samtools_cmd      => 't/bin/aligners/samtools/current/samtools',
                 temp_dir          => $temp_dir,
                 reference         => 'reference.fasta',
                 id_run            => 35,
                 position          => 1
               );

  is( $bam->bwa_cmd(), join(q[/], getcwd(), q[t/bin/bwa]), 'correct aligner command with absolute path' );
  is( $bam->bwa_aln_options(), '-q 15 -t 6', 'default (LSF aware) bwa aln options');

  like( $bam->_generate_bwa_aln_command('reference.fasta', 'input.bam', 1),
        qr{t/bin/bwa aln -q 15 -t 6 reference.fasta -b1 input.bam > $temp_dir/1.sai$},
        'correct (LSF aware) bwa aln command for read 1'
      );
}

{
  local $ENV{CLASSPATH} = $test_classpath;
  my $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 'input.bam',
                 is_paired_read    => 0,
                 non_consent_split => 0,
                 spiked_phix_split => 0,       
                 read_length       => 24,
                 output_prefix     => 'output',
                 temp_dir          => $temp_dir,
                 bwa_cmd           => 't/bin/aligners/bwa/bwa-0.5.8c/bwa',
                 id_run            => 35,
                 position          => 1
              );
               
  isa_ok($bam, 'npg_common::sequence::BAM_Alignment', 'object test');

  ok( $bam->no_alignment(), 'read too short no alignment for output' );
  
  $bam->clear_no_alignment();  
  $bam->read_length(76);
  ok( $bam->no_alignment(), 'no reference given, no alignment for output' );
  
  $bam->clear_no_alignment();
  $bam->reference('reference.fasta');
  $bam->human_reference('human_reference.fasta');
  is( scalar @{$bam->_bwa_aln_commands()}, 1, '1 bwa aln command for single-read input file');
  is( scalar @{$bam->_human_bwa_aln_commands()}, 1, '1 bwa aln command for single-read input file against human');

  like ( $bam->bwa_sam_command(),
         qr{t/bin/bwa samse -t 4 reference.fasta $temp_dir/0.sai input.bam$},
         'correct bwa samse command for single read'
       );
  like ( $bam->_human_bwa_sam_command(),
         qr{t/bin/bwa samse -t 4 human_reference.fasta $temp_dir/human0.sai input.bam$},
         'correct human bwa samse command for single read'
       );

  my @bwa_sam_commands = sort keys %{$bam->_bwa_sam_commands()};
  is(scalar @bwa_sam_commands, 1, 'only one bwa sam command returned without non consented splitting');
  ok(!$bam->_bam_alignment_filter_command(), 'no alignment filter required');
}

{ 
  local $ENV{CLASSPATH} = $test_classpath;
  my $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 'input.bam',
                 is_paired_read    => 0,
                 non_consent_split => 1,
                 spiked_phix_split => 1,       
                 output_prefix     => 'output',
                 temp_dir          => $temp_dir,
                 bwa_cmd           => 't/bin/aligners/bwa/bwa-0.5.8c/bwa',
                 reference         => 'reference.fasta',
                 human_reference   => 'human_reference.fasta',
                 id_run            => 35,
                 position          => 1
               );
               
  isa_ok($bam, 'npg_common::sequence::BAM_Alignment', 'object test');

  ok( !$bam->no_alignment(), 'default with alignment when reference given' );

  my @bwa_sam_commands = sort keys %{$bam->_bwa_sam_commands()};
  is(scalar @bwa_sam_commands, 3, '3 bwa sam command returned including 1 cat-tee command');

  ok( -e "$temp_dir/fifo0.bam", 'one fifo created');
  ok( -e "$temp_dir/human_fifo0.bam", 'one human fifo created');
  is($bwa_sam_commands[2], "cat input.bam | tee $temp_dir/fifo0.bam $temp_dir/human_fifo0.bam > /dev/null", 'one cat tee command');
  like ( $bam->bwa_sam_command(),
         qr{t/bin/bwa samse -t 1 reference.fasta $temp_dir/0.sai $temp_dir/fifo0.bam$},
         'correct bwa samse command for single read using fifo'
       );
  like ( $bam->_human_bwa_sam_command(),
         qr{t/bin/bwa samse -t 1 human_reference.fasta $temp_dir/human0.sai $temp_dir/human_fifo0.bam$},
         'correct human bwa samse command for single read using fifo'
       );
}

my $skip_message = 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run. ';

SKIP: {

  skip $skip_message, 12 if (!$ENV{'TOOLS_INSTALLED'});

  { 
    my $bam;
    my $jar_error;

    $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',
                 is_paired_read    => 1,
                 non_consent_split => 1,
                 spiked_phix_split => 1,       
                 output_prefix     => 'output',
                 temp_dir          => $temp_dir,
                 reference         => 'reference.fasta',
                 human_reference   => 'human_reference.fasta',
                 alignment_metrics_autoqc_command => undef,
                 java_xmx_flag => $java_memory,
                 id_run        => 35,
                 position      => 1,
              );

      my $bwa_aln_pg_lines = $bam->_bwa_aln_pg_lines();
      my $human_bwa_aln_pg_lines = $bam->_human_bwa_aln_pg_lines();

      is( @{ $bam->_original_pg_lines()}, 8, 'correct number PG lines' );
      like ($bwa_aln_pg_lines->[0], qr{^\@PG\tID\:bwa_aln\tPN:bwa\tPP:SplitBamByReadGroup}, 'correct first bwa aln pg');
      like ($human_bwa_aln_pg_lines->[0], qr{^\@PG\tID\:bwa_aln\tPN:bwa\tPP:SplitBamByReadGroup}, 'correct first human bwa aln pg');
      like ($bwa_aln_pg_lines->[1], qr{^\@PG\tID\:bwa_aln_1\tPN:bwa\tPP:bwa_aln}, 'correct second bwa aln pg');
      like ($human_bwa_aln_pg_lines->[1], qr{^\@PG\tID\:bwa_aln_1\tPN:bwa\tPP:bwa_aln}, 'correct human second bwa aln pg');
   
      my $output_bam = $temp_dir.'/test_out.bam';
      my $pg_records = [];
      my $output_bam_fh;
      lives_ok { ($output_bam_fh, $pg_records) = 
                  $bam->generate_output_bam_fh($output_bam, $pg_records, $pg_records)} 'generate output bam stream';
      is( @{$pg_records}, 2, '2 PG returned');
      like ($pg_records->[0], qr{^\@PG\tID\:Picard_SamFormatConverter\tPN\:SamFormatConverter\tVN\:\d\.\d+\(\d+\)}, 'first program is Picard SamFormatConverter');
      like ($pg_records->[1], qr{^\@PG\tID\:samtools_fixmate\tPN:samtools}, 'second program is samtools fixmate');
      isa_ok( $output_bam_fh,  'GLOB', 'output file handle');
      close $output_bam_fh;

      my $command = 'cat t/data/sequence/5431_6#3.sam';
      $output_bam = "$temp_dir/test.bam";
      my $output_bam_md5 = $output_bam.q{.md5};
      lives_ok{ $bam->_bwa_sam_out_to_bam($command, $bwa_aln_pg_lines, $output_bam)} 'covert sam stream from a command to a bam with steps: samtools fixmate and merge input with bam out put';
      ok(-e $output_bam, 'bam generated');
  }
} # end skip no TOOLS_INSTALLED

{
  local $ENV{CLASSPATH} = $test_classpath;
  my $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',  
                 output_prefix     => 'output',
                 temp_dir          => $temp_dir,
                 reference         => 'reference.fasta',
                 human_reference   => 'human_reference.fasta',
                 alignment_metrics_autoqc_command => undef,
                 id_run            => 35,
                 position          => 1,
               );
  
  ok(! $bam->_lims(), 'no id_run given, no lims object');
  ok(!$bam->is_paired_read(), 'is_paired_read is false when no id_run given');
  ok(! $bam->spiked_phix_split(), 'spiked_phix_split is false when no id_run given');
  ok(!$bam->non_consent_split(), 'clear_non_consent_split is false when no id_run given');
}

my $cache = q[t/data/sequence];
local $ENV{NPG_WEBSERVICE_CACHE_DIR} = $cache;

{  
   SKIP: {
     skip "reference repository root $REP_ROOT is not accessible", 14 if !-d $REP_ROOT;

     local $ENV{CLASSPATH} = $test_classpath;
   
     my $bam = npg_common::sequence::BAM_Alignment->new(
                input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 'output',
                 temp_dir          => $temp_dir,
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 1,
                 alignment_metrics_autoqc_command => undef,
               );
      isa_ok($bam->_lims(), 'st::api::lims', 'lims object created');

      ok($bam->is_paired_read(), 'run 5175 is paired read run');
      ok(! $bam->spiked_phix_split(), 'run 5175 lane 1 and plex 1 no need spiked phix split');
      ok(! $bam->non_consent_split(), 'run 5175 lane 1 and plex 1 no non-consented split');

      ok(!$bam->no_alignment(), 'run 5175 lane 1 and plex 1 needs to be aligned');
   
      lives_ok { $bam->reference_info() } 'getting reference info';
      lives_ok { $bam->reference_info() } 'getting reference info';
      is ( $bam->reference(), $REP_ROOT . 'references/Bordetella_bronchiseptica/RB50/all/bwa/B_bronchiseptica_RB50.fasta', 'reference for run 5175 lane 1 plex 1');
      is( $bam->bwa_aln_options(),'-q 15 -t 4', 'correct bwa aln option');
      is( $bam->ref_dict(), $REP_ROOT . 'references/Bordetella_bronchiseptica/RB50/all/picard/B_bronchiseptica_RB50.fasta.dict', 'correct ref dict');
      
      is( $bam->human_reference(), $REP_ROOT . 'references/Homo_sapiens/1000Genomes/all/bwa/human_g1k_v37.fasta', 'correct human reference');
      is( $bam->human_ref_dict(), $REP_ROOT . 'references/Homo_sapiens/1000Genomes/all/picard/human_g1k_v37.fasta.dict', 'correct human reference dict');

      is( $bam->_generate_markduplicates_cmd("$temp_dir/input.bam"), "$current_dir/t/bam_mark_duplicate.pl --input_bam $temp_dir/input.bam --id_run 5175 --position 1 --tag_index 1 --output_bam $temp_dir/input_mk.bam --metrics_json_dir $temp_dir/qc", "correct markduplicates command");   

      is( $bam->_generate_markduplicates_cmd("$temp_dir/input_phix.bam", 'phix'), "$current_dir/t/bam_mark_duplicate.pl --input_bam $temp_dir/input_phix.bam --id_run 5175 --position 1 --tag_index 1 --subset phix --output_bam $temp_dir/input_phix_mk.bam --metrics_json_dir $temp_dir/qc", "correct markduplicates command for phix part bam"); # 88
   }
}

{
  SKIP: {
    skip "reference repository root $REP_ROOT is not available", 2 if !-d $REP_ROOT;

    local $ENV{CLASSPATH} = $test_classpath;
    my $bam = npg_common::sequence::BAM_Alignment->new(
                input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 'output',
                 temp_dir          => $temp_dir,
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 0,
                 alignment_metrics_autoqc_command => undef,
               );
    ok(!$bam->no_alignment(), 'run 5175 lane 1 plex 0 will be aligned');

    $bam = npg_common::sequence::BAM_Alignment->new(
               { input             => 't/data/sequence/6062_1.bam',     
                 output_prefix     => 'output',
                 temp_dir          => $temp_dir,
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 alignment_metrics_autoqc_command => undef,
               });
    ok(!$bam->no_alignment(), 'run 5175 lane 1 will be aligned');
  }
}

{
  local $ENV{CLASSPATH} = $test_classpath;

  my $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 temp_dir          => $temp_dir,
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 1,
                 reference       => $human,
               );
  is ($bam->alignment_metrics_autoqc_command, 'qc --check alignment_filter_metrics --id_run 5175 --position 1 --qc_in t/data/output --qc_out t/data/output --tag_index 1', 'correct alignment_metrics_autoqc_command');

  my $qc_in = getcwd();
  $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => '6062_1#',
                 temp_dir          => $temp_dir,
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 1,
                 reference       => $human,
               );
  is ($bam->alignment_metrics_autoqc_command, "qc --check alignment_filter_metrics --id_run 5175 --position 1 --qc_in $qc_in --qc_out $qc_in --tag_index 1", 'correct alignment_metrics_autoqc_command');

  $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 'output/6062_1#',
                 temp_dir          => $temp_dir,
                 reference         => $human,
                 id_run            => 35,
                 position          => 1,
               );
  is( $bam->alignment_metrics_autoqc_command,
      'qc --check alignment_filter_metrics --id_run 35 --position 1 --qc_in output --qc_out output',
      'alignment metrics autoqc command');
}

{
  local $ENV{CLASSPATH} = $test_classpath;

  my $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 5175,
                 position        => 1,
                 is_paired_read  => 1,
                 tag_index       => 1,
                 non_consent_split => 1,
                 contains_nonconsented_xahuman =>1,
                 no_alignment    => 0,
                 reference       => $human,
               );
  is ($bam->contains_nonconsented_xahuman, 0,
        'if non_consent_split is true contains_nonconsented_xahuman reset to false');
  is ($bam->nonconsented_file_name_suffix, 'human', 'default nonconsented_file_name_suffix');

  $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 1,
                 non_consent_split => 0,
                 contains_nonconsented_xahuman =>1,
                 no_alignment    => 0,
                 reference       => $human,
               );
  is ($bam->contains_nonconsented_xahuman, 1,
        'if non_consent_split is true contains_nonconsented_xahuman is not reset');
  is ($bam->nonconsented_file_name_suffix, 'xahuman', 'nonconsented_file_name_suffix reset to xahuman');
  
  throws_ok {npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 1,
                 non_consent_split => 0,
                 contains_nonconsented_xahuman =>1,
                 no_alignment    => 0,
                 reference       => q[],
               )} qr/contains_nonconsented_xahuman is true, human reference should be used/,
                  'error if reference not given';

  throws_ok {npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 1,
                 non_consent_split => 0,
                 contains_nonconsented_xahuman =>1,
                 no_alignment    => 0,
                 reference       => q[/my/ref/repos/Mouse/strain],
               )} qr/contains_nonconsented_xahuman is true, human reference should be used/,
                  'error if other than human reference is used';

  throws_ok { npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 1,
                 non_consent_split => 0,
                 contains_nonconsented_xahuman =>1,
                 no_alignment    => 1,
                 reference       => $human,
               );} qr/no_alignment and contains_nonconsented_xahuman flags cannot both be true/,
                   'error when both no_alignment and contains_nonconsented_xahuman flags are true';
}

{
  local $ENV{CLASSPATH} = $test_classpath;

  my $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 1,
                 non_consent_split => 1,
                 contains_nonconsented_xahuman =>1,
                 no_alignment    => 0,
                 reference       => $human,
                 java_xmx_flag     => $java_memory,
               );
  my $expected_split_command = qq[$JAVA_CMD $java_memory -jar ] . $current_dir .
           q[/t/bin/aligners/illumina2bam/Illumina2bam-tools-1.00/SplitBamByChromosomes.jar];
  like ($bam->split_by_chr_cmd, qr[$expected_split_command],
      'correct split_by_chromosome command with jar in non-standard location');
}

{
  local $ENV{CLASSPATH} = $test_classpath;
  my $bam = npg_common::sequence::BAM_Alignment->new(
               { input             => 'input.bam',
                 is_paired_read    => 1,
                 non_consent_split => 1,
                 spiked_phix_split => 1,       
                 read_length       => 37,
                 output_prefix     => 'output',
                 bwa_cmd           => 't/bin/aligners/bwa/bwa-0.5.8c/bwa',
                 samtools_cmd      => 't/bin/aligners/samtools/current/samtools',
                 temp_dir          => $temp_dir,
                 reference         => 'reference.fasta',
                 java_xmx_flag     => '-Xmx3000m',
                 id_run            => 35,
                 position          => 1,
               });

  my $illumina2bam_path = abs_path( q[t/bin/aligners/illumina2bam/current] );
  my $expected_alignment_filter_jar = qq[$illumina2bam_path/AlignmentFilter.jar];
  my $expected_filter_cmd = qq{$JAVA_CMD -Xmx3000m -jar $expected_alignment_filter_jar VALIDATION_STRINGENCY=SILENT CREATE_MD5_FILE=true METRICS_FILE=output.bam_alignment_filter_metrics.json IN=input.bam IN=${temp_dir}/human_output_fifo.bam IN=${temp_dir}/output_fifo.bam OUT=output_phix.bam OUT=output_human.bam OUT=output.bam};
  like( $bam->_bam_alignment_filter_command(), qr[$expected_filter_cmd], 'correct filter command with java_xmx_flag explicitly set');
}

{
  local $ENV{CLASSPATH} = $test_classpath;

  throws_ok {npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 1,
                 non_consent_split => 0,
                 contains_nonconsented_xahuman =>1,
                 separate_y_chromosome_data => 1,
                 no_alignment    => 0,
                 reference       => $human,
               )} qr/separate_y_chromosome_data and contains_nonconsented_xahuman flags cannot both be true/,
                  'error if both split by (Y,MT) and split by Y specified';
  
  throws_ok {npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 1,
                 non_consent_split => 0,
                 contains_nonconsented_xahuman =>0,
                 separate_y_chromosome_data => 1,
                 no_alignment    => 0,
                 reference       => q[],
               )} qr/separate_y_chromosome_data is true, human reference should be used/,
                  'error if reference not given';

  throws_ok {npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 1,
                 non_consent_split => 0,
                 contains_nonconsented_xahuman => 0,
                 separate_y_chromosome_data => 1,
                 no_alignment    => 0,
                 reference       => q[/my/ref/repos/Mouse/strain],
               )} qr/separate_y_chromosome_data is true, human reference should be used/,
                  'error if other than human reference is used';

  throws_ok { npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 1,
                 non_consent_split => 0,
                 separate_y_chromosome_data => 1,
                 contains_nonconsented_xahuman =>0,
                 no_alignment    => 1,
                 reference       => $human,
               );} qr/no_alignment and separate_y_chromosome_data flags cannot both be true/,
                  'error when both no_alignment and separate_y_chromosome_data flags are true';
}

{
  local $ENV{CLASSPATH} = $test_classpath;

  throws_ok { npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 5175,
                 is_paired_read  => 1,
                 position        => 1,
                 tag_index       => 1,
                 non_consent_split => 1,
                 contains_nonconsented_xahuman =>1,
                 separate_y_chromosome_data => 1,
                 no_alignment    => 0,
                 reference       => $human,
               );} qr/separate_y_chromosome_data and non_consent_split cannot both be true/,
                  'error when both non_consent_split and separate_y_chromosome_data flags are true';
}

{
  local $ENV{CLASSPATH} = $test_classpath;
  my $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 10371,
                 position        => 1,
                 tag_index       => 33,
                 is_paired_read  => 1,
                 non_consent_split => 0,
                 contains_nonconsented_xahuman => 0,
                 no_alignment    => 0,
                 reference       => $human,
               );

  my $expected_split_command = qq[$JAVA_CMD -Xmx1000m -jar ] . $current_dir . qq[/t/bin/aligners/illumina2bam/Illumina2bam-tools-1.00/SplitBamByChromosomes.jar S=Y V=true];
  like ($bam->split_by_chr_cmd, qr[$expected_split_command], 'correct split_by_chromosome command with no non-consented data for tag 33 in lane 1');
  is ($bam->separate_y_chromosome_data, 1, 'split by Y chromosome is set for tag 33 in lane 1');

  $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 10371,
                 position        => 1,
                 tag_index       => 0,
                 is_paired_read  => 1,
                 non_consent_split => 0,
                 contains_nonconsented_xahuman => 0,
                 no_alignment    => 0,
                 reference       => $human,
               );
  ok ($bam->separate_y_chromosome_data, 'split by Y chromosome is true for tag 0 in lane 1');

  $bam = npg_common::sequence::BAM_Alignment->new(
                 input             => 't/data/sequence/6062_1#0.bam',     
                 output_prefix     => 't/data/output/6062_1#',
                 id_run          => 10371,
                 position        => 1,
                 is_paired_read  => 1,
                 non_consent_split => 0,
                 contains_nonconsented_xahuman => 0,
                 no_alignment    => 0,
                 reference       => $human,
               );
  ok ($bam->separate_y_chromosome_data, 'split by Y chromosome is true for lane 1 with no tag specified');
}

1;
__END__
