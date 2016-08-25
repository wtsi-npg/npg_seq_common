package npg_common::sequence::BAM_MarkDuplicate;

use Moose;
use Carp;
use English qw(-no_match_vars);
use File::Temp qw( tempdir );
use File::Spec::Functions qw(catfile);
use File::Copy qw(move);
use Perl6::Slurp;
use Parallel::ForkManager;
use POSIX qw(mkfifo);
use IPC::Open3;
use File::Basename;
use autodie qw(:all);
use Readonly;

use npg_qc::autoqc::results::bam_flagstats;

with qw/
        MooseX::Getopt
        npg_common::roles::log
        npg_common::roles::software_location
       /;

our $VERSION = '0';

Readonly::Scalar our $FIFO_MODE       => oct '600';
Readonly::Scalar our $EXIT_CODE_SHIFT => 8;

Readonly::Scalar our $STATS1_FILTER       => '0x900';
Readonly::Scalar our $STATS2_FILTER       => '0xB00';

## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::sequence::BAM_MarkDuplicate

=head1 VERSION

=head1 SYNOPSIS

  my $bam = npg_common::sequence::BAM_MarkDuplicate->new(
                 input_bam         => 'input.bam',
                 output_bam        => 'output.bam',
                 metrics_json_dir  => 'qc_dir',
               );
  $bam->process();

=head1 DESCRIPTION

=head1 PUBLIC INTERFACE

=head1 SUBROUTINES/METHODS

=cut

=head2 bammarkduplicates_path

Absolute path to bammarkduplicates executable

=cut

has 'bammarkduplicates_path'  => (
                       is      => 'ro',
                       isa     => 'NpgCommonResolvedPathExecutable',
                       lazy    => 1,
                       builder => '_build_bammarkduplicates_path',
                              );
sub _build_bammarkduplicates_path {
  my $self = shift;
  return catfile($self->biobambam_bin,'bammarkduplicates2');
}

=head2 biobambam_bin

Directory where biobambam family executables are 

=cut 

has 'biobambam_bin' => (
                       isa     => 'NpgTrackingDirectory',
                       is      => 'ro',
                       lazy    => 1,
                       builder => '_build_biobambam_bin',
                            );
sub _build_biobambam_bin {
  my $self = shift;
  return dirname($self->bamsort_cmd());
}

=head2 input_bam

input bam file name with path

=cut

has 'input_bam'      => ( isa             => 'Str',
                          is              => 'rw',
                          required        => 1,
                          documentation   => 'input bam filename with path',
                        );

=head2 output_bam

output bam with path

=cut

has 'output_bam'     => ( isa             => 'Str',
                          is              => 'rw',
                          required        => 1,
                          documentation   => 'output bam filename with path',
                        );

=head2 reference

reference file

=cut

has 'reference'     => ( isa             => 'Str',
                          is              => 'rw',
                          required        => 0,
                          documentation   => 'reference file',
                        );

=head2 metrics_file

output metrics file name

=cut

has 'metrics_file'   => ( isa             => 'Str',
                          is              => 'rw',
                          required        => 0,
                          lazy_build      => 1,
                          documentation   => 'output metrics file name',
                        );
sub _build_metrics_file {
   my $self = shift;
   my $filename = $self->output_bam;
   $filename =~ s/[.]bam$/.markdups_metrics.txt/mxs;
   return $filename;
}

=head2 temp_dir

temp dir

=cut

has 'temp_dir'       => ( isa             => 'Str',
                          is              => 'rw',
                          lazy_build      => 1,
                          documentation   => 'temp dir',
                        );
sub _build_temp_dir {
   my $self = shift;
   my $output = $self->output_bam();
   my ($name,$path,$suffix) = fileparse($output);
   return tempdir(CLEANUP => 1, DIR => $path);
}

=head2 sorted_input_bam_prefix

sorted input bam prefix with path

=cut

has 'sorted_input_bam_prefix'=> ( isa            => 'Str',
                                  is             => 'rw',
                                  lazy_build     => 1,
                                  documentation  => 'sorted input bam prefix with path',
                                );

sub _build_sorted_input_bam_prefix {
   my $self = shift;
   return catfile ($self->temp_dir(), q[sorted]);
}

=head2 metrics_json_dir

output directory for json metrics file

=cut
has 'metrics_json_dir' => ( isa             => 'Str',
                            is              => 'rw',
                            required        => 1,
                            documentation   => 'output directory for json metrics file',
                          );

=head2 id_run

Run id for the lane to be checked.

=cut
has 'id_run'      => (
                       isa            => 'Int',
                       is             => 'rw',
                       required       => 1,
                       documentation  => 'only for metrics json output',
		      );

=head2 position

Lane number. An integer from 1 to 8 inclusive.

=cut
has 'position'    => (isa       => 'Int',
                      is        => 'rw',
                      required  => 1,
                      documentation  => 'only for metrics json output',
                     );

=head2 tag_index

Tag index. An integer from 0 to 10000. Zero means that no tag has been matched.

=cut
has 'tag_index'   => (isa        => 'Int',
                      is         => 'rw',
                      required   => 0,
                      documentation  => 'only for metrics json output',
                     );

=head2 subset

=cut
has 'subset'        => (isa          => 'Str',
                        is           => 'ro',
                        required     => 0,
                        documentation => 'to pass to bam_flagstats autoqc result object',
                       );

=head2 replace_file

replace the input bam file with the output duplicates marked bam file

=cut
has 'replace_file'=> (isa        => 'Bool',
                      is         => 'rw',
                      required   => 0,
                      documentation   => 'replace the input bam file with the output duplicates marked bam file',
                     );

has '_tee_cmd' => ( isa       => 'Str',
                    is        => 'ro',
                    required  => 0,
                    lazy      => 1,
                    builder   => '_build_tee_cmd',
                    documentation => 'entire tee cmd built according to alignment',
                  );

sub _build_tee_cmd {
  my $self = shift;

  my $mark_duplicate_cmd = $self->mark_duplicate_cmd();
  $mark_duplicate_cmd = qq{set -o pipefail;$mark_duplicate_cmd};
  $mark_duplicate_cmd .= ' | tee' ;
  $mark_duplicate_cmd .= q{ } . $self->_bam_md5_fifo_name_mk;
  $mark_duplicate_cmd .= q{ } . $self->_bam_flagstat_fifo_name_mk;
  $mark_duplicate_cmd .= q{ } . $self->_bam_stats1_fifo_name_mk;
  $mark_duplicate_cmd .= q{ } . $self->_bam_stats2_fifo_name_mk;
  $mark_duplicate_cmd .= q{ } . $self->_bam_bschk_fifo_name_mk;
  $mark_duplicate_cmd .= q{ } . $self->_alt_bschk_fifo_name_mk;

  if (!$self->no_alignment) {
    $mark_duplicate_cmd .= q{ } . $self->_bam_index_fifo_name_mk;
    $mark_duplicate_cmd .= q{ } . $self->_bam_pb_cal_fifo_name_mk;
  }
  $mark_duplicate_cmd .= q{ } . $self->_bam_scramble_fifo_name_mk;

  $mark_duplicate_cmd .= ' > ' . $self->output_bam;

  return $mark_duplicate_cmd;
}

=head2 fork_cmds

construct the comands to be run by ForkManager

=cut

sub fork_cmds {
  my $self = shift;

  my @cmds = ();
  my $cmd = q{};

  $cmd = $self->_tee_cmd();
  push @cmds, $cmd;

  my $bamseqchksum_cmd = $self->bamseqchksum_cmd(q{bam});
  $cmd =  'set -o pipefail; cat ' . $self->_bam_bschk_fifo_name_mk . ' | ' . $bamseqchksum_cmd;
  $cmd .=  ' | tee ' . $self->_bam_seqchksum_fifo_name_mk . ' > ' . $self->_bam_seqchksum_file_name_mk;
  push @cmds, $cmd;

  my $altseqchksum_cmd = $self->bamseqchksum_cmd(q{bam});
  $cmd =  'set -o pipefail; cat ' . $self->_alt_bschk_fifo_name_mk . ' | ' . $bamseqchksum_cmd . ' hash=sha512primesums512 > ' . $self->_alt_seqchksum_file_name_mk;
  push @cmds, $cmd;

  $cmd = $self->scramble_cmd() . ' -I bam -O cram';
  if (! $self->no_alignment()) {
    if ($self->reference()) {
      my $refname = $self->reference();
      $refname =~ s{/bwa[\d_\-.]*/}{/fasta/}msx;
      $cmd .= " -r $refname" ;
    }
  }
  $cmd .= ' < ' . $self->_bam_scramble_fifo_name_mk() .  ' | tee ' . $self->_cram_fifo_name_mk ;
  if (! $self->no_alignment()) {
    $cmd .= q{ } . $self->_cram_index_fifo_name_mk ;
  }
  $cmd .= q{ } . $self->_cram_md5_fifo_name_mk . ' > ' . $self->_cram_file_name_mk;
  push @cmds, $cmd;

  $cmd =  'set -o pipefail; cat ' . $self->_cram_fifo_name_mk . ' | ' . $self->bamseqchksum_cmd(q{cram});
  $cmd .=  ' | tee ' . $self->_cram_seqchksum_fifo_name_mk . ' > ' . $self->_cram_seqchksum_file_name_mk;
  push @cmds, $cmd;

  if (! $self->no_alignment()) {
    $cmd =  'set -o pipefail; cat ' . $self->_cram_index_fifo_name_mk(). ' | ' . $self->cram_index_cmd() . q{ - } . $self->_cram_index_file_name_mk();
    push @cmds, $cmd;
  }

  $cmd =  'set -o pipefail; cat ' . $self->_cram_md5_fifo_name_mk(). ' | ' . $self->create_md5_cmd() . q{ > } . $self->_cram_md5_file_name_mk();
  push @cmds, $cmd;

  $cmd = join q{ }, q{diff}, $self->_bam_seqchksum_fifo_name_mk(), $self->_cram_seqchksum_fifo_name_mk();
  push @cmds, $cmd;

  $cmd = 'set -o pipefail; cat ' . $self->_bam_md5_fifo_name_mk . ' | ' . $self->create_md5_cmd() . q{ > } . $self->_md5_file_name_mk;
  push @cmds, $cmd;

  $cmd = 'set -o pipefail; cat ' . $self->_bam_flagstat_fifo_name_mk . ' | ';
  $cmd .= $self->samtools_cmd() . ' flagstat -  > ' . $self->_flagstat_file_name_mk;
  push @cmds, $cmd;

  $cmd = 'set -o pipefail; cat ' . $self->_bam_stats1_fifo_name_mk . ' | ';
  $cmd .= $self->samtools_irods_cmd . ' stats -F ' . $STATS1_FILTER . ' > ' . $self->_stats1_file_name_mk;
  push @cmds, $cmd;

  $cmd = 'set -o pipefail; cat ' . $self->_bam_stats2_fifo_name_mk . ' | ';
  $cmd .= $self->samtools_irods_cmd . ' stats -F ' . $STATS2_FILTER . ' > ' . $self->_stats2_file_name_mk;
  push @cmds, $cmd;

  # we only create an index if we are doing alignment
  if (! $self->no_alignment()) {
    $cmd = 'set -o pipefail; cat ' . $self->_bam_index_fifo_name_mk . ' | ';
    $cmd .= $self->create_index_cmd() . ' > ' . $self->_index_file_name_mk;
    push @cmds, $cmd;

    if ($self->pb_cal_cmd()) {
      my $prefix = $self->output_bam;
      $prefix =~ s/[.]bam$//msx;
      $cmd = 'set -o pipefail; cat ' . $self->_bam_pb_cal_fifo_name_mk . ' | ';
      $cmd .= $self->pb_cal_cmd() . " -p $prefix -filter-bad-tiles 2 -";
      push @cmds, $cmd;
    }
  }

  return \@cmds;
}

=head2 _bam_bschk_fifo_name_mk

Make the bam bschk_fifo name

=cut 

has '_bam_bschk_fifo_name_mk'  => ( isa         => 'Str',
                                    is          => 'ro',
                                    required    => 0,
                                    lazy_build  => 1,
                                  );

sub _build__bam_bschk_fifo_name_mk {
    my $self = shift;

    return $self->output_bam() . q{.bschk.fifo};
}

=head2 _alt_bschk_fifo_name_mk

Make the sha512primesums512 seqchksum generation fifo name

=cut 

has '_alt_bschk_fifo_name_mk'  => ( isa         => 'Str',
                                    is          => 'ro',
                                    required    => 0,
                                    lazy_build  => 1,
                                  );

sub _build__alt_bschk_fifo_name_mk {
    my $self = shift;

    return $self->output_bam() . q{.alt.bschk.fifo};
}

=head2 _bam_md5_fifo_name_mk

Make the bam md5_fifo name

=cut 

has '_bam_md5_fifo_name_mk'  => ( isa         => 'Str',
                                  is          => 'ro',
                                  required    => 0,
                                  lazy_build  => 1,
                                );

sub _build__bam_md5_fifo_name_mk {
    my $self = shift;

    my $bam_md5_fifo_name_mk = $self->output_bam;
    $bam_md5_fifo_name_mk .= q{.md5.fifo};

    return $bam_md5_fifo_name_mk;
}

=head2 _bam_flagstat_fifo_name_mk

Make the bam flagstat_fifo name

=cut 

has '_bam_flagstat_fifo_name_mk'  => (isa         => 'Str',
                                      is          => 'ro',
                                      required    => 0,
                                      lazy_build  => 1,
                                      );

sub _build__bam_flagstat_fifo_name_mk {
    my $self = shift;

    my $bam_flagstat_fifo_name_mk = $self->output_bam;
    $bam_flagstat_fifo_name_mk .= q{.flagstat.fifo};

    return $bam_flagstat_fifo_name_mk;
}

=head2 _bam_index_fifo_name_mk

Make the fifo name for bam index generation

=cut 

has '_bam_index_fifo_name_mk'  => ( isa         => 'Str',
                                    is          => 'ro',
                                    required    => 0,
                                    lazy_build  => 1,
                                  );

sub _build__bam_index_fifo_name_mk {
    my $self = shift;

    return $self->output_bam() . q{.index.fifo};
}

=head2 _cram_index_fifo_name_mk

Make the fifo name for cram index generation

=cut 

has '_cram_index_fifo_name_mk'  => ( isa         => 'Str',
                                    is          => 'ro',
                                    required    => 0,
                                    lazy_build  => 1,
                                  );

sub _build__cram_index_fifo_name_mk {
    my $self = shift;

    return $self->_cram_index_file_name_mk() . q{.fifo};
}

=head2 _bam_stats1_fifo_name_mk

Make the bam stats1_fifo name

=cut 

has '_bam_stats1_fifo_name_mk'  => (isa         => 'Str',
                                      is          => 'ro',
                                      required    => 0,
                                      lazy_build  => 1,
                                      );

sub _build__bam_stats1_fifo_name_mk {
    my $self = shift;

    my $bam_stats1_fifo_name_mk = $self->output_bam;
    $bam_stats1_fifo_name_mk .= q{.stats1.fifo};

    return $bam_stats1_fifo_name_mk;
}

=head2 _bam_stats2_fifo_name_mk

Make the bam stats2_fifo name

=cut 

has '_bam_stats2_fifo_name_mk'  => (isa         => 'Str',
                                      is          => 'ro',
                                      required    => 0,
                                      lazy_build  => 1,
                                      );

sub _build__bam_stats2_fifo_name_mk {
    my $self = shift;

    my $bam_stats2_fifo_name_mk = $self->output_bam;
    $bam_stats2_fifo_name_mk .= q{.stats2.fifo};

    return $bam_stats2_fifo_name_mk;
}

=head2 _bam_pb_cal_fifo_name_mk

Make the bam pb_cal_fifo name

=cut 

has '_bam_pb_cal_fifo_name_mk'  => (  isa         => 'Str',
                                      is          => 'ro',
                                      required    => 0,
                                      lazy_build  => 1,
                                      );

sub _build__bam_pb_cal_fifo_name_mk {
    my $self = shift;

    my $bam_pb_cal_fifo_name_mk = $self->output_bam;
    $bam_pb_cal_fifo_name_mk .= q{.pb_cal.fifo};

    return $bam_pb_cal_fifo_name_mk;
}

=head2 _bam_scramble_fifo_name_mk

Make the bam scramble_fifo name

=cut 

has '_bam_scramble_fifo_name_mk'  => (isa           => 'Str',
                                      is            => 'ro',
                                      required      => 0,
                                      lazy_build    => 1,
                                      );

sub _build__bam_scramble_fifo_name_mk {
    my $self = shift;

    my $bam_scramble_fifo_name_mk = $self->output_bam;
    $bam_scramble_fifo_name_mk .= q{.scramble.fifo};

    return $bam_scramble_fifo_name_mk;
}

=head2 _cram_file_name_mk

Make the cram file name

=cut 

has '_cram_file_name_mk'  => (isa         => 'Str',
                              is          => 'ro',
                              required    => 0,
                              lazy_build  => 1,
                              );

sub _build__cram_file_name_mk {
    my $self = shift;

    my $cram_file_name_mk = $self->output_bam;
    $cram_file_name_mk =~ s/[.]bam$/.cram/mxs;

    return $cram_file_name_mk;
}

=head2 _cram_md5_file_name_mk

Make the cram md5 file name

=cut 

has '_cram_md5_file_name_mk'  => (isa         => 'Str',
                              is          => 'ro',
                              required    => 0,
                              lazy_build  => 1,
                              );

sub _build__cram_md5_file_name_mk {
    my $self = shift;

    return $self->_cram_file_name_mk() . q{.md5};
}

=head2 _cram_md5_fifo_name_mk

Make the fifo name for cram md5 generation

=cut 

has '_cram_md5_fifo_name_mk'  => (isa         => 'Str',
                              is          => 'ro',
                              required    => 0,
                              lazy_build  => 1,
                              );

sub _build__cram_md5_fifo_name_mk {
    my $self = shift;

    return $self->_cram_md5_file_name_mk() . q{.fifo};
}

=head2 _cram_index_file_name_mk

Make the cram index file name

=cut 

has '_cram_index_file_name_mk'  => (isa         => 'Str',
                              is          => 'ro',
                              required    => 0,
                              lazy_build  => 1,
                              );

sub _build__cram_index_file_name_mk {
    my $self = shift;

    return $self->_cram_file_name_mk() . q{.crai};
}

=head2 _cram_fifo_name_mk

Make the cram fifo name

=cut 

has '_cram_fifo_name_mk'  => (isa         => 'Str',
                              is          => 'ro',
                              required    => 0,
                              lazy_build  => 1,
                              );

sub _build__cram_fifo_name_mk {
    my $self = shift;

    my $cram_fifo_name_mk = $self->output_bam;
    $cram_fifo_name_mk =~ s/[.]bam$/.cram.fifo/mxs;

    return $cram_fifo_name_mk;
}

=head2 _cram_seqchksum_fifo_name_mk

Make the cram_seqchksum fifo name

=cut 

has '_cram_seqchksum_fifo_name_mk'  => (isa           => 'Str',
                                        is            => 'ro',
                                        required      => 0,
                                        lazy_build    => 1,
                                        );

sub _build__cram_seqchksum_fifo_name_mk {
    my $self = shift;

    my $cram_seqchksum_fifo_name_mk = $self->output_bam;
    $cram_seqchksum_fifo_name_mk =~ s/[.]bam$/.cram.seqchksum.fifo/mxs;

    return $cram_seqchksum_fifo_name_mk;
}

=head2 _cram_seqchksum_file_name_mk

Make the cram_seqchksum file name

=cut 

has '_cram_seqchksum_file_name_mk'  => (isa         => 'Str',
                                        is          => 'ro',
                                        required    => 0,
                                        lazy_build  => 1,
                                        );

sub _build__cram_seqchksum_file_name_mk {
    my $self = shift;

    my $cram_seqchksum_file_name_mk = $self->output_bam;
    $cram_seqchksum_file_name_mk =~ s/[.]bam$/.cram.seqchksum/mxs;

    return $cram_seqchksum_file_name_mk;
}

=head2 _bam_seqchksum_file_name_mk

Make the bam_seqchksum file name

=cut 

has '_bam_seqchksum_file_name_mk'  => ( isa         => 'Str',
                                        is          => 'ro',
                                        required    => 0,
                                        lazy_build  => 1,
                                        );

sub _build__bam_seqchksum_file_name_mk {
    my $self = shift;

    return $self->output_bam . q{.seqchksum};
}

=head2 _alt_seqchksum_file_name_mk

Make the alt_seqchksum file name

=cut 

has '_alt_seqchksum_file_name_mk'  => ( isa         => 'Str',
                                        is          => 'ro',
                                        required    => 0,
                                        lazy_build  => 1,
                                        );

sub _build__alt_seqchksum_file_name_mk {
    my $self = shift;

    return $self->output_bam . q{.sha512primesums512.seqchksum};
}

=head2 _bam_seqchksum_fifo_name_mk

Make the bam_seqchksum fifo name

=cut 

has '_bam_seqchksum_fifo_name_mk'  => ( isa         => 'Str',
                                        is          => 'ro',
                                        required    => 0,
                                        lazy_build  => 1,
                                      );

sub _build__bam_seqchksum_fifo_name_mk {
    my $self = shift;

    my $bam_seqchksum_fifo_name_mk = $self->output_bam;
    $bam_seqchksum_fifo_name_mk .= q{.seqchksum.fifo};

    return $bam_seqchksum_fifo_name_mk;
}

=head2 _index_file_name_mk

Make the index file name

=cut 

has '_index_file_name_mk'  => (isa        => 'Str',
                               is         => 'ro',
                               required   => 0,
                               lazy_build => 1,
                              );

sub _build__index_file_name_mk {
    my $self = shift;

    my $index_file_name_mk = $self->output_bam;
    $index_file_name_mk =~ s/[.]bam$/.bai/mxs;

    return $index_file_name_mk;
}

=head2 _flagstat_file_name_mk

Make the flagstat file name

=cut 

has '_flagstat_file_name_mk'  => (isa         => 'Str',
                                  is          => 'ro',
                                  required    => 0,
                                  lazy_build  => 1,
                                  );

sub _build__flagstat_file_name_mk {
    my $self = shift;

    my $flagstat_file_name_mk = $self->output_bam;
    $flagstat_file_name_mk =~ s/[.]bam$/.flagstat/mxs;

    return $flagstat_file_name_mk;
}

=head2 _stats1_file_name_mk

Make the stats1 file name

=cut 

has '_stats1_file_name_mk'  => (isa         => 'Str',
                                is          => 'ro',
                                required    => 0,
                                lazy_build  => 1,
                               );

sub _build__stats1_file_name_mk {
    my $self = shift;

    my $stats1_file_name_mk = $self->output_bam;
    $stats1_file_name_mk =~ s/[.]bam$/_F$STATS1_FILTER.stats/mxs;

    return $stats1_file_name_mk;
}

=head2 _stats2_file_name_mk

Make the stats2 file name

=cut 

has '_stats2_file_name_mk'  => (isa         => 'Str',
                                is          => 'ro',
                                required    => 0,
                                lazy_build  => 2,
                               );

sub _build__stats2_file_name_mk {
    my $self = shift;

    my $stats2_file_name_mk = $self->output_bam;
    $stats2_file_name_mk =~ s/[.]bam$/_F$STATS2_FILTER.stats/mxs;

    return $stats2_file_name_mk;
}

=head2 _md5_file_name_mk

Make the md5 file name

=cut 

has '_md5_file_name_mk'  => (isa        => 'Str',
                             is         => 'ro',
                             required   => 0,
                             lazy_build => 1,
                             );

sub _build__md5_file_name_mk {
    my $self = shift;

    my $md5_file_name_mk = $self->output_bam;
    $md5_file_name_mk =~ s/[.]bam$/.bam.md5/mxs;

    return $md5_file_name_mk;
}

=head2 sort_cmd

construct the whole biobambam sort command

=cut

sub sort_cmd {
   my $self = shift;

   my $cmd = $self->bamsort_cmd();
   $cmd .= q{ tmpfile=}.$self->temp_dir().q{/};
   $cmd .= q{ < } . $self->input_bam();
   $cmd .= q{ > } . $self->sorted_input_bam_prefix() . q {.bam};

   return $cmd;
}

=head2 mark_duplicate_cmd

construct the whole biobambam markduplicate command

=cut

sub mark_duplicate_cmd {
   my $self = shift;

   my $cmd = $self->bammarkduplicates_path;
   if( !$self->no_alignment() ){
      $cmd .= q{ I=} . $self->sorted_input_bam_prefix() . q{.bam};
   }else {
      $cmd .= q{ I=}.$self->input_bam();
   }
   $cmd .= q{ O=/dev/stdout};
   $cmd .= q{ tmpfile=}.$self->temp_dir().q{/};
   $cmd .= q{ M=}.$self->metrics_file();

   return $cmd;
}


has 'create_index_cmd' => ( isa           => 'Str',
                            is            => 'ro',
                            lazy_build    => 1,
                            documentation => 'Command to create BAM index',
                          );

sub _build_create_index_cmd {
   my $self = shift;
   return  $self->samtools_cmd .  q[ index /dev/stdin /dev/stdout];
}

has 'create_md5_cmd' => ( isa           => 'Str',
                          is            => 'ro',
                          default       => 'md5sum -b | tr -d "\\n *-" ', #note \\ squashes down to \ even in this 'unevaluated' string
                          documentation =>'Command to create MD5 checksum',
                        );

=head2 bamseqchksum_cmd

a command to execute bamseqchksum for a bam or cram file to create file_name.[b|cr]am.seqchksum 

=cut

sub bamseqchksum_cmd {
  my $self = shift;
  my $file_type = shift;

  my $chk_command = $self->bamseqchksum().qq( verbose=0 inputformat=$file_type);
  my $output = q{};

  if ($file_type eq q(cram) ) {
    if ($self->reference()) {
      my $reference = $self->reference();
      $reference =~ s{/bwa[\d_\-.]*/}{/fasta/}msx;
      $chk_command .= qq{ reference=$reference};
    }
  }

  return $chk_command;
}

=head2 bamsort_cmd

biobambam sort command for the input bam

=cut
has 'bamsort_cmd'   => ( is      => 'ro',
                         isa     => 'NpgCommonResolvedPathExecutable',
                         coerce  => 1,
                         default => 'bamsort',
                       );


=head2 bamcheck_flags

bamcheck flags for bamcheck command

=cut
has 'bamcheck_flags'   => ( is      => 'ro',
                         isa     => 'Str',
                         default => q{},
                       );


=head2 scramble_cmd

scramble command for the input bam

=cut
has 'scramble_cmd'   => ( is      => 'ro',
                         isa     => 'NpgCommonResolvedPathExecutable',
                         coerce  => 1,
                         default => 'scramble',
                       );


=head2 cram_index_cmd

cram_index command for the input bam

=cut
has 'cram_index_cmd'   => ( is      => 'ro',
                         isa     => 'NpgCommonResolvedPathExecutable',
                         coerce  => 1,
                         default => 'cram_index',
                       );


=head2 calibration command

calibration command for the input bam

=cut
has 'pb_cal_cmd'   => ( is      => 'ro',
                         isa     => 'NpgCommonResolvedPathExecutable',
                         coerce  => 1,
                         default => 'calibration_pu',
                       );


=head2 bamseqchksum command

biobambam bamseqchksum command. Default checksum uses hash=crc32prod. Other checksum to use is sha512primesums512

 
=cut
has 'bamseqchksum' => ( is  => 'ro',
                            isa     => 'NpgCommonResolvedPathExecutable',
                            coerce  => 1,
                            default => 'bamseqchksum',
                          );


=head2 no_alignment

to indicate input bam with alignment or not

=cut
has 'no_alignment'=> (isa             => 'Bool',
                      is              => 'rw',
                      required        => 0,
                      lazy_build      => 1,
                      documentation   => 'indicate input bam with alignment or not',
                     );
sub _build_no_alignment {
   my $self = shift;

   my $no_alignment = 1;

   my $samtools_view_cmd = $self->samtools_cmd() . q{ view -H } . $self->input_bam();

   $self->log($samtools_view_cmd);

   my $pid = open3( undef, my $input_header_fh, undef, $samtools_view_cmd );

   while( my $line = <$input_header_fh>){
      if($line =~ /^\@SQ\t/mxs){
         $no_alignment = 0;
         $self->log('SQ tag found in input bam file, the input bam with alignment: ' . $self->input_bam());
         last;
      }
   }

   close $input_header_fh;

   if( $no_alignment ){
      $self->log('NO SQ tag found in input bam file, the input bam without alignment:' . $self->input_bam());
   }

   return $no_alignment;
}

=head2 process

main method to call

=cut

sub process {
  my $self = shift;

  if( ! -e $self->input_bam() ){
     croak 'Input bam file does not exist to mark duplicate: '.$self->input_bam();
  }

  if( !$self->no_alignment() ){
    $self->log('Sort input bam file');
    my $sort_cmd = $self->sort_cmd();
    $self->log("sort_command: $sort_cmd");
    my $sort_rs = system $sort_cmd;
    if( $sort_rs ) {
       croak 'Biobambam sort failed!';
    }
  }

  my @fifos = ();
  push @fifos,  $self->_bam_md5_fifo_name_mk;
  push @fifos,  $self->_bam_flagstat_fifo_name_mk;
  push @fifos,  $self->_bam_stats1_fifo_name_mk;
  push @fifos,  $self->_bam_stats2_fifo_name_mk;
  push @fifos,  $self->_bam_bschk_fifo_name_mk;
  push @fifos,  $self->_alt_bschk_fifo_name_mk;

  if (! $self->no_alignment()) {
    push @fifos,  $self->_bam_pb_cal_fifo_name_mk;
    push @fifos,  $self->_bam_index_fifo_name_mk;
    push @fifos,  $self->_cram_index_fifo_name_mk;
  }

  push @fifos,  $self->_bam_seqchksum_fifo_name_mk;
  push @fifos,  $self->_bam_scramble_fifo_name_mk;
  push @fifos,  $self->_cram_seqchksum_fifo_name_mk;
  push @fifos,  $self->_cram_md5_fifo_name_mk;
  push @fifos,  $self->_cram_fifo_name_mk;

  foreach my $fifo (@fifos) {
    $self->log("Making fifo $fifo");
    mkfifo $fifo, $FIFO_MODE or croak "Failed to make fifo $fifo ($ERRNO)";
  }

  $self->log('Setting up fork commands');

  my $fork_cmds = $self->fork_cmds();

  my $pm = Parallel::ForkManager->new(scalar @{$fork_cmds});
  setpgrp;

  $pm->run_on_finish(
    sub { my ($pid, $exit_code, $ident) = @_;
      if ($exit_code){
        carp "PID $pid returned exit code: $exit_code. Fail: $ident";
        carp 'Killing all sibling processes';
        local $SIG{TERM} = 'IGNORE';
        kill TERM => 0;
        croak 'Killed all sibling processes';
      } else {
        $self->log( "PID $pid and exit code: $exit_code. Success: $ident") ;
      }
    }
  );

  $pm->run_on_start(
    sub { my ($pid,$ident)=@_;
      $self->log( "Job $pid started: $ident" );
    }
  );

  foreach my $command (@{$fork_cmds}) {
    $pm->start($command) and next;
    exec $command or
      $self->log("$command does not exist"); # only checks the first command in the pipe
    $pm->finish($CHILD_ERROR >> $EXIT_CODE_SHIFT);
  }

  $pm->wait_all_children;

  if ( !$self->_tee_cmd() ) {
    $self->log('WEIRDNESS WARNING: mark_duplicate_cmd is NOT returned');
  }
  $self->_finalise_output(\@fifos);
  $self->_bam_flagstats_qc();

  $self->log('Finished in BAM_Markduplicate!');

  return 1;
}

sub _finalise_output {
  my ($self, $fifos) = @_;
  $self->log('Replacing input file with duplicates marked file');

  $self->_move_file($self->output_bam, $self->input_bam);

  my $to_move = {
      $self->metrics_file                 => '.markdups_metrics.txt',
      $self->_flagstat_file_name_mk       => '.flagstat',
      $self->_index_file_name_mk          => '.bai',
      $self->_stats1_file_name_mk         => '_F' . $STATS1_FILTER . '.stats',
      $self->_stats2_file_name_mk         => '_F' . $STATS2_FILTER . '.stats',
      $self->_md5_file_name_mk            => '.bam.md5',
      $self->_alt_seqchksum_file_name_mk  => '.sha512primesums512.seqchksum',
      $self->_bam_seqchksum_file_name_mk  => '.seqchksum',
      $self->_cram_file_name_mk           => '.cram',
      $self->_cram_index_file_name_mk     => '.cram.crai',
      $self->_cram_md5_file_name_mk       => '.cram.md5',
  };
  foreach my $suf (qw(_quality_cycle_caltable.txt _quality_cycle_surv.txt _quality_error.txt)) {
    my $file_name_mk = $self->output_bam;
    $file_name_mk =~ s/[.]bam$/$suf/mxs;
    $to_move->{$file_name_mk} = $suf;
  }

  my $root = $self->input_bam;
  $root =~ s/[.]bam\Z//xms;
  while ( my ($f, $ext) = each %{$to_move} ) {
    $self->_move_file($f, $root.$ext);
  }

  if (-e $self->_cram_seqchksum_file_name_mk) {
    $self->_remove_file($self->_cram_seqchksum_file_name_mk);
  }

  $self->log('Removing fifos');
  foreach my $fifo ( @{$fifos} ) {
    $self->_remove_file($fifo);
  }

  return;
}

sub _remove_file {
  my ($self, $source) = @_;
  $self->log("Deleting $source");
  unlink $source or croak "Failed to delete $source ($ERRNO)";
  return;
}

sub _move_file {
  my ($self, $source, $des) = @_;
  if (-e $source) {
    $self->log("Moving $source to $des");
    move $source, $des or croak "Failed to move $source to $des ($ERRNO)";
  }
  return;
}

sub _bam_flagstats_qc {
  my $self = shift;

  my $result = {};
  $result->{'sequence_file'} = $self->input_bam;
  $result->{'id_run'}        = $self->id_run;
  $result->{'position'}      = $self->position;
  if(defined $self->tag_index){
    $result->{'tag_index'}   = $self->tag_index;
  }
  if($self->subset){
    $result->{'subset'}      = $self->subset;
  }

  my $robj = npg_qc::autoqc::results::bam_flagstats->new($result);
  $self->log('Calling bam_flagstats execute() method');
  $robj->execute();
  $robj->set_info('Samtools', $self->current_version($self->samtools_cmd) || q[not known]);
  $self->log('Serializing bam_flagstats object to ' . $self->metrics_json_dir);
  $robj->store($self->metrics_json_dir);

  return $robj;
}


no Moose;

1;
__END__

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose

=item MooseX::Getopt

=item Carp

=item English -no_match_vars

=item Readonly

=item Perl6::Slurp

=item File::Temp

=item File::Basename

=item File::Spec::Functions

=item File::Copy

=item autodie qw(:all)

=item IPC::Open3

=item POSIX

=item npg_qc::autoqc::results::bam_flagstats

=item npg_common::roles::log

=item npg_common::roles::software_location

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Guoying Qi E<lt>gq1@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015 GRL

This file is part of NPG.

NPG is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
