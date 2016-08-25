package npg_common::sequence::BAM_Alignment;

use Carp;
use autodie qw(:all);
use English qw(-no_match_vars);
use Moose;
use MooseX::StrictConstructor;
use File::Spec::Functions qw(catfile splitdir catdir);
use File::Temp qw(tempfile tempdir);
use File::Spec;
use File::Path qw(remove_tree);
use Cwd qw(cwd);
use File::Basename;
use Parallel::ForkManager;
use Perl6::Slurp;
use FindBin qw($Bin);
use File::Which;

use st::api::lims;
use npg_tracking::data::reference;
use npg_tracking::data::reference::info;
use npg_tracking::util::abs_path qw(abs_path);

with qw/
        MooseX::Getopt
        npg_common::roles::log
        npg_common::roles::software_location
        npg_tracking::glossary::lane
        npg_tracking::glossary::run
       /;
with 'npg_tracking::illumina::run::folder';
with 'npg_tracking::illumina::run::long_info';

use Readonly;

our $VERSION = '0';

Readonly::Scalar my $PICARD_SAM_FORMAT_CONVERTER_JAR
                                               => q[SamFormatConverter.jar];
Readonly::Scalar my $ILLUMINA2BAM_BAM_MERGER_JAR
                                               => q[BamMerger.jar];
Readonly::Scalar my $ILLUMINA2BAM_ALIGNMENT_FILTER_JAR
                                               => q[AlignmentFilter.jar];
Readonly::Scalar my $ILLUMINA2BAM_SPLIT_BAM_BY_CHROMOSOMES_JAR
                                               => q[SplitBamByChromosomes.jar];
Readonly::Scalar my $QC_EXECUTABLE_NAME        => q[qc];

Readonly::Scalar my $BAM_MARK_DUPLICATES_CMD  => q[bam_mark_duplicate.pl];

Readonly::Scalar my $DEFAULT_JAVA_XMX         => q{-Xmx1000m};
Readonly::Scalar my $JAVA_XMX_DUPLICATES      => 6000;

Readonly::Scalar my $DEFAULT_BWA_TRIM_QUALITY => q{-q 15};
Readonly::Scalar my $DEFAULT_BWA_THREADS      => 4;
Readonly::Scalar my $DEFAULT_BWA_SAM_MAX_THREADS      => 6;
Readonly::Scalar my $MIN_READ_LENGTH_TO_ALIGN => 31;

Readonly::Scalar my $EXIT_CODE_SHIFT          => 8;

## no critic (Documentation::RequirePodAtEnd)

=head1 NAME

npg_common::sequence::BAM_Alignment

=head1 VERSION

=head1 SYNOPSIS

  my $bam = npg_common::sequence::BAM_Alignment->new({
                 input             => 'input.bam',
                 is_paired_read    => 1,
                 non_consent_split => 1,
                 spiked_phix_split => 1,
               });
  $bam->generate();

=head1 DESCRIPTION

=head1 PUBLIC INTERFACE

=head1 SUBROUTINES/METHODS

=cut

=head2 input

input bam file, probably with phix alignment

=cut
has 'input'           => (isa           => q{Str},
                          is            => q{ro},
                          required      => 1,
                          documentation => q{input bam file, probably with phix alignment},
                         );
around 'is_paired_read' => sub {
  my $orig = shift;
  my $self = shift;

  if (!$self->has_is_paired_read) {
    return;
  }
  return $self->$orig();
};

=head2 id_run

=head2 position

=head2 tag_index 	 
	  	 
optional field tag_index to get lims information
	  	 
=cut

has 'tag_index'  => (isa           => q{Int},
	             is            => q{rw},
	             required      => 0,
	             documentation => q{tag index},
	            );

=head2 _lims

lims access object

=cut
has  '_lims'    =>  ( isa             => 'Maybe[st::api::lims]',
                      is              => 'ro',
                      init_arg        => undef,
                      lazy_build      => 1,
                    );
sub _build__lims {
    my $self = shift;
    my $lims;
    eval{
        $lims = st::api::lims->new(id_run    => $self->id_run,
                                   position  => $self->position,
                                   tag_index => $self->tag_index,
                                  );
        # Creating an object does not necessary mean that LIMs
        # data are available; getting the data reveals problems.
        $lims->study_id();
        1;
    } or do {
        carp "Failed to create st::api::lims object\n$EVAL_ERROR";
        $lims = undef;
    };

    return $lims;
}

=head2 spiked_phix_split

indicate the input file contains spiked phix alignment, which are to be split out

=cut
has 'spiked_phix_split'=> (isa           => q{Bool},
                           is            => q{ro},
                           required      => 0,
                           lazy_build    => 1,
                           documentation => q{indicate the input file contains spiked phix alignment to be split out},
                          );
sub _build_spiked_phix_split {
  my $self = shift;

  my $spiked_phix_tag_index = $self->_lims() ? $self->_lims()->spiked_phix_tag_index() : 0;
  if( $spiked_phix_tag_index ){
    if( defined $self->tag_index() && $self->tag_index() == $spiked_phix_tag_index ){
      $self->log('The input is spiked phix');
      return 0;
    }
    return 1;
  }
  return 0;
}

=head2 separate_y_chromosome_data

indicate the input file is to be split into Y and non Y 

=cut
has 'separate_y_chromosome_data'=> (isa           => q{Bool},
                           is            => q{ro},
                           required      => 0,
                           lazy_build    => 1,
                           writer        => '_set_separate_y_chromosome_data',
                           documentation => q{indicate the input file is to be split into Y and non Y},
                          );
sub _build_separate_y_chromosome_data{
  my $self = shift;

  return $self->_lims() ? $self->_lims()->separate_y_chromosome_data() : 0;
}

=head2 non_consented_split

indicate the input file contains non consented human data or not

=cut
has 'non_consent_split'=> (isa           => q{Bool},
                           is            => q{ro},
                           required      => 0,
                           lazy_build    => 1,
                           writer        => '_set_non_consent_split',
                           documentation => q{indicate the input file contains non consented human data or not},
                          );
sub _build_non_consent_split{
  my $self = shift;
  return $self->_lims() ? $self->_lims()->contains_nonconsented_human() : 0;
}


=head2 contains_nonconsented_xahuman

indicates whether the input contains nonconcented X and autosomal human

=cut
has 'contains_nonconsented_xahuman' => (isa           => q{Bool},
                                        is            => q{ro},
                                        required      => 0,
                                        lazy_build    => 1,
                                        writer        => '_set_contains_nonconsented_xahuman',
                                        documentation => q{indicates whether the input contains nonconcented X and autosomal human},
                                       );
sub _build_contains_nonconsented_xahuman {
   my $self = shift;
   return $self->_lims() ? $self->_lims()->contains_nonconsented_xahuman() : 0;
}

=head2 human_reference

A path to human reference for bwa alignment to split non-consented reads

=cut
has 'human_reference'   => (isa           => 'Str',
                            is            => 'rw',
                            required      => 0,
                            lazy_build    => 1,
                            documentation => 'A path to human reference for bwa alignment to split non-consented reads',
                           );
sub _build_human_reference {
   my $self = shift;

   my $ruser = Moose::Meta::Class->create_anon_class(
          roles => [qw/npg_tracking::data::reference::find/])
          ->new_object({species => q{Homo_sapiens}});

   my $human_ref;
   eval {
      $human_ref = $ruser->refs->[0];
      1;
   } or do{
      carp $EVAL_ERROR;
   };

   return $human_ref;
}

=head2 phix_reference

A path to phix reference

=cut
has 'phix_reference'   => (isa           => 'Str',
                            is            => 'rw',
                            required      => 0,
                            lazy_build    => 1,
                            documentation => 'A path to phix reference',
                           );
sub _build_phix_reference {
   my $self = shift;

   my $ruser = Moose::Meta::Class->create_anon_class(
          roles => [qw/npg_tracking::data::reference::find/])
          ->new_object({species => q{PhiX}});

   my $phix_ref;
   eval {
      $phix_ref = $ruser->refs->[0];
      1;
   } or do{
      carp $EVAL_ERROR;
   };

   return $phix_ref;
}

=head2 human_ref_dict

A path to human reference picard dictionary file for non-consented bam output

=cut
has 'human_ref_dict'    => (isa           => 'Str',
                            is            => 'rw',
                            required      => 0,
                            lazy_build    => 1,
                            documentation => 'A path to human reference picard dictionary file for non-consented bam output',
                           );
sub _build_human_ref_dict {
   my $self = shift;

   my $ruser = Moose::Meta::Class->create_anon_class(
          roles => [qw/npg_tracking::data::reference::find/])
          ->new_object({species => q{Homo_sapiens}, aligner => q{picard} } );

   my $human_ref;
   eval {
      $human_ref = $ruser->refs->[0];
      $human_ref .= q{.dict};
      1;
   } or do{
      carp $EVAL_ERROR;
   };

   return $human_ref;
}
=head2 read_length

read length in input file

=cut
has 'read_length'     => (isa           => q{Int},
                          is            => q{rw},
                          required      => 0,
                          documentation => q{read lenght in input file},
                         );
#TODO lazy build based on input bam file


=head2 reference_info

an object of reference_info including bwa options and reference prefix for BWA alignment

=cut

has 'reference_info'  => (isa           => 'Maybe[npg_tracking::data::reference::info]',
                          is            => 'rw',
                          required      => 0,
                          lazy_build    => 1,
                          documentation => 'Reference info object including bwa options and reference prefix for BWA alignment',
                         );

sub _build_reference_info {
  my $self = shift;

  if (!$self->_lims()) { return; }

  my $ref_bwa = npg_tracking::data::reference->new(
             lims      => $self->_lims(),
             id_run    => $self->id_run(),
             position  => $self->position(),
             tag_index => $self->tag_index(),
            );

  my $ref_info;
  eval {
    $ref_info = $ref_bwa->ref_info();
    1;
  } or do {
    carp qq{Error: reference bwa options cannot be retrieved; $EVAL_ERROR};
  };

  return $ref_info;
}


=head2 reference

A path to reference for bwa alignment

=cut
has 'reference'   => (isa           => 'Str',
                      is            => 'rw',
                      required      => 0,
                      lazy_build    => 1,
                      documentation => 'A path to reference for bwa alignment',
                     );
sub _build_reference {
  my $self = shift;

  my $ref_info = $self->reference_info();
  if(!$ref_info){
     return q{};
  }
  my($filename, $directories, $suffix) = fileparse($ref_info->ref_path());
  return catfile(abs_path($directories), $filename);
}

=head2 ref_dict

A path to reference picard dictionary for bam output

=cut
has 'ref_dict'    => (isa           => 'Str',
                      is            => 'rw',
                      required      => 0,
                      lazy_build    => 1,
                      documentation => 'A path to reference picard dictionary for bam output',
                     );

sub _build_ref_dict {
  my $self = shift;

  if (!$self->_lims()) { return; }

  my $ref_picard = npg_tracking::data::reference->new(
             lims      => $self->_lims,
             id_run => $self->id_run(),
             position  => $self->position(),
             aligner   => q{picard},
             tag_index => $self->tag_index(),
            );
  my $refs;
  eval {
    $refs = $ref_picard->refs();
    1;
  } or do {
    my $error_message = qq{Error: binary reference fasta cannot be retrieved; $EVAL_ERROR};
    carp $error_message;
    return q{};
  };

  if (scalar @{$refs} > 1) {
    carp q[Run ] . $self->id_run . q[ lane ] . $self->position . q[ - multiple references for a pooled lane: ] . (join q[;], @{$refs});
    return q{};
  }

  if(! scalar @{$refs}){
    return q{};
  }

  my $reference_dict = pop @{$refs};
  $reference_dict .= q{.dict};

  if(! -e $reference_dict){
    croak q[Run ] . $self->id_run . q[ lane ] . $self->position . q[ - picard dict missing; ].$reference_dict;
  }

  $self->log('Reference Picard dictionary found: ' . $reference_dict );
  return abs_path( $reference_dict );
}

=head2 no_alignment

skip the alignment, generate bam file without alignment

=cut
has 'no_alignment' => (isa           => q{Bool},
                       is            => q{ro},
                       required      => 0,
                       lazy_build    => 1,
                       documentation => q{skip alignment},
                      );

sub _build_no_alignment {
  my $self = shift;

  if( $self->read_length() && $self->read_length() < $MIN_READ_LENGTH_TO_ALIGN ){
    $self->log('Read length too short for alignment');
    #TODO don't do this checking here, add option to bwa aln to change seed length to do alignment
    return 1;
  }elsif( !$self->reference() ){
    $self->log('No reference given for alignment');
    return 1;
  }
  if ($self->_lims() && $self->_lims()->study_id) {
    return !$self->_lims()->alignments_in_bam;
  }
  return 0; # this means we will try to align at lane level for plexes and tag zero
            # regardless of the options set for individual tags
}


=head2 output_prefix

prefix for output bam file

=cut
has 'output_prefix'   => (isa           => q{Str},
                          is            => q{ro},
                          required      => 1,
                          documentation => q{prefix for output bam file},
                         );

=head2 spiked_phix_output

output bam file name with path for spiked phix

=cut
has 'spiked_phix_output'=> (isa          => q{Str},
                           is            => q{ro},
                           required      => 0,
                           lazy_build    => 1,
                           documentation => q{output bam file name with path for spiked phix},
                          );
sub _build_spiked_phix_output{
   my $self = shift;
   return $self->output_prefix() . q{_phix.bam};
}

=head2 nonconsented_file_name_suffix

suffix to be appended to file names for nonconsented output

=cut
has 'nonconsented_file_name_suffix' => (isa           => q{Str},
                                        is            => q{rw},
                                        required      => 0,
                                        default       => 'human',
                                        documentation => q{suffix to be appended to file names for nonconsented output},
                                       );

=head2 non_consented_output

output bam file name with path for non consented data

=cut
has 'non_consented_output' => (isa           => q{Str},
                               is            => q{ro},
                               required      => 0,
                               lazy_build    => 1,
                               documentation => q{output bam file name with path for non consented data},
                              );
sub _build_non_consented_output{
   my $self = shift;
   return $self->output_prefix() . q{_} . $self->nonconsented_file_name_suffix . q{.bam};
}

=head2 output

output bam file name with path for target data

=cut
has 'output'               => (isa           => q{Str},
                               is            => q{ro},
                               required      => 0,
                               lazy_build    => 1,
                               documentation => q{output bam file name with path for target data},
                              );
sub _build_output{
   my $self = shift;
   return $self->output_prefix() . q{.bam};
}

=head2 bamcheck_flags

bamcheck flags for bamcheck command - simply passed through to bam_mark_duplicate.pl, which handles any validation

=cut
has 'bamcheck_flags'   => ( is      => 'ro',
                         isa     => 'Str',
                         default => q{},
                       );

=head2 _bam_output_fifo

output bam fifo name from bwa

=cut
has '_bam_output_fifo'     => (isa           => 'Str',
                               is            => 'ro',
                               required      => 0,
                               lazy_build    => 1,
    );

sub _build__bam_output_fifo{
   my $self = shift;
   return $self->temp_dir(). q{/output_fifo.bam};
}

=head2 _human_bam_output_fifo

output bam fifo name from bwa against human reference

=cut
has '_human_bam_output_fifo' => (isa           => 'Str',
                                 is            => 'ro',
                                 required      => 0,
                                 lazy_build    => 1,
    );

sub _build__human_bam_output_fifo{
   my $self = shift;
   return $self->temp_dir(). q{/human_output_fifo.bam};
}

=head2 bwa_aln_options

bwa aln command options, for eaxmple, '-q 15 -t 6'
 
=cut
has 'bwa_aln_options' => (isa             => 'Str',
                          is              => 'rw',
                          required        => 0,
                          lazy_build      => 1,
                          documentation   => 'aligner bwa aln command options, for eaxmple, -q 15 -t 6',
                         );
sub _build_bwa_aln_options {
  my $self = shift;

  my $ref_info = $self->reference_info();
  my $aligner_options = $ref_info ? $ref_info->aligner_options() : $DEFAULT_BWA_TRIM_QUALITY;
  if( $aligner_options eq $npg_tracking::data::reference::find::NPG_DEFAULT_ALIGNER_OPTION ){
    $aligner_options = $DEFAULT_BWA_TRIM_QUALITY;
  }
  return $aligner_options.q{ -t } . $self->bwa_aln_threads();
}

=head2 bwa_aln_threads

The number of threads for bwa aln command
 
=cut
has 'bwa_aln_threads' => (isa             => 'Int',
                          is              => 'ro',
                          lazy_build      => 1,
                          documentation   => 'The number of threads for bwa aln command',
                         );
sub _build_bwa_aln_threads {
    if (my$mnslot = $ENV{LSB_MCPU_HOSTS}) { return $1 if $mnslot=~/\S+\s+(\d+)/smx; }
    return $DEFAULT_BWA_THREADS;
}

=head2 bwa_sam_threads

The number of threads for bwa sam{se,pe} command
 
=cut
has 'bwa_sam_threads' => (isa             => 'Int',
                          is              => 'ro',
                          lazy_build      => 1,
                          documentation   => q[The number of threads for bwa sam{se,pe} command (0 for don't use -t option)],
                         );
sub _build_bwa_sam_threads {
  my($self) = @_;
  my $bwa_sam_help = slurp $self->bwa_cmd().' sampe 2>&1 |';
  if (not $bwa_sam_help =~ /-t \s+ INT \s+ number \s+ of \s+ threads/smx ) {return 0;}
  my $nthreads = $self->bwa_aln_threads;
  if ($nthreads > $self->bwa_sam_max_threads) { $nthreads = $self->bwa_sam_max_threads;}
  if($self->non_consent_split) {$nthreads /= (2|1);}# perlcritic does not like '3'
  if($nthreads ==1) {return 0;} #may as well not specify threading if num of threads = 1 ! 
  return int $nthreads;
}


=head2 bwa_sam_max_threads

The max number of threads for bwa sam{se,pe} command if bwa_sam_threads is not given
 
=cut
has 'bwa_sam_max_threads' => (isa             => 'Int',
                               is              => 'ro',
                               lazy_build      => 1,
                               documentation   => q[The maximum number of threads for bwa sam{se,pe} command if bwa_sam_max_threads is not given. default:].$DEFAULT_BWA_SAM_MAX_THREADS,
                         );
sub _build_bwa_sam_max_threads {
  return $DEFAULT_BWA_SAM_MAX_THREADS;
}


=head2 alignment_metrics_autoqc_command
 
=cut
has 'alignment_metrics_autoqc_command' => (isa             => 'Maybe[Str]',
                                           is              => 'ro',
                                           lazy_build      => 1,
                                           documentation   => 'Alignment metrics autoqc command',
                                          );
sub _build_alignment_metrics_autoqc_command {
    my $self = shift;

    my @path = splitdir $self->output_prefix;
    pop @path;
    my $in =  @path ? catdir @path : cwd;
    my $qc_out = catfile($in, q[qc]);
    if (!-d $qc_out) {
        $qc_out = $in;
    }

    my $command = sprintf '%s --check alignment_filter_metrics --id_run %i --position %i --qc_in %s --qc_out %s',
                               $QC_EXECUTABLE_NAME,
                               $self->id_run,
                               $self->position,
                               $in,
                               $qc_out;
    if (defined $self->tag_index) {
        $command .= ' --tag_index ' . $self->tag_index;
    }
    return $command;
}


=head2 bwa_sampe_options

bwa sampe command options
 
=cut
has 'bwa_sampe_options' => (isa             => 'Str',
                            is              => 'rw',
                            required        => 0,
                            documentation   => 'aligner bwa sampe command options',
                           );

=head2 java_xmx_flag

The flag specifying max memory for java, e.g. -Xmx1000m
 
=cut
has 'java_xmx_flag' => (isa             => 'Str',
                          is              => 'ro',
                          lazy_build      => 1,
                          documentation   => 'flag specifying max memory for java',
                         );
sub _build_java_xmx_flag {
    return $DEFAULT_JAVA_XMX;
}

=head2 sam_format_converter_jar

Picard  SamFormatConverter jar file

=cut

has 'sam_format_converter_jar' => (
                 is            => 'ro',
                 isa           => 'NpgCommonResolvedPathJarFile',
                 default       => $PICARD_SAM_FORMAT_CONVERTER_JAR,
                 coerce        => 1,
                 documentation => 'Picard SamFormatConverter jar file',
    );

=head2 _picard_converter_cmd

Picard  SamFormatConverter command, java with jar file

=cut

sub _picard_converter_cmd {
   my ($self) = @_;
   return $self->java_cmd . q[ ] . $self->java_xmx_flag . q[ -jar ]
        . $self->sam_format_converter_jar;
}

=head2 bam_merger_jar

=cut

has 'bam_merger_jar'     => (
           is            => 'ro',
           isa           => 'NpgCommonResolvedPathJarFile',
           coerce        => 1,
           default       => $ILLUMINA2BAM_BAM_MERGER_JAR,
           documentation => 'illumina2bam Bam merger jar file',
);

=head2 _bam_merger_cmd

illumina2bam bam merger command, java with the jar file

=cut

sub _bam_merger_cmd {
   my $self = shift;
   return $self->java_cmd . q[ ] . $self->java_xmx_flag .  q[ -jar ] . $self->bam_merger_jar;
}

=head2 alignment_filter_jar

=cut

has 'alignment_filter_jar'  => (
              is            => 'ro',
              isa           => 'NpgCommonResolvedPathJarFile',
              coerce        => 1,
              default       => $ILLUMINA2BAM_ALIGNMENT_FILTER_JAR,
              documentation => 'Illumina AlignmentFilter jar file',
);

=head2 _alignment_filter_cmd

illumina2bam AlignmentFilter cmd, java with the jar file

=cut

sub _alignment_filter_cmd {
   my $self = shift;
   return $self->java_cmd . q[ ] . $self->java_xmx_flag . q[ -jar ] . $self->alignment_filter_jar;
}

=head2 split_bam_by_chromosomes_jar

=cut

has 'split_bam_by_chromosomes_jar'     => (
              is            => 'ro',
              isa           => 'NpgCommonResolvedPathJarFile',
              coerce        => 1,
              default       => $ILLUMINA2BAM_SPLIT_BAM_BY_CHROMOSOMES_JAR,
              documentation => 'Illumina SplitBamByChromoses jar file',
);

=head2 split_by_chr_cmd

a command to execute illumina2bam SplitBamByChromosomes jar file, either with the default of X and MT, or with a specified chromosome Y

=cut

sub split_by_chr_cmd {
  my $self = shift;

  my $split_command=$self->java_cmd . q[ ] . $self->java_xmx_flag . q[ -jar ] . $self->split_bam_by_chromosomes_jar;
  if ($self->separate_y_chromosome_data()) { $split_command .= q[ S=Y V=true]; }

  return $split_command;
}


=head2 temp_dir

temp directory

=cut
has 'temp_dir'         => (isa             => 'Str',
                           is              => 'ro',
                           required        => 0,
                           default         => sub {  tempdir( CLEANUP => 1 ); },
                           documentation   => 'temp directory',
                         );

=head2 _original_pg_lines

original PG lines from input bam file

=cut

has '_original_pg_lines'   => (isa             => 'ArrayRef',
                               is              => 'rw',
                               required        => 0,
                               lazy_build      => 1,
                              );
sub _build__original_pg_lines {
   my $self = shift;

   my @pg_lines = ();

   my $samtools_view_header_cmd = $self->samtools_cmd() . q{ view -H } . $self->input();

   open my $input_bam_head_fh, q{-|}, $samtools_view_header_cmd;

   while(my $line = <$input_bam_head_fh> ){
      chomp $line;
      if( $line =~ /^[@]PG\t/mxs ){
         push @pg_lines, $line;
      }
   }
   close $input_bam_head_fh;

   return \@pg_lines;
}

=head2 _bwa_aln_commands

actual bwa aln command list

=cut
has '_bwa_aln_commands'    => (isa             => 'ArrayRef',
                               is              => 'ro',
                               required        => 0,
                               lazy_build      => 1,
                              );
sub _build__bwa_aln_commands {
   my $self = shift;

   my $command_list = [];

   if( $self->is_paired_read() ) {
      push @{$command_list}, $self->_generate_bwa_aln_command($self->reference(), $self->input(), 1);
      push @{$command_list}, $self->_generate_bwa_aln_command($self->reference(), $self->input(), 2);
   }else {
      push @{$command_list}, $self->_generate_bwa_aln_command($self->reference(), $self->input(), 0);
   }

   return $command_list;
}


=head2 _bwa_aln_pg_lines

bwa aln program pg lines for bam header

=cut
has '_bwa_aln_pg_lines'    => (isa             => 'ArrayRef',
                               is              => 'ro',
                               required        => 0,
                               lazy_build      => 1,
                              );
sub _build__bwa_aln_pg_lines {
   my $self = shift;

   my @pg_header = ();
   my @original_pg_header = @{$self->_original_pg_lines()};

   my $bwa_aln_commands = $self->_bwa_aln_commands();

   foreach my $command ( @{$bwa_aln_commands} ){
        my $pg_id    = make_unique_pg_id( [@original_pg_header, @pg_header], 'bwa_aln' ) || 'bwa_aln' ;
        my $pp_field = prev_command_field( [@original_pg_header, @pg_header] );
        my $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:bwa} . qq{$pp_field\t}
             . q{VN:} . $self->current_version($self->bwa_cmd) . qq{\t}
             . q{CL:} . $command;
        push @pg_header, $pg_line;
   }

   return \@pg_header;
}

=head2 bwa_sam_command

actual bwa sampe or samse command, with lazy build methods

The lazy built value could be rewritten using fifo if doing default and human alignments at the same time

=cut
has 'bwa_sam_command'       => (isa             => 'Str',
                                is              => 'rw',
                                required        => 0,
                                lazy_build      => 1,
                                documentation   => 'actual bwa sampe or samse command',
                             );
sub _build_bwa_sam_command {
   my $self = shift;

   my $inputs = [$self->input()];
   if( $self->is_paired_read() ){
      push @{$inputs}, $self->input();
   }

   return $self->_generate_bwa_sam_command($self->reference(), $inputs);
}


=head2 _human_bwa_aln_commands

actual bwa aln command list against human reference

=cut
has '_human_bwa_aln_commands'=> (isa             => 'ArrayRef',
                                 is              => 'ro',
                                 required        => 0,
                                 lazy_build      => 1,
                                );
sub _build__human_bwa_aln_commands {
   my $self = shift;

   my $command_list = [];
   if( $self->is_paired_read() ) {
      push @{$command_list}, $self->_generate_bwa_aln_command($self->human_reference(), $self->input(), 1, 1);
      push @{$command_list}, $self->_generate_bwa_aln_command($self->human_reference(), $self->input(), 2, 1);
   }else {
      push @{$command_list}, $self->_generate_bwa_aln_command($self->human_reference(), $self->input(), 0, 1);
   }

   return $command_list;
}

=head2 _human_bwa_aln_pg_lines

bwa aln program pg lines for bam header against human reference

=cut
has '_human_bwa_aln_pg_lines'=> (isa             => 'ArrayRef',
                                 is              => 'ro',
                                 required        => 0,
                                 lazy_build      => 1,
                                );
sub _build__human_bwa_aln_pg_lines {
   my $self = shift;

   my @pg_header = ();
   my @original_pg_header = @{$self->_original_pg_lines()};

   my $bwa_aln_commands = $self->_human_bwa_aln_commands();

   foreach my $command ( @{$bwa_aln_commands} ){
        my $pg_id    = make_unique_pg_id( [@original_pg_header, @pg_header], 'bwa_aln' ) || 'bwa_aln' ;
        my $pp_field = prev_command_field( [@original_pg_header, @pg_header] );
        my $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:bwa} . qq{$pp_field\t}
             . q{VN:} . $self->current_version($self->bwa_cmd) . qq{\t}
             . q{CL:} . $command ;
        push @pg_header, $pg_line;
   }

   return \@pg_header;
}

=head2 _human_bwa_sam_command

actual bwa sampe or samse command against human reference, with lazy build methods

The lazy built value could be rewritten using fifo if doing default and human alignments at the same time

=cut
has '_human_bwa_sam_command' => (isa             => 'Str',
                                 is              => 'rw',
                                 required        => 0,
                                 lazy_build      => 1,
                             );
sub _build__human_bwa_sam_command {
   my $self = shift;

   my $inputs = [$self->input()];
   if( $self->is_paired_read() ){
      push @{$inputs}, $self->input();
   }
   return $self->_generate_bwa_sam_command($self->human_reference(), $inputs, 1);
}

=head2 BUILD

Runs just before returning the instance reference by ->new().
Checks the consistency of arguments and re-sets them if necessary.

=cut

sub BUILD {
  my $self = shift;

  $self->log('The input '. ($self->no_alignment ? 'does not need' : 'needs') .' alignment');
  $self->log('phix '. ($self->spiked_phix_split ? 'will' : 'will not') .' be split out');
  $self->log('Y chromosome '. ($self->separate_y_chromosome_data ? 'will' : 'will not') .' be split out');
  $self->log('The input '. ($self->non_consent_split ? 'contains':'does not contain') .' nonconsented human data');
  $self->log('The input '. ($self->is_paired_read ? 'is' : 'is not') .' paired read');

  if ($self->non_consent_split && $self->separate_y_chromosome_data() ) {
      croak('separate_y_chromosome_data and non_consent_split cannot both be true.');
  }
  if ($self->non_consent_split && $self->contains_nonconsented_xahuman) {
    $self->log('Both non_consent_split and contains_nonconsented_xahuman flags are true.
                Setting contains_nonconsented_xahuman to false');
    $self->_set_contains_nonconsented_xahuman(0);

  }

  my $reference = $self->reference;
  $self->log('Reference' . ($reference ? ": $reference" : ' not found'));

  if ($self->contains_nonconsented_xahuman) {
    if (!$reference || $reference !~ /Homo_sapiens/smx) {
      croak 'contains_nonconsented_xahuman is true, human reference should be used';
    }
    if ($self->separate_y_chromosome_data() ) {
      croak('separate_y_chromosome_data and contains_nonconsented_xahuman flags cannot both be true.');
    }
    if ($self->no_alignment) {
      croak('no_alignment and contains_nonconsented_xahuman flags cannot both be true.');
    }
    $self->log('The input contains nonconsented X and autosomal human,
                setting nonconsented_file_name_suffix to "xahuman"');
    $self->nonconsented_file_name_suffix('xahuman');
  }

  if ($self->separate_y_chromosome_data() ) {
    if (!$reference || $reference !~ /Homo_sapiens/smx) {
      croak 'separate_y_chromosome_data is true, human reference should be used';
    }
    if ($self->no_alignment) {
      croak('no_alignment and separate_y_chromosome_data flags cannot both be true.');
    }
    $self->log('Setting suffix to yhuman');
    $self->nonconsented_file_name_suffix('yhuman');
  }
  $self->log("************************************************************\n");

  return;
}


sub _generate_bwa_aln_command {
   my ($self, $ref, $input, $read, $non_consented_split_alignment) = @_;

   my $command = $self->bwa_cmd()
                . q{ aln };

   if($self->bwa_aln_options()){
     $command  .=  $self->bwa_aln_options();
   }

   $command    .= qq{ $ref}
                . qq{ -b$read $input}
                . q{ > } . $self->temp_dir() .q{/};
   if($non_consented_split_alignment){
     $command .= q{human};
   }
   $command .= $read . q{.sai};

   return $command;
}

sub _generate_bwa_sam_command {
   my ($self, $ref, $inputs, $non_consented_split_alignment) = @_;

   my $command = $self->bwa_cmd();

   if( scalar @{$inputs} == 2 ){

      $command .= q{ sampe};
      if( my $nthread = $self->bwa_sam_threads ){
         $command .= qq{ -t $nthread};
      }
      if( $self->bwa_sampe_options() ){
         $command .= q{ } . $self->bwa_sampe_options();
      }
      if( $non_consented_split_alignment ){
         $command .= q{ } . $ref;
         $command .= q{ } . $self->temp_dir() . q{/human1.sai};
         $command .= q{ } . $self->temp_dir() . q{/human2.sai};
      }else{
         $command .= q{ } . $ref;
         $command .= q{ } . $self->temp_dir() . q{/1.sai};
         $command .= q{ } . $self->temp_dir() . q{/2.sai};
      }
   }elsif( scalar @{$inputs} == 1 ){

      $command .= q{ samse};
      if( my $nthread = $self->bwa_sam_threads ){
         $command .= qq{ -t $nthread};
      }
      if( $non_consented_split_alignment ){
         $command .= q{ } . $ref;
         $command .= q{ } . $self->temp_dir() . q{/human0.sai};
      }else{
         $command .= q{ } . $ref;
         $command .= q{ } . $self->temp_dir() . q{/0.sai};
      }
   }else {
      croak 'The number of input bam is not 1 or 2' . (join qq{\n}, @{$inputs} );
   }

   $command .= q{ } . join q{ }, @{$inputs};

   return $command;
}

=head2 generate

MAIN METHOD to CALL
 
=cut

sub generate {
  my $self = shift;

  if(! -e $self->input() ){
    croak 'The input bam file does not exist: ' . $self->input();
  }

  $self->log( q{Starting doing alignments} );
  my $number_aln = $self->_run_bwa_alns();
  $self->log( "$number_aln  bwa aln processes finished" );

  $self->log( q{Starting doing bwa sampe or samse, and filter output by alignments} );
  $self->_run_bwa_sam();

  if($self->no_alignment() &&  !$self->non_consent_split() && ! $self->spiked_phix_split() ) {
    $self->_soft_link_output();
  }

  my $qc_command = $self->alignment_metrics_autoqc_command;
  if ( ($self->non_consent_split() || $self->spiked_phix_split()) && $qc_command ) {
    if (which($QC_EXECUTABLE_NAME)) {
      $self->log("Executing $qc_command");
      system $qc_command;
      if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
        croak "Alignment filter autoqc check failed: $ERRNO\n$qc_command";
      }
    } else {
      carp 'autoqc modules not available, not parsing the output of the alignment filter'
    }
  }

  remove_tree( $self->temp_dir(), {keep_root => 1, verbose => 1} );

  $self->log( q{Sorting output and mark duplicates now} );

  my @markduplicates_commands = ();
  push @markduplicates_commands, $self->_generate_markduplicates_cmd( $self->output(), undef, $self->reference() );

  if( $self->non_consent_split() ) {
    # the sample reference is not human, non consented output aligned to the default human reference
    push @markduplicates_commands, $self->_generate_markduplicates_cmd($self->non_consented_output(), $self->nonconsented_file_name_suffix, $self->human_reference());
  }

  if( $self->contains_nonconsented_xahuman() || $self->separate_y_chromosome_data() ) {
    # the sample reference is human, non consented output aligned to the sample reference
    push @markduplicates_commands, $self->_generate_markduplicates_cmd($self->non_consented_output(), $self->nonconsented_file_name_suffix, $self->reference());
  }

  if( $self->spiked_phix_split() ) {
    push @markduplicates_commands, $self->_generate_markduplicates_cmd($self->spiked_phix_output(), q{phix}, $self->phix_reference());
  }

  foreach my $mk_cmd ( @markduplicates_commands ){
    $self->log( $mk_cmd );
    system $mk_cmd;
    if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
      croak "Mark duplicates command failed: $ERRNO";
    }
  }
  $self->log('Finished in BAM_Alignment!');

  return;
}

sub _soft_link_output{
    my $self = shift;

    $self->log('The input file does not need alignment and filtering, create a soft link now');

    my $cwd = cwd();

    my ($output_vol, $out_dir, $output_file) = File::Spec->splitpath($self->output());
    chdir $out_dir;

    my $input_rel_path = File::Spec->abs2rel($self->input());

    if( -e $output_file ){
       unlink $output_file;
    }
    symlink $input_rel_path, $output_file;

    my $output_md5 = $output_file . q{.md5};
    my $input_md5  = $input_rel_path . q{.md5};

    if( -e $output_md5 ) {
       unlink $output_md5;
    }
    symlink $input_md5, $output_md5;

    chdir $cwd;

    return;
}

sub _run_bwa_alns {
  my $self = shift;

  my @aln_commands = ();

  if( !$self->no_alignment() ){
    push @aln_commands, @{$self->_bwa_aln_commands()};
  }

  if( $self->non_consent_split() ){
    push @aln_commands, @{$self->_human_bwa_aln_commands()};
  }

  if( scalar @aln_commands == 0 ){
    $self->log('No alignment needed');
    return 0;
  }

  foreach my $command ( @aln_commands ) {

    $self->log($command);
    system $command;
    if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
      croak "bwa align failed: $ERRNO\n$command";
    }
  }

  return scalar @aln_commands;
}

sub _run_bwa_sam { ##no critic (Subroutines/ProhibitExcessComplexity)
  my $self = shift;

  my $sam_commands = $self->_bwa_sam_commands();

  my $fifo2chr_split = ( $self->contains_nonconsented_xahuman() || $self->separate_y_chromosome_data() )? catfile($self->temp_dir, 'fifo2split_by_chr.bam') : undef;

  my $alignment_filter_cmd = $self->_bam_alignment_filter_command($fifo2chr_split);
  if( $alignment_filter_cmd ) {
    $sam_commands->{$alignment_filter_cmd} = q{alignment_filter};
  } elsif ( $fifo2chr_split ) {
    my $mkfifo_cmd = join q[ ], 'mkfifo', $self->_bam_output_fifo();
    $self->log($mkfifo_cmd);
    system $mkfifo_cmd;
    if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
      croak "$mkfifo_cmd failed: $ERRNO";
    }
    $fifo2chr_split = $self->_bam_output_fifo();
  }

  if ($fifo2chr_split) {
    my $call = q();
    $call = sprintf  '%s VALIDATION_STRINGENCY=SILENT I=%s TARGET_PATH=%s EXCLUDED_PATH=%s',
                         $self->split_by_chr_cmd,
                         $fifo2chr_split,
                         $self->output,
                         $self->non_consented_output;
    $sam_commands->{$call} = 'chr_split';
  }

  my $pm = Parallel::ForkManager->new( scalar keys %{$sam_commands} );

  $pm->run_on_finish(
    sub { my ($pid, $exit_code, $ident) = @_;
       #exit code shift 8?
       if($exit_code){
          croak "PID $pid and exit code: $exit_code. Fail: $ident";
       }else{
          $self->log( "PID $pid and exit code: $exit_code. Success: $ident") ;
       }
    }
  );

  $pm->run_on_start(
    sub { my ($pid,$ident)=@_;
      $self->log( "Job $pid started: $ident" );
    }
  );

  my $count = 0;
  foreach my $command (sort keys %{$sam_commands} ){

      $count++;
      $self->log("Fork $count: $command");
      $pm->start("Fork $count $command") and next;

      my $program =  $sam_commands->{$command};
      ##no critic (ProhibitCascadingIfElse)
      if( $program eq q{cat} ){ #cat | tee
# YET TO BE TESTED, was: my $sam_cmd_rs = system qq[set -o pipefail; $command];
          my $sam_cmd_rs = system qq[/bin/bash -c "set -o pipefail; $command"];
          if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
             croak "bwa samse or sampe or creating input pipe failed: $ERRNO\n$command";
          }
      } elsif ( $program eq q{alignment} ){
          my $output_bam = ($alignment_filter_cmd || $fifo2chr_split) ? $self->_bam_output_fifo() : $self->output();
          my $bwa_aln_pg_lines = $self->_bwa_aln_pg_lines();
          my $ref_dict = $self->ref_dict();
          my $ref_dict_contents;
          if($ref_dict && -e $ref_dict){
            $ref_dict_contents = slurp $ref_dict;
          }
          $self->_bwa_sam_out_to_bam($command, $bwa_aln_pg_lines, $output_bam,  $ref_dict_contents);
      } elsif ( $program eq q{human_alignment} ){
          my $output_bam = $alignment_filter_cmd ? $self->_human_bam_output_fifo() : $self->non_consented_output();
          my $bwa_aln_pg_lines = $self->_human_bwa_aln_pg_lines();
          my $human_ref_dict = $self->human_ref_dict();
          my $human_ref_dict_contents;
          if($human_ref_dict && -e $human_ref_dict){
            $human_ref_dict_contents = slurp $human_ref_dict;
          }
          $self->_bwa_sam_out_to_bam($command, $bwa_aln_pg_lines, $output_bam, $human_ref_dict_contents);
      } elsif ($program eq q{alignment_filter}){
          system $command;
          if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
             croak "Alignment filter command failed: $ERRNO\n$command";
          }
      } elsif ($program eq 'chr_split') {
	  system $command;
          if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
             croak "Split chr command failed: $ERRNO\n$command";
          }
      }

      $pm->finish;
  }
  $pm->wait_all_children;

  return 1;
}

sub _bwa_sam_out_to_bam {
    my ($self, $command, $bwa_aln_pg_lines, $output_bam, $ref_dict) = @_;

    #make new bwa_sam PG ID based on all existed PG
    my @old_pg_lines = @{$self->_original_pg_lines()};
    push @old_pg_lines, @{$bwa_aln_pg_lines};
    my $pg_id    = make_unique_pg_id(\@old_pg_lines, 'bwa_sam' ) || 'bwa_sam' ;

    my $pp_field = prev_command_field( $bwa_aln_pg_lines);

    my $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:bwa} . qq{$pp_field\t}
             . q{VN:} . $self->current_version($self->bwa_cmd) . qq{\t}
             . q{CL:} . $command;

    #only keep bwa aln and bwa sam PG together
    my @pg_lines = @{$bwa_aln_pg_lines};
    push @pg_lines, $pg_line;
    push @old_pg_lines, $pg_line;

    open my $bwa_sam_output_fh, q{-|}, $command; ## no critic (InputOutput::RequireBriefOpen);

    #pass all PG lines, and the list of new PG lines to get more new PG lines
    my ($output_bam_fh, $new_pg_lines) = $self->generate_output_bam_fh( $output_bam, \@pg_lines, \@old_pg_lines );

    if( $ref_dict ){
       print {$output_bam_fh} $ref_dict or croak $OS_ERROR;
    }

    #write out all bwa sam output header except PG lines
    my $line = <$bwa_sam_output_fh>;
    while ( $line && $line =~ /^[@]/mxs ) {
        if( ( $line !~ /^[@]PG/mxs ) && !$ref_dict ){
           print {$output_bam_fh} $line or croak $OS_ERROR;
        }
        $line = <$bwa_sam_output_fh>;
    }


    #write out all new pg lines, the old pg lines will be merged in from the input bam file later
    foreach my $pg_line ( @{$new_pg_lines} ){

            print {$output_bam_fh} $pg_line, "\n" or croak $OS_ERROR;;
    }

    #write out the alignment records
    while( $line ){
            print {$output_bam_fh} $line or croak $OS_ERROR;
            $line = <$bwa_sam_output_fh>;
    }

    close $bwa_sam_output_fh;# or croak "Couldn't close output sam: $OS_ERROR";
    close $output_bam_fh;    # or croak "Couldn't close output bam: $OS_ERROR";
    return;
}

sub _bwa_sam_commands{
  my $self = shift;
     my %sam_commands = ();
     if( $self->non_consent_split() && !$self->no_alignment() ) {
        my  @input_pipes = ();

        if( $self->is_paired_read() ){

           my $fifo1 = $self->temp_dir().q{/fifo1.bam};
           push @input_pipes, $fifo1;
           my $fifo2 = $self->temp_dir().q{/fifo2.bam};
           push @input_pipes, $fifo2;

           my $human_fifo1 = $self->temp_dir().q{/human_fifo1.bam};
           push @input_pipes, $human_fifo1;
           my $human_fifo2 = $self->temp_dir().q{/human_fifo2.bam};
           push @input_pipes, $human_fifo2;

           my $cat_command1 = q{cat } . $self->input() . qq{ | tee $fifo1 $human_fifo1 > /dev/null};
           $sam_commands{$cat_command1} = q{cat};

           my $cat_command2 = q{cat } . $self->input() . qq{ | tee $fifo2 $human_fifo2 > /dev/null};
           $sam_commands{$cat_command2} = q{cat};

           my $sam_command = $self->_generate_bwa_sam_command($self->reference(),
                                                              [$fifo1, $fifo2]
                                                             );

           $self->bwa_sam_command($sam_command);

           $sam_commands{$sam_command} = q{alignment};

           my $sam_command_human = $self->_generate_bwa_sam_command($self->human_reference(),
                                                              [$human_fifo1, $human_fifo2],
                                                              1
                                                            );
           $self->_human_bwa_sam_command( $sam_command_human );
           $sam_commands{$sam_command_human} = q{human_alignment};

        }else{

           my $fifo0 = $self->temp_dir().q{/fifo0.bam};
           push @input_pipes, $fifo0;

           my $human_fifo0 = $self->temp_dir().q{/human_fifo0.bam};
           push @input_pipes, $human_fifo0;

           my $cat_command = q{cat } . $self->input() . qq{ | tee $fifo0 $human_fifo0 > /dev/null};
           $sam_commands{$cat_command} = q{cat};

           my $sam_command = $self->_generate_bwa_sam_command($self->reference(),
                                                              [$fifo0]
                                                             );
           $self->bwa_sam_command( $sam_command );
           $sam_commands{$sam_command} = q{alignment};

           my $sam_command_human = $self->_generate_bwa_sam_command($self->human_reference(),
                                                              [$human_fifo0],
                                                              1
                                                             );
           $self->_human_bwa_sam_command($sam_command_human);
           $sam_commands{$sam_command_human} = q{human_alignment};
        }

        if( scalar  @input_pipes ){
             my $mkfifo_cmd = "mkfifo @input_pipes";
             $self->log($mkfifo_cmd);
             system $mkfifo_cmd;

             if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
                 croak "mkfifo failed: $ERRNO";
             }
        }

   }elsif( $self->non_consent_split() ){
      $sam_commands{$self->_human_bwa_sam_command()} = q{human_alignment};
   }elsif( !$self->no_alignment() ){
      $sam_commands{$self->bwa_sam_command()} = q{alignment};
   }

   return \%sam_commands;
}

sub _bam_alignment_filter_command {
    my ($self, $fifo_output4target) = @_;

    my $metrics = $self->output() . '_alignment_filter_metrics.json';
    my $filter_cmd = $self->_alignment_filter_cmd() .
                            q{ VALIDATION_STRINGENCY=SILENT CREATE_MD5_FILE=true} .
                           qq{ METRICS_FILE=$metrics};
    my @inputs = ();
    my @outputs = ();
    my @fifos = ();

    if( $self->spiked_phix_split() ){
       push @inputs, $self->input();
       push @outputs, $self->spiked_phix_output();
    }

    if( $self->non_consent_split() ){
       push @inputs, $self->_human_bam_output_fifo();
       push @outputs, $self->non_consented_output();
       push @fifos, $self->_human_bam_output_fifo();
    }

    if ( scalar @inputs == 0 ){
       return;
    }

    if( $self->no_alignment() ){
       $filter_cmd .= q{ UNALIGNED=} . $self->output();
    }else{
       push @inputs, $self->_bam_output_fifo();
       push @fifos, $self->_bam_output_fifo();
       if ($fifo_output4target) {
          push @outputs, $fifo_output4target;
	  push @fifos, $fifo_output4target;
       } else {
          push @outputs, $self->output();
       }
    }

    foreach my $input (@inputs){
       $filter_cmd .= qq{ IN=$input};
    }
    foreach my $output (@outputs){
       $filter_cmd .= qq{ OUT=$output};
    }

    if( scalar @fifos == 0 ){
       return $filter_cmd;
    }

    eval{
       unlink @fifos;
       1;
    } or do {
     carp $EVAL_ERROR;
    };

    my $mkfifo_cmd = "mkfifo @fifos";
    $self->log($mkfifo_cmd);
    system $mkfifo_cmd;

    if( $CHILD_ERROR >> $EXIT_CODE_SHIFT ){
        croak "mkfifo failed: $ERRNO";
    }

    return $filter_cmd;
}

=head2 generate_output_bam_fh

=cut

sub generate_output_bam_fh {
    my ($self, $output_bam, $pg_records, $old_pg_records) = @_;

    my @pg_header = @{$pg_records};
    my @old_pg_header = @{$old_pg_records};

    my $sam2bam_command = $self->_picard_converter_cmd()
                       . q{ VALIDATION_STRINGENCY=SILENT}
                       . q{ INPUT=/dev/stdin OUTPUT=/dev/stdout COMPRESSION_LEVEL=0};

    my $pg_id    = make_unique_pg_id( \@old_pg_header, 'Picard_SamFormatConverter' ) || 'Picard_SamFormatConverter' ;
    my $pp_field = prev_command_field( \@pg_header );
    my $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:SamFormatConverter} . qq{$pp_field\t}
             . q{VN:} . $self->current_version($self->sam_format_converter_jar) . qq{\t}
             . q{CL:} . $sam2bam_command;
    push @pg_header, $pg_line;
    push @old_pg_header, $pg_line;

    my $bam_fixmate_command = $self->samtools_cmd() . q{ fixmate - -};

    $pg_id     = make_unique_pg_id( \@old_pg_header, 'samtools_fixmate' ) || 'samtools_fixmate';
    $pp_field  = prev_command_field( \@pg_header );
    $pg_line = q{@} . qq{PG\tID:$pg_id\tPN:samtools} . qq{$pp_field\t}
                . q{VN:} . $self->current_version($self->samtools_cmd) . qq{\t}
                . q{CL:} . $bam_fixmate_command ;
    push @pg_header, $pg_line;
    push @old_pg_header, $pg_line;

    my $bam_merger_cmd = $self->_bam_merger_cmd()
                      . q{ ALIGNED=/dev/stdin I=} . $self->input()
                      . qq{ O=$output_bam}
                      . q{ PG=NULL KEEP=true KEEP_PG=true CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT};

    my $command = $sam2bam_command . q{ | } . $bam_fixmate_command . q{ | } . $bam_merger_cmd;
    $self->log($command);

    open my $bam_out_fh, q[|-], qq[/bin/bash -c "set -o pipefail; $command"];

    return ($bam_out_fh, \@pg_header);
}

=head2 prev_command_field

find the last PG ID in the list

=cut
sub prev_command_field {
    my ($pg_lines) = @_;

    my @lines = validate_pg_arg($pg_lines);
    return q{} if scalar @lines == 0;

    my $last_line = $lines[-1];

    my $pp_name;
    if ( $last_line =~ m/ ID: ( [^\t\n]+ ) /msx ) {
        $pp_name = $1;
    }
    else {
        croak "Error inferring PP field from $last_line" if !$pp_name;
    }

    return qq{\tPP:$pp_name};
}

=head2 make_unique_pg_id

=cut

sub make_unique_pg_id {
    my ( $pg_lines, $base_id ) = @_;

    croak 'Base id argument required' if !$base_id;
    my @lines = validate_pg_arg($pg_lines);
    return $base_id if !scalar @lines;

    my %id_list = ();
    foreach my $pg_line (@lines) {

       my ($id) = $pg_line =~ /ID\:($base_id?\w+)?[\t\n]/mxs;
       if($id){
          $id_list{$id} = 1;
       }
    }

    my $new_id = $base_id;
    my $count = 0;
    while($id_list{$new_id}){

      $count++;
      $new_id = $base_id . q{_} . $count
    }

    return $new_id;
}

=head2 validate_pg_arg

=cut
sub validate_pg_arg {
    my ($pg_lines) = @_;

    croak 'Arrayref arguments only' if ref $pg_lines ne 'ARRAY';

    return @{$pg_lines};
}

sub _generate_markduplicates_cmd {
   my ($self, $bam, $alignment_filter, $reference) = @_;

   FindBin::again();

   my $cmd = catfile($Bin, $BAM_MARK_DUPLICATES_CMD) . qq{ --input_bam $bam};
   if($self->bamcheck_flags) {
      $cmd .= q{ --bamcheck_flags "} . $self->bamcheck_flags . q{"}; # note: raw string, no validation. Should be used carefully
   }

   if ($reference) {
      $cmd .= q{ --reference } . $reference;
   }
   $cmd .= q{ --id_run } . $self->id_run();
   $cmd .= q{ --position } . $self->position();
   if(defined $self->tag_index() ){
      $cmd .= q{ --tag_index } . $self->tag_index();
   }

   my($filename, $path, $suffix) = fileparse($bam, q{.bam});

   my $output_bam  = catfile( $path, $filename.q{_mk.bam} );
   my $qc_dir      = catfile( $path, q{qc} );
   if(!-d $qc_dir){
      mkdir $qc_dir;
   }

   if( $alignment_filter ) {
      $cmd .= qq{ --subset $alignment_filter};
   }

   $cmd .= qq{ --output_bam $output_bam};
   $cmd .= qq{ --metrics_json_dir $qc_dir};

   return $cmd;
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

=item MooseX::StrictConstructor

=item File::Spec::Functions

=item File::Spec

=item File::Temp

=item File::Path

=item File::Which

=item Cwd

=item File::Basename

=item Carp

=item English -no_match_vars

=item autodie

=item Perl6::Slurp

=item use FindBin qw($Bin)

=item Readonly

=item Parallel::ForkManager

=item npg_common::roles::log

=item npg_common::roles::software_location

=item npg_tracking::glossary::run

=item npg_tracking::glossary::lane

=item npg_tracking::util::abs_path

=item st::api::lims

=item npg_tracking::data::reference

=item npg_tracking::data::reference::info

=item npg_tracking::illumina::run::folder

item npg_tracking::illumina::run::long_info

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
