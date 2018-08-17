#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib ( -d "$Bin/../lib/perl5" ? "$Bin/../lib/perl5" : "$Bin/../lib" );
use Carp;
use English qw(-no_match_vars);
use File::chdir;
use File::Path qw(make_path remove_tree);
use File::Basename;
use File::Copy;
use IO::File;
use Getopt::Long;
use List::Util qw(sum none);
use autodie qw(:all);
use Readonly;
use Pod::Usage qw(pod2usage);
use feature qw(say);
use POSIX qw(strftime);
use npg_tracking::util::abs_path qw(abs_path);

our $VERSION = '0';

Readonly::Scalar my $DIR_PERM => oct 755;
Readonly::Scalar my $SKIP_STRING => "\n*** Skipping to next aligner ***\n\n";
Readonly::Scalar my $LATEST_SALMON => "salmon0_10";
Readonly::Scalar my $DEFAULT_ANNOT_DIR => q[gtf];
Readonly::Scalar my $ONE_HUNDRED => 100;
Readonly::Scalar my $MIN_RATE => 0.5;
Readonly::Scalar my $MAX_LINES => 10;
Readonly::Scalar my $FASTA_FILE_REGEX => qr{[.]f(?:ast|n)?a\Z}imsx;
Readonly::Scalar my $GTF_FILE_REGEX => qr{[.]gtf\Z}imsx;
Readonly::Scalar my $GTF_GFF_FILE_REGEX => qr{[.]g(?:t|f)f3?\Z}imsx;
Readonly::Scalar my $BT2_FILE_REGEX => qr{(?<![.]rev)[.][1-9]+[.]bt2\Z}imsx;
Readonly::Scalar my $DICT_FILE_REGEX => qr{[.]dict\Z}imsx;
Readonly::Scalar my $GENCODE_ATTR_REGEX => qr{\A(?:\s*(?:.+)\s\"(?:.+)\";)+\Z}ismx;

# All the tools that can be passed as options must be listed here.
# Explicit order is necessary for some tools (e.g. Salmon: oldest -> newest)
Readonly::Array  my @TOOLS => qw{
                        cellranger
                        fasta
                        gtf
                        rna_seqc
                        salmon0_8
                        salmon0_10
                        salmon
                        tophat2
                    };


my ($help, $verbose);
my ($working_dir, $ref_genome_dir, $annotation_dir, $bt2_dir, $dict_dir);
my $save_discarded = 0;

my %build = map { $_ => 0 } @TOOLS;
pod2usage(2) if (! @ARGV);
GetOptions(
    \%build,
    keys %build,
    'annotation=s' => \$annotation_dir,
    'bowtie2=s'=> \$bt2_dir,
    'dictionary=s' => \$dict_dir,
    'genome=s' => \$ref_genome_dir,
    'help',
) || pod2usage(2);

pod2usage(0) if $help;
pod2usage(-message => q[Mandatory argument `-genome' is missing], -verbose => 99,
          -sections => [qw{OPTIONS}], -exitval => 1) unless (defined $ref_genome_dir);
pod2usage(-message => q[Mandatory argument `-annotation' is missing], -verbose => 99,
          -sections => [qw{OPTIONS}], -exitval => 1) unless (defined $annotation_dir);

# if no tools were specified do all
if (sum( values %build ) == 0) {
    foreach my $k (keys %build) { $build{$k} = 1; }
}

my $status = main();

exit $status;


sub main {
    my @genome = inspect_dir($ref_genome_dir, $FASTA_FILE_REGEX);
    if (! @genome) {
        croak qq[Cannot access $ref_genome_dir: no such directory];
    } elsif (scalar @genome != 1){
        croak q[One and only one fasta file expected.];
    }
    my $ref_genome_file_abs_path = abs_path(qq[$ref_genome_dir/$genome[0]]);

    my @annotation = inspect_dir($annotation_dir, $GTF_GFF_FILE_REGEX);
    if (! @annotation) { 
        croak qq[Cannot access $annotation_dir: no such directory];
    } elsif (scalar @annotation != 1){
        croak q[One and only one annotation file expected.];
    }
    my $annot_file_abs_path = abs_path(qq[$annotation_dir/$annotation[0]]);

    $working_dir = $CWD;

    my (%arguments, %index_command, @failed, @passed, %subs);

    # Build the indexing arguments
    $subs{'cellranger'} = \&cellranger;
    $subs{'fasta'}      = \&fasta;
    $subs{'rna_seqc'}   = \&rna_seqc;
    $subs{'salmon0_8'}  = \&salmon;
    $subs{'salmon0_10'} = \&salmon;
    $subs{'salmon'}     = \&salmon_latest;
    $subs{'tophat2'}    = \&tophat2;

    # Create common annotation in local gtf directory first then
    # remove it from the build list in case it's there via CLI
    if (! gtf('gtf', $ref_genome_file_abs_path, $annot_file_abs_path)) {
        croak q[Failed to create custom annotation gtf.];
    } else {
        push @passed, 'gtf';
    }
    delete $build{'gtf'};
    $CWD = $working_dir;

    # Iterate over subs for aligners/tools
    foreach my $tool (@TOOLS){
        if ($build{$tool}) {
            if (!$subs{$tool}->($tool, $ref_genome_file_abs_path)) {
                push @failed, $tool;
            } else {
                push @passed, $tool;
            }
            $CWD = $working_dir;
        }
    }

    return 1;
}


sub clean_slate {
    my ($dirname) = @_;
    ( -e $dirname ) && remove_tree($dirname);
    make_path( $dirname, { mode => $DIR_PERM } );
    return;
}


sub create_symlink{
    my($source, $target, $tool) = @_;
    if (!-e $source || (lstat $source && ! stat $source)) {
        carp qq[**** symlinking $source failed: No such file or directory or broken link];
        return 0;
    }
    eval {
        (-e $target) && (unlink $target);
        symlink $source, $target;
        1;
    } or do {
        carp qq[**** symlinking $target file failed: $EVAL_ERROR $SKIP_STRING];
        return 0;
    };
    return 1;
}


sub execute_command {
    my ($tool, $tool_command) = @_;
    eval { 
        system $tool_command;
        1;
    } or do {
        carp qq[*** $tool execution failed: $EVAL_ERROR $SKIP_STRING];
        return 0;
    };
    return 1;
}


sub inspect_dir {
    my ($dir_name, $file_pattern) = @_;
    my $dh;
    if (! -d $dir_name){
        return;
    }
    opendir $dh, $dir_name;
    my @files = grep {m/$file_pattern\z/imsx} readdir $dh;
    closedir $dh;
    return @files;
}


sub log_timestamp {
    my ($state, $job) = @_;
    my $message = sprintf "%s: %s %% %s\n",
                  uc $state, $job, strftime("%Y-%m-%d %H:%M:%S\n", localtime);
    print {*STDOUT} $message or carp $OS_ERROR;
    return;
}

sub get_common_gtf {
    # get annotation from local gtf directory
    $CWD = $working_dir;
    my @gtf = inspect_dir(q[gtf], $GTF_FILE_REGEX);
    if (! @gtf) {
        croak qq[Cannot access common annotation file in $CWD];
    }
    my $gtf_abs_path = abs_path(qq[gtf/$gtf[0]]);
    return $gtf_abs_path;
}
####################
# Tool methods
####################


sub cellranger {
    my ($tool, $genome_abs_path) = @_;

#    $CWD = $working_dir;

    log_timestamp('start', $tool);

    my $annotation_abs_path = get_common_gtf;
    my $gtf_name = fileparse($annotation_abs_path, $GTF_FILE_REGEX);

    ( -e '10X' ) && remove_tree('10X');

    say q[*** Cellranger index creation: creating tmp filtered gtf];
    my $cellranger_command = qq[cellranger mkgtf $annotation_abs_path ].
        qq[tmp.$gtf_name-filtered.gtf ].
        qq[--attribute=gene_biotype:protein_coding];

    if(execute_command('cellranger - mkgtf', $cellranger_command)){
        say q[*** Cellranger index creation: creating ref index];
        $cellranger_command = qq[cellranger mkref ].
            qq[--genome=10X ].
            qq[--fasta=$genome_abs_path ].
            qq[--genes=tmp.$gtf_name-filtered.gtf ].
            qq[--nthreads=16 --memgb=32]; #TODO: this will be documented in the POD but
                                          #      there should be a better way than hard-coded
        if(execute_command('cellranger - mkgtf', $cellranger_command)){
            say q[*** Cellranger index creation: clean up];
            make_path(q[10X/logs], { mode => $DIR_PERM });
            move(q[Log.out], q[10X/logs/cellranger_mkref_Log.out]);
            move(q[cellranger_mkref.log], q[10X/logs/]);
            move(q[cellranger_mkgtf.log], q[10X/logs/]);
            unlink qq[tmp.$gtf_name-filtered.gtf];
            chmod $DIR_PERM, q[10X];
        } else {
            log_timestamp('abort', 'cellranger - mkgtf');
            return 0;
        }
    } else {
        log_timestamp('abort', 'cellranger - mkref');
        return 0;
    };

    log_timestamp('stop', $tool);
    return 1;
}


sub fasta {
    my ($tool, $genome_abs_path) = @_;

#    $CWD = $working_dir;
    log_timestamp('start', $tool);
    clean_slate('fasta');

    my $annotation_abs_path = get_common_gtf;
    my $genome_name = fileparse($genome_abs_path, $FASTA_FILE_REGEX);

    my $fasta_command = qq[gffread -w fasta/$genome_name.transcripts.fa ].
                        qq[-g $genome_abs_path $annotation_abs_path];

    if(! execute_command('fasta', $fasta_command)){
        log_timestamp('abort', $tool);
        return 0;
    }
    
    log_timestamp('stop', $tool);
    return 1;
}


sub gtf {
    my ($tool, $genome_abs_path, $annotation_abs_path) = @_;

    log_timestamp('start', $tool);

    my $annotation_name = fileparse($annotation_abs_path, $GTF_GFF_FILE_REGEX);
    my $putative_gtf_abs_path = abs_path(qq[gtf/$annotation_name.gtf]) // q[];
    if($annotation_abs_path eq $putative_gtf_abs_path) {
        carp qq[*** Ref annotation gtf creation: ].
             qq[Source and target annotation files are the same: nothing else to do here];
        log_timestamp('stop', $tool);
        return 1;
    }

    clean_slate('gtf');
    $CWD = q[./gtf];

    # verify whether annotation file
    # is in GTF format or not
    my $gtf_format = 0;
    my $num_lines = 0;
    my $fh_probe_annot = IO::File->new($annotation_abs_path, 'r');
    while (my $line = $fh_probe_annot->getline) {
        next if ($line =~ /\A\#/smx);
        my @columns = split /\t/smx, $line;
        if (scalar @columns >= 9){
            if ($columns[8] =~ $GENCODE_ATTR_REGEX){
                $gtf_format += 1;
            }
        }
        $num_lines += 1;
        last if ($num_lines >= $MAX_LINES || ($gtf_format >= $MAX_LINES));
    }
    $fh_probe_annot->close();

    # if the annotation file is in GTF format,
    # copy it and move on ...
    if ($gtf_format >= $MAX_LINES){
        carp qq[*** Ref annotation gtf creation: ].
             qq[Source is in gtf format: copying file];
        copy($annotation_abs_path, q[.]);
        log_timestamp('stop', $tool);
        return 1;
    }

    # ... otherwise, convert it to GTF
    my $gtf_command = qq[gffread $annotation_abs_path ].
                      qq[-T -o $annotation_name.gtf];
    if(! execute_command('gtf', $gtf_command)){
        log_timestamp('abort', $tool);
        return 0;
    }

    log_timestamp('stop', $tool);
    return 1;
}


sub rna_seqc {
    my ($tool, $genome_abs_path) = @_;

    # validate annotation file is RNA-SeQC-compliant:
    #- ensembl gtf format assumed (gene type in column 2)
    #- all entries must have transcript id
    #- all entries must be on contig in reference sequence dict.
    #- files must be called *.gtf
    log_timestamp('start', $tool);

    my $annotation_abs_path = get_common_gtf;
    my $gtf_name = fileparse($annotation_abs_path, $GTF_FILE_REGEX);

    # a reference genome dictionary is required
    # with matching contig names. load it onto a hash
    my $dict_file_abs_path;
    if (defined $dict_dir) {
        my @piccard_dict = inspect_dir($dict_dir, $DICT_FILE_REGEX);
        if (! @piccard_dict) {
            carp q[*** RNA-SeQC annotation creation: ].
                 qq[Cannot access $dict_dir: no such directory];
            log_timestamp('abort', $tool);
            return 0;
        } elsif (scalar @piccard_dict != 1){
            carp q[*** RNA-SeQC annotation creation: ].
                 q[One and only one Piccard dictionary file expected.];
            log_timestamp('abort', $tool);
            return 0;
        }
        $dict_file_abs_path = abs_path(qq[$dict_dir/$piccard_dict[0]]);
    } else {
        carp q[*** RNA-SeQC annotation creation: ].
             q[Path to picard dictionary file required ];
        log_timestamp('abort', $tool);
        return 0;
    }
    my %ref_dict = ();
    my $fh_dict = IO::File->new($dict_file_abs_path, 'r');
    while (my $line = $fh_dict->getline) {
        chomp $line;
        if ($line =~ /\A\@SQ\sSN:(.+)\sLN/smx){
            $ref_dict{$1} = 1;
        }
    }
    $fh_dict->close();

    clean_slate('RNA-SeQC');
    $CWD = 'RNA-SeQC';

    my $fh_out = IO::File->new(qq[> $gtf_name.gtf]);
    my $fh_gtf = IO::File->new($annotation_abs_path, 'r');

    my ($valid_entries, $all_entries,
        $no_transcript_id, $invalid_entries,
        $not_in_dict, $no_gencode_format) = (0, 0, 0, 0, 0, 0);

    while (my $line = $fh_gtf->getline) {
        chomp $line;
        if($line =~ /\A\#/smx){
            print $fh_out qq[$line\n];
            next;
        }
        my $valid = 1;
        $all_entries += 1;
        my @columns = split /\t/smx, $line;
        my $chr = $columns[0] // q[];
        my $attr = $columns[8] // q[];
        if (! exists $ref_dict{$chr}) {
            $not_in_dict += 1;
            $valid = 0;
        }
        if ($attr =~ $GENCODE_ATTR_REGEX) {
            if ($attr !~ /transcript_id/ismx){
                $no_transcript_id += 1;
                $valid = 0;
            }
        } else {
            $no_gencode_format += 1;
            $valid = 0;
        }
        if ($valid) {
            $valid_entries += 1;
            print $fh_out qq[$line\n];
        }
    }

    $fh_gtf->close();
    $fh_out->close();

    my $total_discarded = $not_in_dict + $no_gencode_format + $no_transcript_id;
    my $error_rate = sprintf '%.5f', $valid_entries / $all_entries;
    say qq[*** Rate of valid/invalid records: $error_rate:\n].
        qq[****** Total records in source: $all_entries\n].
        qq[****** Total records in target: $valid_entries\n].
        qq[*** Discarded: \n].
        qq[****** Not in reference dict: $not_in_dict\n].
        qq[****** Not in Gencode GTF format: $no_gencode_format\n].
        qq[****** No transcript_id present: $no_transcript_id\n];
    if ($error_rate < $MIN_RATE) {
        carp q[*** RNA-SeQC annotation creation: ].
             q[too many records were discarded (>50%).].
             q[Ignore this message if that was expected.];
    }

    log_timestamp('stop', $tool);
    return 1;
}


sub salmon {
    my ($salmon_version, $genome_abs_path) = @_;

    log_timestamp('start', $salmon_version);
    clean_slate($salmon_version);

    my $annotation_abs_path = get_common_gtf;
    my $genome_name = fileparse($genome_abs_path, $FASTA_FILE_REGEX);
    my $ref_transcriptome_exists = 1;
    my @fasta = inspect_dir('fasta', $FASTA_FILE_REGEX);
    if (! @fasta) {
        say q[*** $salmon_version index creation: ].
            q[Cannot access fasta: no such directory];
        $ref_transcriptome_exists = 0;
    } else {
        if (scalar @fasta == 0) {
            say qq[*** $salmon_version index creation: ].
                qq[no transcriptome reference fasta file present.];
            $ref_transcriptome_exists = 0;
        } elsif (scalar @fasta > 1) {
            carp qq[*** $salmon_version index creation failed: ].
                 qq[One and only one transcriptome reference fasta file expected.];
            log_timestamp('abort', $salmon_version);
            return 0;
        }
    }

    if (! $ref_transcriptome_exists) {
        say qq[*** $salmon_version index creation: ].
            qq[attempting to create a reference transcriptome fasta file.];
        if (! fasta(qq[fasta (from $salmon_version)], $genome_abs_path)){
            carp qq[*** $salmon_version index creation failed: ].
                 qq[create transcriptome reference fasta file.];
            log_timestamp('abort', $salmon_version);
            return 0;
        } else {
             # remove it from the build list in case it's there via CLI
            delete $build{'fasta'};
        }
    }

    $CWD = $salmon_version;
    my $executable = $salmon_version;
    my $salmon_command = qq[$executable --no-version-check index ].
                         qq[-t ../fasta/$genome_name.transcripts.fa ].
                         qq[-i . --type quasi -k 31];

    if(! execute_command($salmon_version, $salmon_command)){
        log_timestamp('abort', $salmon_version);
        return 0;
    }

    log_timestamp('stop', $salmon_version);
    return 1;
}


sub salmon_latest {
    my $tool = qq[salmon (latest: $LATEST_SALMON)];
    log_timestamp('start', $tool);
    if (! create_symlink($LATEST_SALMON, 'salmon', $tool)) {
        log_timestamp('abort', $tool);
        return 0;
    }
    log_timestamp('stop', $tool);
    return 1;
}


sub tophat2 {
    my ($tool, $genome_abs_path) = @_;

    log_timestamp('start', $tool);
    clean_slate('tophat2');

    my $bt2_index_name;
    if (defined $bt2_dir) {
        my @bt2_files = inspect_dir($bt2_dir, $BT2_FILE_REGEX);
        if (! @bt2_files) {
            carp q[*** TopHat2 index creation: ].
                 qq[Cannot access $bt2_dir: no such directory];
            log_timestamp('abort', $tool);
            return 0;
        } else {
            $bt2_index_name = fileparse($bt2_files[0], $BT2_FILE_REGEX);
        }
    } else {
        carp q[*** TopHat2 index creation: ].
             q[Path to bowtie2 index required when building TopHat2 index.];
        log_timestamp('abort', $tool);
        return 0;
    }

    my $annotation_abs_path = get_common_gtf;
    my $genome_name = fileparse($genome_abs_path, $FASTA_FILE_REGEX);
    my $tophat2_command = qq[tophat2 --output-dir tophat2 --GTF=$annotation_abs_path ].
                          qq[--transcriptome-index=tophat2/$genome_name.known ].
                          qq[$bt2_dir/$bt2_index_name];

    if(! execute_command('tophat2', $tophat2_command)){
        log_timestamp('abort', $tool);
        return 0;
    }

    log_timestamp('stop', $tool);
    return 1;
}


__END__


=head1 NAME

Ref_Maker_T - loop through all the transcriptome tools we use to build
an index or an auxiliary file.

=head1 VERSION


=head1 SYNOPSIS

Change CWD to the output ofach tool will be written in the
transcriptome repository. A path to the annotation file in GTF2/GFF3
format - used to build the gtf file that will be used by the other
tools - and a path to the reference genome in fasta format are always
required.

  perl Ref_Maker_T --genome=/path/to/ref/ --annotation=/path/to/annotation/ [options]

If you pass tool names as arguments only the indexes or auxiliary files
of those will be built.

  perl Ref_Maker_T --genome=/path/to/ref/ --annotation=/path/to/annotation/ --salmon --salmon0_8

Some tools need to know the path to specific files to work. E.g.
Tophat2, requires a path where appropriate bowtie2-index files exist.

  perl Ref_Maker_T --genome=/path/to/ref/ --annotation=/path/to/annotation/ --bowtie2=/path/to/bt2/ --tophat2

Under LSF (tl;dr) if running all of the tools, reserve at least 35GB
of memory and 16 cores.

To run under LSF you need to reserve memory as well as multiple cores
for multi-threading. However, resource requirements vary a lot
from tool to tool. If running all of the tools at once, follow the
above recommendation. When running specific tools, those that build
indexes (e.g. Cellranger), requesting memory 12-15 times the size of
the fasta file should be safe. But for tools that build auxiliary
files only (e.g. RNA-SeQC, gtf) then memory requirements are minimal,
4Gb should be enough. E.g. for a do-all run of a mouse genome, syntax
would be:

  bsub -n 16 -M35000 -R'select[mem>35000] rusage[mem=35000] span[hosts=1]' perl Ref_Maker_T --genome=/path/to/ref/

=head1 SUBROUTINES/METHODS

=head2 clean_slate

Delete and previous version/attempt and create a fresh directory.

=head2 log_timestamp

Log a timestamp for the start or end (normal or aborted) of a job.

=head2 create_symlink

Create a symlink between a source and target.

=head1 REQUIRED ARGUMENTS

=head1 OPTIONS

=over

=item B<-annotation=<path to annotation directory>>

Mandatory. Path to the directory where reference genome annotation
is located. Supported formats: GFF3 and GTF2.

=item B<-bowtie2=<path to bowtie2 index directory>>

Required when building TopHat2's known transcripts. Path to the
directory where the bowtie2-generated index of the reference genome
can be found.

=item B<-dictionary=<path to picard dictionary directory>>

Required when building RNA-SeQC's annotation. Path to the directory
where a picard-generated dictionary of the reference genome is
located.

=item B<-genome=<path to ref genome directory>>

Mandatory. Path to the directory where the reference genome file in
fasta format is stored. Suffixes supported: .fa, .fasta.

=item B<-help>

Print this documentation and exit.

=item B<-<tool>> [B<-<tool>>]

Build this <tool>'s files. Accepts: -cellranger, -fasta, -gtf
-rna_seqc, -salmon0_8, -salmon0_10, -salmon, -tophat2. Omit to
execute all of them.

=back

=head1 EXIT STATUS

=head1 DIAGNOSTICS

=head1 DEPENDENCIES

=over

=item strict

=item warnings

=item autodie

=item Carp

=item English

=item Readonly

=item File::Basename

=item File::chdir

=item File::Path

=item File::Copy

=item FindBin

=item Getopt::Long

=item List::Util

=item POSIX

=item npg_tracking::util::abs_path

=item IO::File

=item feature

=item Pod::Usage

=back

=head1 BUGS AND LIMITATIONS

Presence of external programs (executables) is not evaluated.

=head1 DIAGNOSTICS

=head1 INCOMPATIBILITIES

=head1 CONFIGURATION

=head1 AUTHOR

Ruben Bautista

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2018 GRL

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

