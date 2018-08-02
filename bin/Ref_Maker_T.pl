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
Readonly::Scalar my $LATEST_SALMON => "salmon0_8";
Readonly::Scalar my $DEFAULT_ANNOT_DIR => q{gtf};
Readonly::Scalar my $ONE_HUNDRED => 100;
Readonly::Scalar my $MIN_RATE => 0.5;
Readonly::Scalar my $FASTA_FILE_REGEX => qr{[.]f(?:ast|n)?a\Z}imsx;
Readonly::Scalar my $GTF_FILE_REGEX => qr{[.]gtf\Z}imsx;
Readonly::Scalar my $GTF_GFF_FILE_REGEX => qr{[.]g(?:tf|ff\Z)}imsx;
Readonly::Scalar my $BT2_FILE_REGEX => qr{(?<![.]rev)[.][1-9]+[.]bt2\Z}imsx;
Readonly::Scalar my $DICT_FILE_REGEX => qr{[.]dict\Z}imsx;
Readonly::Array  my @GENCODE_GENE_TYPE => qw{gene transcript  exon  CDS
                                             UTR  start_codon  stop_codon};

my ($help, $verbose);
my ($working_dir, @files, $dh, $ref_genome_dir, $bt2_dir, $dict_dir);
my $annotation_dir = $DEFAULT_ANNOT_DIR;
my $rate_valid = $MIN_RATE;
my $save_discarded = 0;
my ($ref_genome_name, $ref_genome_file_abs_path, $bt2_index_name, $ref_dict_namex);


my %build = (
    'rna_seqc'   => 0,
    'tophat2'    => 0,
    'fasta'      => 0,
    'salmon'     => 0,
    'cellranger' => 0,
    'salmon0_10' => 0,
    'salmon0_8'  => 0,
);

if (! @ARGV) { pod2usage(2) };
GetOptions(
    \%build,
    keys %build,
    'annotation=s'           => \$annotation_dir,
    'bowtie2=s'              => \$bt2_dir,
    'dictionary=s'           => \$dict_dir,
    'genome=s'               => \$ref_genome_dir,
    'help'                   => \$help,
    'rate_valid_rnaseqc=s'   => \$rate_valid,
    'save_discarded_rnaseqc' => \$save_discarded,
) || pod2usage(2);
if ($help) {pod2usage(0);}

my $status = main();

exit $status;


sub main {
    my @genome = inspect_dir($ref_genome_dir, $FASTA_FILE_REGEX);
    if (! @genome) {
        croak 'No reference genome directory found.';
    } elsif (scalar @genome != 1){
        croak 'One and only one fasta file expected.';
    }
    $ref_genome_file_abs_path = abs_path(qq{$ref_genome_dir/$genome[0]});
    $ref_genome_name = fileparse($ref_genome_file_abs_path, $FASTA_FILE_REGEX);

    $working_dir = $CWD;

    if ( sum( values %build ) == 0 ) {
        foreach my $k ( keys %build ) { $build{$k} = 1; }
    }

    my (%arguments, %index_command, @failed, @passed, %subs);
    # Build the indexing arguments
    $subs{'fasta'}      = \&fasta;
    $subs{'rna_seqc'}   = \&rna_seqc;
    $subs{'salmon0_8'}  = \&salmon;
    $subs{'salmon0_10'} = \&salmon;
    $subs{'salmon'}     = \&salmon_latest;
    $subs{'cellranger'} = \&cellranger;
    $subs{'tophat2'}    = \&tophat2;

    # Iterate over subs for aligners/tools that
    # require more than a standard command
    foreach my $tool (keys %subs) {
        if ($build{$tool}) {
            if (!$subs{$tool}->($tool)) {
                push @failed, $tool;
            } else {
                push @passed, $tool;
            }
            $CWD = $working_dir;
        }
    }

    if (@failed) {
        my $message = 'WARNING: Task(s) failed, see Ref_Maker_T output for '.
            'details: '.join(', ', @failed)."\n";
        print {*STDERR} $message or carp $OS_ERROR;
        return 1;
    } else {
        return 0;
    }

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
        carp qq{**** symlinking $source failed: No such file or directory or broken link};
        return 0;
    }
    eval {
        (-e $target) && (unlink $target);
        symlink $source, $target;
        1;
    } or do {
        carp qq{**** symlinking $target file failed: $EVAL_ERROR $SKIP_STRING};
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
        carp qq{*** $tool execution failed: $EVAL_ERROR $SKIP_STRING};
        log_timestamp('abort', $tool);
        return 0;
    };
    log_timestamp('stop', $tool);
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


####################
# Tool methods
####################


sub cellranger {
    log_timestamp('start', 'Cellranger');
    clean_slate('10X');

    my @gtf = inspect_dir($annotation_dir, $GTF_FILE_REGEX);
    if (! @gtf) { 
        carp 'Nor gtf or a different annotation directory were found.';
        log_timestamp('abort', 'Cellranger');
        return 0;
    } elsif (scalar @gtf != 1) {
        carp 'One and only one file in gtf format file is expected.';
        log_timestamp('abort', 'Cellranger');
        return 0;
    }
    my $gtf_file_abs_path = abs_path(qq{$annotation_dir/$gtf[0]});
    my $gtf_name = fileparse($gtf_file_abs_path, $GTF_FILE_REGEX);

    say q{*** Cellranger index creation: creating tmp filtered gtf};
    my $cellranger_command = qq{cellranger mkgtf $gtf_file_abs_path }.
        qq{tmp.$gtf_name-filtered.gtf }.
        qq{--attribute=gene_biotype:protein_coding};

    if(execute_command('cellranger - mkgtf', $cellranger_command)){
        say q{*** Cellranger index creation: creating ref index};
        $cellranger_command = qq{cellranger mkref }.
            qq{--genome=10X }.
            qq{--fasta=$ref_genome_file_abs_path }.
            qq{--genes=tmp.$gtf_name-filtered.gtf}.
            qq{--nthreads=16 --memgb=32}; #TODO: this will be documented in the POD but
                                          #      there should be a better way than hard-coded
        if(execute_command('cellranger - mkgtf', $cellranger_command)){
            say q{*** Cellranger index creation: clean up};
            make_path(q{10X/logs}, { mode => $DIR_PERM });
            move(q{Log.out}, q{10X/logs/cellranger_mkref_Log.out});
            move(q{cellranger_mkref.log}, q{10X/logs/});
            move(q{cellranger_mkgtf.log}, q{10X/logs/});
            unlink qq{tmp.$gtf_name-filtered.gtf};
        } else {
            return 0;
        }
    } else {
        return 0;
    };
    log_timestamp('stop', 'Cellranger');

    return 1;
}


sub fasta {
    log_timestamp('start', 'fasta');

    my @gtf_or_gff = inspect_dir($annotation_dir, $GTF_GFF_FILE_REGEX);
    if (! @gtf_or_gff) { 
        carp q{Nor gtf or a different annotation directory were found.};
        log_timestamp('abort', 'fasta');
        return 0;
    }
    my $annot_file_abs_path = abs_path(qq{$annotation_dir/$gtf_or_gff[0]});
    
    clean_slate('fasta');
    $CWD = q{./fasta};

    my $executable = q{gffread};
    my $fasta_command = qq{$executable -w $ref_genome_name.transcripts.fa }.
                        qq{-g $ref_genome_file_abs_path $annot_file_abs_path};
    return execute_command('fasta', $fasta_command);
}


sub rna_seqc{
    # validate annotation file is RNA-SeQC-compliant:
    #- ensembl gtf format assumed (gene type in column 2)
    #- all entries must have transcript id
    #- all entries must be on contig in reference sequence dict.
    #- files must be called *.gtf
    log_timestamp('start', 'RNA-SeQC');

    my @gtf = inspect_dir($annotation_dir, $GTF_FILE_REGEX);
    if (! @gtf) { 
        carp 'Nor gtf or a different annotation directory were found.';
        log_timestamp('abort', 'RNA-SeQC');
        return 0;
    } elsif (scalar @gtf != 1) {
        carp 'One and only one file in gtf format file is expected.';
        log_timestamp('abort', 'RNA-SeQC');
        return 0;
    }
    my $gtf_file_abs_path = abs_path(qq{$annotation_dir/$gtf[0]});
    my $gtf_name = fileparse($gtf_file_abs_path, $GTF_FILE_REGEX);

    my @dict = inspect_dir($dict_dir, $DICT_FILE_REGEX);
    if (! @dict) {
        carp q{Path to reference dictionary file required }.
             q{when creating RNA-SeQC annotation.};
        log_timestamp('abort', 'RNA-SeQC');
        return 0;
    }

    clean_slate('RNA-SeQC');
    $CWD = 'RNA-SeQC';

    my %ref_dict = ();
    my $fh_dict = IO::File->new(qq{$dict_dir/$dict[0]}, 'r');
    while (my $line = $fh_dict->getline) {
        chomp $line;
        if ($line =~ /\A\@SQ\sSN:(.+)\sLN/smx){
            $ref_dict{$1} = 1;
        }
    }
    $fh_dict->close();

    my $fh_out_discarded;
    if ($save_discarded) { 
        $fh_out_discarded = IO::File->new(qq{> $gtf_name.discarded})
    };
    my $fh_out = IO::File->new(qq{> $gtf_name.gtf});
    my $fh_gtf = IO::File->new($gtf_file_abs_path, 'r');

    my ($valid_entries, $all_entries, $no_transcript_id,
        $not_gencode_gene_type, $invalid_entries, $not_in_dict,
        $fewer_cols) = (0, 0, 0, 0, 0, 0, 0);
    while (my $line = $fh_gtf->getline) {
        chomp $line;
        if($line =~ /\A\#/smx){
            print $fh_out qq{$line\n};
            next;
        }
        $all_entries += 1;
        my @columns = split /\t/smx, $line;
        if (scalar @columns < 9) {
            $fewer_cols += 1;
            if ($save_discarded) { print $fh_out_discarded qq{$line\n} };
            next;
        }
        my ($chr, $type, $attributes) = ($columns[0], $columns[2], $columns[8]);
        if (! exists $ref_dict{$chr}) {
            $not_in_dict += 1;
            if ($save_discarded) { print $fh_out_discarded qq{$line\n} };
            next;
        }
        if ($attributes !~ /transcript_id/ismx){
            $no_transcript_id += 1;
            if ($save_discarded) { print $fh_out_discarded qq{$line\n} };
            next;
        }
        if (none {$_ eq $type} @GENCODE_GENE_TYPE) {
            $not_gencode_gene_type += 1;
            # TODO: should these be discarded?
        }
        $valid_entries += 1;
        print $fh_out qq{$line\n};
    }
    $fh_gtf->close();
    $fh_out->close();
    if ($save_discarded) { $fh_out_discarded->close() };
    my $error_rate = sprintf '%.2f', $valid_entries / $all_entries;

    say qq{*** Rate of valid/invalid records: $error_rate:\n}.
        qq{****** With missing columns: $fewer_cols\n}.
        qq{****** Not in ref dict: $not_in_dict\n}.
        qq{****** No transcript_id present: $no_transcript_id\n}.
        qq{****** Gene_type not a Genecode type (not discarded): }.
        qq{$not_gencode_gene_type\n};

    if ($error_rate < $rate_valid) {
        carp qq{*** RNA-SeQC annotation creation: }.
             qq{far too many records discarded.\n}.
             qq{*** Minimum acceptable rate = $rate_valid. }.
             qq{To adjust, use option --rate=<minimum rate>};
        log_timestamp('abort', 'RNA-SeQC');
        return 0;
    }
    log_timestamp('stop', 'RNA-SeQC');

    return 1;
}


sub salmon {
    my $salmon_version = shift;

    log_timestamp('start', $salmon_version);
    clean_slate($salmon_version);

    my $ref_transcriptome_exists = 1;
    my @fasta = inspect_dir('fasta', $FASTA_FILE_REGEX);

    if (! @fasta) {
        say qq{*** $salmon_version index creation: }.
            qq{no fasta directory found.};
        $ref_transcriptome_exists = 0;
    } else {
        if (scalar @fasta == 0) {
            say qq{*** $salmon_version index creation: }.
                qq{no transcriptome reference fasta file present.};
            $ref_transcriptome_exists = 0;
        } elsif (scalar @fasta > 1) {
            carp qq{*** $salmon_version index creation failed: }.
                 qq{One and only one transcriptome reference fasta file expected.};
            log_timestamp('abort', $salmon_version);
            return 0;
        }
    }

    if (! $ref_transcriptome_exists) {
        say qq{*** $salmon_version index creation: }.
            qq{attempting to create a reference transcriptome fasta file.};
        fasta();
        $CWD = $working_dir;
    }

    $CWD = $salmon_version;
    my $executable = $salmon_version;
    my $salmon_command = qq{$executable --no-version-check index }.
                         qq{-t ../fasta/$ref_genome_name.transcripts.fa }.
                         qq{-i . --type quasi -k 31};

    return execute_command($salmon_version, $salmon_command);
}


sub salmon_latest {
    log_timestamp('start', 'salmon (latest)');

    make_path($LATEST_SALMON, { mode => $DIR_PERM });

    if (create_symlink('salmon', $LATEST_SALMON, 'salmon (latest)')) {
        log_timestamp('stop', 'salmon (latest)');
        return 1;
    } else {
        log_timestamp('abort', 'salmon (latest)');
        return 0;
    }
}


sub tophat2 {
    log_timestamp('start', 'TopHat2');

    my @gtf_or_gff = inspect_dir($annotation_dir, $GTF_GFF_FILE_REGEX);
    if (! @gtf_or_gff) { 
        say 'Nor gtf or a different annotation directory were found.';
        log_timestamp('abort', 'TopHat2');
        return 0;
    }
    my $annot_file_abs_path = abs_path(qq{$annotation_dir/$gtf_or_gff[0]});

    my @bt2_files = inspect_dir($bt2_dir, $BT2_FILE_REGEX);
    if (! @bt2_files) {
        say q{*** TopHat2 index creation: }.
            q{Path to bowtie2 index required when building TopHat2 index.};
        log_timestamp('abort', 'TopHat2');
        return 0;
    } else {
        $bt2_index_name = fileparse($bt2_files[0], $BT2_FILE_REGEX);
    }

    clean_slate('tophat2');
    $CWD = './tophat2';

    my $executable = 'tophat2';
    my $tophat2_command = qq{$executable -G $annot_file_abs_path }.
                          qq{--transcriptome-index=$ref_genome_name.known }.
                          qq{$bt2_dir/$bt2_index_name};

    return execute_command($executable, $tophat2_command);
}


__END__


=head1 NAME

Ref_Maker_T - loop through all the transcriptome tools we use to build
an index or an auxiliary file.

=head1 VERSION


=head1 SYNOPSIS

Change cwd to parent directory of 'gtf' subdirectory in the transcriptome
repository. Alternatively, provide a path to the annotation file in
GTF22/GFF3 format. Most tools require a path to the reference genome in
fasta format so this is a required argument.

    C<perl Ref_Maker_T --genome=/path/to/ref/ [options]>
    C<perl Ref_Maker_T --genome=/path/to/ref/ --annotation=/path/to/annotation/ [options]>

If you specify tools as arguments only those index/annotation files
will be built.

    C<perl Ref_Maker_T --genome=/path/to/ref/ --salmon>

Some tools need to know the path to specific files to work. E.g.
Tophat2, needs the path where to find approriate bowtie2-index files.

    C<perl Ref_Maker_T --genome=/path/to/ref/ --bowtie2=/path/to/bt2/ --tophat2>

=head1 DESCRIPTION

Produce various transcriptome index or auxiliary files required by NPG.

=head1 USAGE

Tl;dr: Under LSF reserve at least 35GB of memory and 16 cores.

To run under LSF you need to reserve memory as well as multiple cores for
multi-threading. Requesting memory 5.5 times the size of the fasta
file should be safe. E.g. for a 1Gb genome, the syntax is:

bsub -n 16 -M5500 -R'select[mem>5500] rusage[mem=5500] span[hosts=1]' perl Ref_Maker_T

Most vertebrate genomes should be run on the 'long' queue. Smaller genomes
will be fine on 'normal'.

=head1 SUBROUTINES/METHODS

=head2 clean_slate

Delete and previous version/attempt and create a fresh directory.

=head2 log_timestamp

Log a timestamp for the start or end (normal or aborted) of a job.

=head2 create_symlink

Create a symlink between a source and target.

=head1 REQUIRED ARGUMENTS

=head1 OPTIONS

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

=item FindBin

=item Getopt::Long

=item List::Util

=item POSIX

=item npg_tracking::util::abs_path

=item IO::File

=item feature

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

