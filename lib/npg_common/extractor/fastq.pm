package npg_common::extractor::fastq;

use strict;
use warnings;
use Carp;
use English qw(-no_match_vars);
use Exporter qw(import);
use IO::Tee;
use Fcntl ':seek';
use Readonly;

our $VERSION = '0';

our @EXPORT_OK = qw(
                    read_count
                    split_reads
                   );

Readonly::Scalar my $LINES_PER_READ => 4;

sub read_count {
    my $fname = shift;

    my $exe = qq[wc -l $fname];
    open my $fh, q[-|], $exe or croak $ERRNO;
    my $line_count = <$fh>;
    close $fh or carp $ERRNO;
    ($line_count) = $line_count =~ /(\d+)/smx;
    if (defined $line_count) {
        $line_count = int $line_count/$LINES_PER_READ;
    } else {
        croak "Failed to get a line count for file $fname";
    }

    return $line_count;
}


sub split_reads { ##no critic (ProhibitExcessComplexity)
    my ($fq, $num_bases, $new_fqs) = @_;

    if (!$fq) {
        croak q[Input file name should be given];
    }
    if (!$num_bases) {
        croak q[Read lengths for the target files should be given as an array reference];
    }

    ## no critic (RequireBriefOpen ProhibitMagicNumbers ProhibitTwoArgOpen ProhibitDeepNests)
    open my $source, q[<], $fq or croak qq[Cannot open $fq for reading.];

    my $num_bases1;
    my $num_bases2;

    if (!@{$num_bases}) {
        my $second_line = <$source>;
        if ($second_line) {
            $second_line = <$source>;

            if ($second_line) {
                my $read_length = (length $second_line) - 1;
                if ($read_length % 2 != 0) {
                    croak "Odd number of bases in $fq that should be split in halves";
                }
                $num_bases1 = $read_length / 2;
            } else {
                croak qq[Only one line in $fq];
            }
        } else {
            $num_bases1 = 1;
	}
        seek $source, 0, SEEK_SET;
        $num_bases2 = $num_bases1;
    } else {
        if ($num_bases->[0] <= 0 || (scalar @{$num_bases} > 1 && $num_bases->[1] <= 0)) {
	    croak q[Target read length should be positive];
        }
        $num_bases1 = int $num_bases->[0];
        $num_bases2 = scalar @{$num_bases} > 1 ? int $num_bases->[1] : 0;
    }

    my $total_wanted = $num_bases1 + $num_bases2;

    my $dest1;
    my $dest2;

    my $fq1_fh;
    my $fq2_fh;
    my @destinations = ();

    my $fqe1 = $new_fqs && @{$new_fqs} ? $new_fqs->[0] : $num_bases1 . q[_1.fastq];
    open $fq1_fh, q[>], $fqe1 or croak qq[Cannot open $fqe1 for writing.];
    push @destinations, $fq1_fh;
    $dest1 = IO::Tee->new(@destinations);

    if ($num_bases2) {
        @destinations = ();
        my $fqe2 = $new_fqs && @{$new_fqs} > 1 ? $new_fqs->[1] : $num_bases2 . q[_2.fastq];
        open $fq2_fh, q[>], $fqe2 or croak qq[Cannot open $fqe2 for writing.];
        push @destinations, $fq2_fh;
        $dest2 = IO::Tee->new(@destinations);
    }

    my $count = 0;
    while (my $line = <$source>) {
        if ($line eq qq[\n]) { next; }
        my $remainder = $count % $LINES_PER_READ;
        if ($remainder == 1 || $remainder == 3) {
	    if (length $line >= $total_wanted + 1) {
                print {$dest1} substr $line, 0, $num_bases1 or croak $ERRNO;
                print {$dest1} qq[\n] or croak $ERRNO;
                if ($num_bases2) {
                    print {$dest2} substr $line, $num_bases1, $num_bases2 or croak $ERRNO;
                    print {$dest2} qq[\n] or croak $ERRNO;
	        }
	    } else {
                my $l = $count+1;
                croak qq[Line number $l in $fq is too short];
	    }
        } elsif ($remainder == 2) {
            print {$dest1} $line or croak $ERRNO;
            if ($num_bases2) {
                print {$dest2} $line or croak $ERRNO;
	    }
	} else {
            if (!$num_bases2) {
                print {$dest1} $line or croak $ERRNO;
	    } else {
                my $l = (length $line) - 1;
                my $s = substr $line, 0, $l;
                print {$dest1} $s . '/1' . "\n" or croak $ERRNO;
                print {$dest2} $s . '/2' . "\n" or croak $ERRNO;
	    }
	}
        $count++;
    }

    close $source or croak $ERRNO;
    close $fq1_fh or croak $ERRNO;
    $dest1->close;

    if (defined $dest2) {
        close $fq2_fh or croak $ERRNO;
        $dest2->close;
    }

    return;
}

1;

__END__

=head1 NAME

npg_common::extractor::fastq

=head1 VERSION

=head1 SYNOPSIS

This module is for extracting parts of fastq files.

=head1 DESCRIPTION

=head1 SUBROUTINES/METHODS

=head2 read_count - returns a number of reads in a fastq file.

  my $count = line_count($my_path);

=head2 split_reads - depending on input, either trims the number of bases to the requested number
or, if two numbers are given, one output file has the trimmed reads and another has the bases that
start after the first part and extend for as long as required. There should be enough bases in a read.
 If a ref to an empty array of length is given as an argument, the
source file is split in two halves.

For a read like this:

@IL14_1008:1:1:470:276
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNANNNNNNNNNANAAANNANNNNNNNNGNANNNNN
+
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<&&*&&&&&&&&&<&<<;&&;&&&&&&&&(&<&&&&&

  split_read("input.fastq", [37]);

gives output

@IL14_1008:1:1:470:276
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<

written to file 37_1.fastq in the current directory

and

  split_read("input.fastq", [37,37]);

gives output

@IL14_1008:1:1:470:276
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<

and

@IL14_1008:1:1:470:276/2
ANNANNNNNNNNNANAAANNANNNNNNNNGNANNNNN
+
<&&*&&&&&&&&&<&<<;&&;&&&&&&&&(&<&&&&&

written to files 37_1.fastq and 37_2.fastq in the current directory.

Files to output can be given as a third argument
  split_read("input.fastq", [37,37], ['out1', 'out2']);

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item strict

=item warnings

=item Carp

=item English

=item Readonly

=item Exporter

=item Fcntl

=item IO::Tee

=back

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

Marina Gourtovaia

=head1 LICENSE AND COPYRIGHT

This program is free software: you can redistribute it and/or modify
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


