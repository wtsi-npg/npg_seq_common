package npg_common::roles::software_location;

use Moose::Role;
use Moose::Util::TypeConstraints;
use Carp;
use File::Spec::Functions qw(catfile);
use File::Which qw(which);
use IPC::Open3;
use Perl6::Slurp;
use Readonly;

use npg_tracking::util::abs_path qw(abs_path);

our $VERSION = '0';

Readonly::Array my @TOOLS => qw/bwa bwa0_6 samtools samtools_irods bowtie java star/;

subtype 'NpgCommonResolvedPathExecutable'
      => where { (abs_path($_) eq $_) && ( -x ) },
      => as 'Str',
      => message { ($_ || q[]). ' is not an executable' };
coerce 'NpgCommonResolvedPathExecutable',
      from 'Str',
      via { /\//sxm ? (abs_path($_) || croak "'$_' is an invalid path")
                    : ! $_ ? croak 'missing name of executable'
                    : which($_) ? abs_path( (which($_))[0] )
                    : croak "no '$_' executable is on the path" };

foreach my $tool ( @TOOLS ) {
    my $attribute_name = qq[${tool}_cmd];
    has $attribute_name     => (
       is                   => 'ro',
       isa                  => 'NpgCommonResolvedPathExecutable',
       lazy_build           => 1,
       coerce               => 1,
       documentation        => qq[${tool} command, returned resolved to an absolute path to an executable],
    );
   my $build_method = qq[_build_${attribute_name}];
   ##no critic (ProhibitNoStrict ProhibitNoWarnings)
   no strict 'refs';
   no warnings 'redefine';
   *{$build_method}= sub{ return $tool; };
}

subtype 'NpgCommonResolvedPathJarFile'
      => where { ( -r ) && (abs_path($_) eq $_) },
      => as 'Str';
coerce 'NpgCommonResolvedPathJarFile',
      from 'Str',
      via {/\//sxm ? ( abs_path($_) || croak "'$_' is an invalid path" )
                   : _find_jar($_)  || croak "no such file on CLASSPATH: $_"};

sub _find_jar {
    my $name = shift;
    my $jar_path = $ENV{CLASSPATH} || croak qq[Can't find '$_' because CLASSPATH is not set];
    my @search_path = split /\:/smx, $jar_path;
    foreach my $directory (@search_path) {
        my $jar = catfile($directory, $name);
        return abs_path($jar) if (-e $jar);
    }
    return;
}

sub resolved_paths {
    my $self = shift;
    my $predicate_hash = {};
    foreach my $tool (@TOOLS) {
        my $accessor = qq[${tool}_cmd];
        my $method = qq[has_$accessor];
        if ($self->$method) {
            $predicate_hash->{$accessor} = $self->$accessor;
        }
    }
    return %{$predicate_hash};
}

sub current_version {
    my ( $self, $cmd ) = @_;

    croak 'Tool command required as argument' if !$cmd;
    croak "'$cmd' not found" if !-e $cmd;
    my $version;
    if ($cmd =~ /[.]jar$/smx) {
        $cmd = join q[ ], $self->java_cmd, q[-Xmx64m], q[-jar], $cmd, q[--version];
        $version = _get_jar_version($cmd);
    } else {
        my $regex = qr{^(?: $cmd )? \s*
                       version [:]? \s+
                       (\S+ (?: [ \t]+ \S+ )? )
                      }imsx;
        foreach my $arg ( q{}, '--version', '-v', '-V', 'version' ) {
            my $out;
            my $pid = open3( undef, $out, $out, "$cmd $arg" );
            waitpid $pid, 0;
            my $output = slurp($out);
            ($version) = $output =~ m/$regex/igmsx;
            last if defined $version;
        }
    }
    return $version;
}

sub _get_jar_version {
    my $cmd = shift;
    my $out;
    my $pid = open3( undef, $out, $out, $cmd);
    waitpid $pid, 0;
    my $version = slurp($out);
    ##no critic (ErrorHandling::RequireCarping)
    warn qq[Version string for command '$cmd': $version\b];
    ##use critic
    if ($version !~ /^\d+/smx) {
        return;
    } else {
        $version =~ s/\s$//gsmx;
    }
    return $version;
}

no Moose::Util::TypeConstraints;
no Moose::Role;

1;
__END__


=head1 NAME

npg_common::roles::software_location

=head1 VERSION

=head1 SYNOPSIS

  use Moose;
  with 'npg_common::roles::software_location';

=head1 DESCRIPTION

Heuristic for finding at run time installed third-party tools 

=head1 SUBROUTINES/METHODS

=head2 samtools_cmd

samtools command resolved to an absolute path to an executable;
defaults to "samtools" found on the path
 
=head2 samtools_irods_cmd
 
samtools_irods command resolved to an absolute path to an executable;
defaults to "samtools_irods" found on the path

=head2 bwa_cmd

bwa command resolved to an absolute path to an executable;
defaults to "bwa" found on the path

=head2 bwa0_6_cmd

bwa0_6 resolved to an absolute path to an executable;
defaults to "bwa0_6" found on the path
Represents bwa version 0.6 or above.

=head2 bowtie_cmd

bowtie command resolved to an absolute path to an executable;
defaults to "bowtie" found on the path

=head2 bcftools_cmd

bcftools command resolved to an absolute path to an executable;
defaults to "bcftools" found on the path

=head2 java_cmd

java command resolved to an absolute path

=head2 find_jar

find a named jar on the current jar_path

=head2 resolved_paths

returns a hash with accessors which are set 

=head2 current_version

given a tool command, returns the version of the tool
returns undefined if cannot get the version

  my $version = $obj->current_version(q[mypath/bwa]);

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose::Role

=item Moose::Util::TypeConstraints

=item Carp

=item File::Spec::Functions

=item File::Which

=item IPC::Open3

=item Perl6::Slurp

=item Readonly

=item npg_tracking::util::abs_path

=back

=head1 INCOMPATIBILITIES

=head1 DIAGNOSTICS

=head1 BUGS AND LIMITATIONS

Please contact the author with any found.

=head1 AUTHOR

Eduard J. Zuiderwijk, E<lt>ejz@sanger.ac.ukE<gt>

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2017 Genome Research Ltd

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

