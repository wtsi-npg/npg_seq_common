package npg_common::roles::software_location;

use Moose::Role;
use MooseX::Role::Parameterized;
use Moose::Util::TypeConstraints;
use Carp;
use File::Spec::Functions qw(catfile splitdir);
use File::Which qw(which);
use IPC::Open3;
use Perl6::Slurp;
use Readonly;

use npg_tracking::util::abs_path qw(abs_path);

our $VERSION = '0';

Readonly::Array my @TOOLS => qw/bwa bwa0_6 samtools samtools_irods bowtie java star minimap2/;

subtype 'NpgCommonResolvedPathExecutable'
      => where { ((abs_path($_) || q[]) eq $_) && ( -x ) },
      => as 'Str',
      => message { ($_ || q[]). ' is not an executable' };
coerce 'NpgCommonResolvedPathExecutable',
      from 'Str',
      via { /\//sxm ? (abs_path($_) || croak "'$_' is an invalid path")
                    : ! $_ ? croak 'missing name of executable'
                    : which($_) ? abs_path( (which($_))[0] )
                    : croak "no '$_' executable is on the path" };

parameter tools => (
      isa      => 'ArrayRef',
      required => 1,
      default  => sub { return [@TOOLS]; },
);

role {
  my $p = shift;

  foreach my $tool ( @{$p->tools} ) {
    my $attribute_name = qq[${tool}_cmd];
    has $attribute_name     => (
       is                   => 'ro',
       isa                  => 'NpgCommonResolvedPathExecutable',
       lazy_build           => 1,
       coerce               => 1,
       documentation        => qq[${tool} command, returned resolved to an absolute path to an executable],
    );
    method qq[_build_${attribute_name}] => sub { return $tool; };
  }
};

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
        if (not defined $version) { # fallback to pulling version from path
            my @path = splitdir $cmd;
            pop @path;
            my $pver = pop @path;
            if ($pver eq q(bin)) { $pver = pop @path; }
            # presume path based version must have a digit followed by a "."
            if ($pver =~ m{[[:digit:]][.]}smxg) { $version = $pver;}
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

Default use

  use Moose;
  with 'npg_common::roles::software_location';

  $self->samtools_cmd(); #OK
  $seld->bwa_cmd();      #OK

Specifying the tools

  use Moose;
  with 'npg_common::roles::software_location' =>
    { tools => [qw/samtools/] };

  $self->samtools_cmd(); #OK
  $seld->bwa_cmd();      #Error, attribute does not exist

  use Moose;
  with 'npg_common::roles::software_location' =>
    { tools => [qw/samtools my_tool/] };

  $self->samtools_cmd(); #OK
  $seld->my_tool_cmd();  #OK  

=head1 DESCRIPTION

Heuristic for finding at run time installed third-party tools.

=head1 SUBROUTINES/METHODS

Attributes for individial tools listed below are available by default.
If an array of tools is specified via the tools parameter, only
the commands for tools in this array are available.

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
defaults to "bwa0_6" found on the path.
Represents bwa version 0.6 or above.

=head2 bowtie_cmd

bowtie command resolved to an absolute path to an executable;
defaults to "bowtie" found on the path.

=head2 java_cmd

java command resolved to an absolute path

=head2 star_cmd

star command resolved to an absolute path to an executable;
defaults to "star" found on the path.

=head2 minimap2_cmd

minimap2 command resolved to an absolute path to an executable;
defaults to "minimap2" found on the path.

=head2 current_version

Given a full path tool command, returns the version of the tool.
Returns undefined if cannot get the version.

  my $version = $obj->current_version(q[mypath/bwa]);

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=over

=item Moose::Role

item MooseX::Role::Parameterized

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

=over

=item Eduard J. Zuiderwijk

=item David K. Jackson

=item Marina Gourtovaia

=back

=head1 LICENSE AND COPYRIGHT

Copyright (C) 2018 Genome Research Ltd

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

