use strict;
use warnings;
use Test::More tests => 108;
use Test::Exception;
use Test::Deep;
use Cwd qw(abs_path cwd);
use File::Temp qw(tempdir);
use File::Spec::Functions qw(catfile);
use File::Which qw(which);
use Moose::Meta::Class;

use npg_tracking::util::abs_path qw(abs_path);

use_ok('npg_common::roles::software_location');

package class_with_lazy_jar;
use strict;
use Moose;
with qw/npg_common::roles::software_location/;
has 'jar' => (is   => 'rw',isa  => 'NpgCommonResolvedPathJarFile',
     coerce     => 1,lazy_build => 1,);
sub _build_jar {return 'MyJar.jar';}

package class_with_default_jar;
use strict;
use Moose;
with qw/npg_common::roles::software_location/;
has 'jar' => ( is   => 'rw',isa  => 'NpgCommonResolvedPathJarFile',
     coerce     => 1,default    => 'MyJar.jar',);

package main;

my @tools = qw/bwa0_6 samtools samtools_irods bowtie bwa/;

my $temp_dir = tempdir(CLEANUP => 1);

sub _obj {
    return Moose::Meta::Class->create_anon_class(
                roles => [qw/npg_common::roles::software_location/],
           )->new_object(@_);
}

{
    throws_ok { _obj( bwa_cmd => q[]) }
        qr/missing name of executable/,
        "error when the name of the bwa executable is set to an empty string";

    throws_ok { _obj( bwa_cmd => undef) }
        qr/is not an executable/,
        "error when the name of the bwa executable is explicitly undefined";

    throws_ok { _obj( bwa_cmd => q[bwa bowtie]) }
        qr/no 'bwa bowtie' executable is on the path/,
        "error when the name of the bwa executable is set to a complex string";

    throws_ok { _obj( bwa_cmd => [qw(bwa bowtie)]) }
        qr/is not an executable/,
        "error when the name of the bwa executable is set to an array reference";

    throws_ok { _obj( bwa_cmd => {'bwa'=>1,}) }
        qr/is not an executable/,
        "error when the name of the bwa executable is set to a hash reference";

    throws_ok { _obj( bwa_cmd => _obj()) }
        qr/is not an executable/,
        "error when the name of the bwa executable is set to an object reference";
}

my ($abs_path, $software);

{
    local $ENV{PATH} = qq[$temp_dir];
    lives_ok {$software = _obj()} "can create object";

    foreach my $tool ( @tools ) {

        my $method = "${tool}_cmd";
        
        throws_ok {$software->$method}
            qr/no '$tool' executable is on the path/,
            "error when $tool is not on the path";
        
        if ($tool ne 'bwa') {
          system "echo 'mock $tool' > $temp_dir/$tool";
          chmod 0755, "$temp_dir/$tool";
        } else {
          symlink "bwa0_6", "$temp_dir/bwa" or die 'cannot symlink';
        }
        
        $abs_path = abs_path(catfile($temp_dir, $tool));
        lives_ok {$software->$method} 
            "no error when $tool is on the path and is executable";
        is ($software->$method, $abs_path, "returns correct absolute path to $tool");
        lives_ok {$software = _obj($method => qq[$tool]) }
            "$tool is on the path, no error setting it as '$tool'";
        is ($software->$method, $abs_path, "returns correct absolute path");
        
        if ($tool ne 'bwa') {
            chmod 0644, "$temp_dir/$tool";
            throws_ok { _obj($method => qq[$tool]) } 
              qr/no '$tool' executable is on the path/, 
              "error when $tool does exists on the path but is not executable";
            chmod 0755, "$temp_dir/$tool";
        }
    }
}

{
    foreach my $tool ( grep {$_ ne 'bwa'} @tools ) {

        my $method = "${tool}_cmd";
        $abs_path = catfile( abs_path($temp_dir), qq[$tool]);

        chmod 0755, "$temp_dir/$tool";

        lives_ok {$software = _obj($method => $abs_path) }
        "$tool can be specified with absolute path";
        is ($software->$method(), $abs_path, "returns correct absolute path");


        my $cwd = cwd;
        my $rel_path = (split '/',$temp_dir)[2];
        chdir(qw[/tmp]);
        lives_ok { $software = _obj($method => qq[$rel_path/$tool],) }
        "$tool can be specified with relative path";
        is ($software->$method(), $abs_path, "returns correct absolute path");
        chdir($cwd);

        mkdir "$temp_dir/$tool-0.1.2";
        system "echo 'mock $tool' > $temp_dir/$tool-0.1.2/$tool";
        chmod 0755, "$temp_dir/$tool-0.1.2/$tool";
        system "/bin/ln -s $temp_dir/$tool-0.1.2/$tool $temp_dir/link-$tool";
        symlink "$temp_dir/$tool-0.1.2/$tool", "$temp_dir/link-$tool";

        lives_ok {$software = _obj($method => qq[$temp_dir/link-$tool]) }
            "$tool can be specified with a symlink";
        my $version_pattern = qw[\d+\.\d+\.\d+];
        is ($software->$method(), catfile( $temp_dir, qq[$tool-0.1.2/$tool]),
        "returns correctly expanded absolute path");

        throws_ok {$software = _obj($method => qq[$temp_dir/t/$tool]) }
            qr['$temp_dir/t/$tool' is an invalid path],
            "error when $tool is specified with an invalid path";
    }
}

{
    my $test_java = 't/data/java';
    my $abs_test_java = abs_path($test_java);

    is (_obj(java_cmd => $test_java)->java_cmd, $abs_test_java, 
    'a full path to test java command when a relative path is given');  

    local $ENV{PATH} = join q[:], 't/data', $ENV{PATH};
    is (_obj()->java_cmd, $abs_test_java, 
    'a full path to test java command when it is on the path');   
}

{
    local $ENV{PATH} = join q[:], 't/bin', $ENV{PATH};
    my $bin = abs_path 't/bin';
    lives_ok {$software = _obj(samtools_cmd => q[samtools], 
                               bwa_cmd => q[bwa], )}
        'object created for predicate test with samtools & bwa';
    my %paths_hash = $software->resolved_paths();
    is ($paths_hash{samtools_cmd}, catfile($bin, 'samtools'), 'samtools_cmd defined correctly');
    is ($paths_hash{bwa_cmd}, catfile($bin, 'bwa'), 'bwa_cmd defined correctly');
    ok (!$paths_hash{bowtie_cmd}, 'does not have bowtie_cmd defined');
    ok (!$paths_hash{samtools_irods_cmd}, 'does not have samtools_irods_cmd defined');
    ok (!$paths_hash{java_cmd}, 'does not have java_cmd defined');

    my $h =  {
              bowtie_cmd => q[/bin/false],
              java_cmd => q[/bin/true]
             };
    my %actual = _obj($h)->resolved_paths();
    is_deeply (\%actual, $h, 'resolved fields as set');    
}

{ # testing find jar location
    my $obj;
    `mkdir $temp_dir/jar_path`;
    `touch $temp_dir/jar_path/MyJar.jar`;
    `touch $temp_dir/jar_path/MyOtherJar.jar`;

    local $ENV{CLASSPATH} = q[/tmp];
    lives_ok { $obj = class_with_lazy_jar->new()} 
         q[build lazy object without specified jar accessor];
    throws_ok { $obj->jar( )} 
         qr/no such file on CLASSPATH: MyJar.jar/, 
         q[lazy build of jar fails when jar not on the classpath];    
    throws_ok { $obj->jar( 'MyJar.jar' )} 
         qr/no such file on CLASSPATH: MyJar.jar/, 
         q[setting jar fails when jar not on the classpath];    
    throws_ok { class_with_lazy_jar->new( jar => 'MyJar.jar',) }
         qr/no such file on CLASSPATH: MyJar.jar/,
         q[build fails with jar accessor specified and jar not on the classpath];

    local $ENV{CLASSPATH} = qq[$temp_dir/jar_path];
    lives_ok { $obj->jar( ) } 
         q[lazy build of jar succeeds with jar on the classpath];
    is( $obj->jar(), qq[$temp_dir/jar_path/MyJar.jar], 
         q[correct jar MyJar.jar] ); 
    lives_ok { $obj->jar( 'MyOtherJar.jar') } 
         q[setting other jar succeeds with jar on the classpath];
    is( $obj->jar(), qq[$temp_dir/jar_path/MyOtherJar.jar], 
         q[correct jar MyOtherJar.jar] ); 
    lives_ok { class_with_lazy_jar->new( jar => 'MyJar.jar',)} 
         q[build object with specified jar accesser and jar is on the classpath];

    lives_ok { $software = class_with_default_jar->new() }
         q[default build succeeds with jar on classpath];
    is ($software->jar, qq[$temp_dir/jar_path/MyJar.jar], q[correct abs_path to jar]);

    local $ENV{CLASSPATH} = q[/tmp];
    throws_ok { class_with_default_jar->new()}
         qr/no such file on CLASSPATH: MyJar.jar/, 
         q[build fails with default jar not on the classpath];

    my $cwd = cwd;
    chdir($temp_dir);
    lives_ok { $software = class_with_default_jar->new( 
                                 jar => abs_path(q[jar_path/MyJar.jar]),);
             } q[build succeeds with absolute path to jar];
    is ($software->jar(), qq[$temp_dir/jar_path/MyJar.jar], q[correct absolute path]);

    lives_ok { $software = class_with_default_jar->new( 
                                 jar => q[jar_path/MyJar.jar],);
             } q[build succeeds with relative path to jar];
    is ($software->jar(), qq[$temp_dir/jar_path/MyJar.jar], q[correct absolute path]);

    chdir(qq[$temp_dir/jar_path]);
    throws_ok { class_with_default_jar->new( jar => q[MyJar.jar],) }
               qr[no such file on CLASSPATH: MyJar.jar],
               q[fail with jar by name only and jar in current directory, but not on the classpath];

    local $ENV{CLASSPATH} = qq[$temp_dir/jar_path];
    throws_ok { $software = class_with_default_jar->new( 
                                 jar => q[/tmp/jar_path/MyJar.jar],
                           ); }
              qr[/tmp/jar_path/MyJar.jar' is an invalid path],
              q[fail when invalid path to jar is given];

    chdir($cwd);
}

{
    my $test = _obj();
    throws_ok { $test->current_version(); }
              qr/Tool command required as argument/ms,
              'Require tool command for version method';
    throws_ok { $test->current_version('t/some'); }
              qr/'t\/some' not found/ms,
              'Aligner command for version method should be a file';
    my $bin = abs_path 't/data/aligners/bin';
    is ( $test->current_version("$bin/bwa"),       q{0.5.8c (r1536)},
         'Return version string for bwa' );
    is ( $test->current_version("$bin/maq"),       q{0.7.1},
         'Return version string for maq' );
    is ( $test->current_version("$bin/samtools"),  q{0.1.11 (r851)},
         'Return version string for samtools' );
    is ( $test->current_version("$bin/smalt"),     q{0.4},
         'Return version string for smalt' );
    is ( $test->current_version(abs_path 't/data/aligners/0.010/tool_not_in_bin_with_no_version'),
          q{0.010}, 'Return version from path without a bin element');
    is ( $test->current_version(abs_path 't/data/aligners/0.010/bin/tool_in_bin_with_no_version'),
          q{0.010}, 'Return version from path with a bin element');
    is ( $test->current_version(abs_path 't/data/aligners/bin/tool_with_no_version'),
          undef, 'Undefined if we cannot get a version successfully');
}

package test_class;
use Moose;
with qw/npg_common::roles::software_location/;

has 'something' => (
    is  => 'rw',
    isa => 'Str',
);
sub clone_with_propagation {
    my ($self) = @_;
    return new test_class ($self->resolved_paths);
}
sub clone_without_propagation {
    my ($self) = @_;
    return new test_class ();
}
sub clone_with_propagation_and_something {
    my ($self) = @_;
    return new test_class ($self->resolved_paths, something => 'who knows?');
}
sub clone_without_propagation_and_something {
    my ($self) = @_;
    return new test_class (something => 'god knows!');
}

package main;

local $ENV{PATH} = join q[:], 't/bin', $ENV{PATH};
my $bwa = which('bwa');

my $test = new test_class(bwa_cmd => q[bwa]);

local $ENV{PATH} = q[];

my $clone = $test->clone_with_propagation();

lives_ok {$clone->bwa_cmd} q[bwa_cmd not on the path and inherited];

my $dummy = $test->clone_without_propagation();

throws_ok {$dummy->bwa_cmd} 
          qr/no 'bwa' executable is on the path/,
          q[bwa_cmd not on the path and not inherited];

$clone = $test->clone_with_propagation_and_something();
lives_ok {$clone->bwa_cmd && $clone->something} q[both bwa_cmd and attribute are inherited];

is ($clone->bwa_cmd, abs_path($bwa), q[correct absolute path returned for bwa]);
is ($clone->something, q[who knows?], q[has correct attribute set]);

$clone = $test->clone_without_propagation_and_something();
throws_ok {$clone->bwa_cmd} 
          qr/no 'bwa' executable is on the path/,
          q[bwa_cmd not on the path and not inherited];
is ($clone->something, q[god knows!], q[has correct attribute set]);


1;
