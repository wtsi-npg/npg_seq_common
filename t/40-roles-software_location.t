use strict;
use warnings;
use Test::More tests => 101;
use Test::Exception;
use Cwd qw(cwd);
use File::Temp qw(tempdir);
use File::Spec::Functions qw(catfile);
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

{
    my $testdir = "$temp_dir/test_tools";
    mkdir $testdir or die "Failed to create directory $testdir";
    my $file = "$testdir/samtools";
    open my $fh, '>', $file or die "Failed to open $file for writing";
    close $fh;
    chmod 755, $file;
    local $ENV{'PATH'} = join q[:], $testdir, $ENV{'PATH'};

    my $obj = Moose::Meta::Class->create_anon_class(
        roles => ['npg_common::roles::software_location' => { tools => [qw/samtools bowtie/] }],
    )->new_object();
    ok ($obj->can('samtools_cmd'), 'samtools_cmd method exists');
    ok ($obj->can('bowtie_cmd'), 'bowtie_cmd method exists');
    ok (!$obj->can('bwa_cmd'), 'bowtie_cmd method does not exist');

    is ($obj->samtools_cmd(), $file, 'samtools found');

    $obj = Moose::Meta::Class->create_anon_class(
        roles => ['npg_common::roles::software_location' => { tools => [qw/some_test_tool/] }],
    )->new_object();
    ok ($obj->can('some_test_tool_cmd'), 'some_test_tool_cmd method exists');
    throws_ok { $obj->some_test_tool_cmd() } qr/no 'some_test_tool' executable is on the path/,
        'error when the tool is not on the path';

    $file = "$testdir/some_test_tool";
    open $fh, '>', $file or die "Failed to open $file for writing";
    close $fh;
    chmod 755, $file;
    is ($obj->some_test_tool_cmd(), $file, 'tool found');
}

1;
