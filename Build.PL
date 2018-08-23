use strict;
use warnings;
use npg_tracking::util::build;

my $requires = {
                'autodie'                         => '2.10',
                'Carp'                            => '1.04',
                'Cwd'                             => 0,
                'Digest::MD5'                     => 0,
                'English'                         => '1.02',
                'Exporter'                        => '5.63',
                'Fatal'                           => 0,
                'Fcntl'                           => 0,
                'File::Basename'                  => 0,
                'File::chdir'                     => 0,
                'File::Compare'                   => 0,
                'File::Copy'                      => 0,
                'File::Path'                      => 0,
                'File::Slurp'                     => '9999.13',
                'File::Spec'                      => 0,
                'File::Spec::Functions'           => 0,
                'File::stat'                      => 0,
                'File::Temp'                      => 0,
                'File::Which'                     => 0,
                'FindBin'                         => '1.47',
                'IO::File'                        => 0,
                'IO::Tee'                         => 0,
                'IPC::Open3'                      => 0,
                'IPC::Cmd'                        => 0,
                'lib'                             => '0.5565',
                'List::Util'                      => '1.21',
                'Math::Round'                     => '0.06',
                'MIME::Lite'                      => '3.024',
                'Moose'                           => '0.93',
                'Moose::Meta::Class'              => '0',
                'Moose::Role'                     => '0.93',
                'MooseX::Role::Parameterized'     => '0',
                'Moose::Util::TypeConstraints'    => 0,
                'MooseX::ClassAttribute'          => 0,
                'MooseX::Storage'                 => 0,
                'Parallel::ForkManager'           => 0,
                'Perl6::Slurp'                    => '0.03',
                'Pod::Usage'                      => '0',
                'POSIX'                           => '1.09',
                'Readonly'                        => '1.03',
                'Readonly::XS'                    => 0,
                'strict'                          => '1.03', 
                'warnings'                        => '1.05',

                'npg_qc::Schema'                            => 0,
                'npg_tracking::data::reference::list'       => 0,
                'npg_tracking::util::abs_path'              => 0,         
};

my $builder =  npg_tracking::util::build->new(

          'module_name'         => 'npg_seq_common',
          'dist_author'         => q(wtsi-npg <npg@sanger.ac.uk>),
          'dist_abstract'       =>
  'General purpose scripts and modules used by WTSI production pipeline and manual data processing',
          'license'             => 'gpl',
          'dist_version'        => npg_tracking::util::build->git_tag(),

          'configure_requires' => {
                'ExtUtils::CBuilder'              => 0,
                'npg_tracking::util::build'       => 0,
          },

          'build_requires' => {
                'ExtUtils::CBuilder'              => 0,
                'npg_tracking::util::build'       => 0,
                'Test::Cmd'                       => 0,
                'Test::Compile'                   => 0,
                'Test::Distribution'              => 0,
                'Test::Deep'                      => '0.103',
                'Test::Exception'                 => '0.27',
                'Test::More'                      => '0.86',
                'Test::Pod'                       => 0,
                'Test::Pod::Coverage'             => 0,
                'Test::Perl::Critic'              => 0,
          },

          'requires'     => $requires,

          'dist'         => { COMPRESS => 'gzip', SUFFIX => 'gz', },
);

$builder->create_build_script();

1;
