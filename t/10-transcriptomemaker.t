use Test::More tests => 107;
use Test::Exception;
use strict;
use warnings;
use autodie;
use Cwd qw/getcwd/;
use IO::File;
use File::Copy;
use File::Path qw/make_path/;
use File::Slurp;
use File::Spec::Functions;
use File::Temp qw/tempdir/;
use Digest::MD5;
use Readonly;
use JSON;


# test the Transcriptome_Maker script by building auxiliary files for E coli
# confirm md5 checksum of expected output files
SKIP: {
    skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
        107 unless ($ENV{'TOOLS_INSTALLED'});

    my $start_dir = getcwd();

    my $reference_dir = 't/data/references/E_coli/O127_H6_str_E2348_69';
    unless (-d $reference_dir) { die qq[Cannot find master directory: $reference_dir\n] };


    my $fasta_master = qq[$reference_dir/fasta/FM180568.fa];
    my $dict_master = qq[$reference_dir/fasta/FM180568.fa.dict];
    my $gff_master = qq[$reference_dir/annotation/E_coli_o127_h6_str_e2348_69.gff3];

    foreach my $file (($fasta_master, $dict_master, $gff_master)) {
        unless (-e $file) { die qq[Cannot find master file: $file\n] };
    }

    my $tmp = tempdir('Transcriptome_Maker_test_XXXXXX', CLEANUP => 1, DIR => '/tmp');
    note qq[Created temporary directory: $tmp];

    my $tmp_reference_dir = qq[$tmp/reference];
    foreach my $sub_dir ((qw/annotation bowtie2 fasta/)) {
        make_path(qq[$tmp_reference_dir/$sub_dir]);
        my @files = glob(qq[$reference_dir/$sub_dir/*]);
        foreach my $file (@files){
            copy($file, qq[$tmp_reference_dir/$sub_dir]);
        }
    }

    local $ENV{'PATH'} = join q[:], join(q[/], $start_dir, 'bin'), $ENV{'PATH'};

    my $tmp_transcriptome_dir = qq[$tmp/transcriptome];
    make_path($tmp_transcriptome_dir);
    chdir($tmp_transcriptome_dir);

    my $transcriptome_maker = qq[$start_dir/bin/Transcriptome_Maker];
    my $transcriptome_maker_log = qq[$tmp/transcriptome_maker.log];

    # Run Transcriptome_Maker, redirecting stdout and stderr to .log file
    is(system(qq[$transcriptome_maker > $transcriptome_maker_log 2>&1]),
        512, 'script fails when given no arguments');

    # Test script's options
    is(system(qq[$transcriptome_maker --genome=/path/to/fasta/ ].
              qq[>> $transcriptome_maker_log.errored 2>&1]),
        256, 'script fails when missing path to annotation');
    is(system(qq[$transcriptome_maker --annotation=/path/to/annotation/ ].
              qq[>> $transcriptome_maker_log.errored 2>&1]),
        256, 'script fails when missing path to reference genome');
    is(system(qq[$transcriptome_maker --annotation=$tmp_reference_dir/annotation ].
              qq[--genome=$tmp_reference_dir/fasta --rna_seqc ].
              qq[>> $transcriptome_maker_log.errored 2>&1]),
        256, 'tool: rna_seqc - script fails when missing path to dictionary');
    is(system(qq[$transcriptome_maker --annotation=$tmp_reference_dir/annotation ].
              qq[--genome=$tmp_reference_dir/fasta --tophat2 ].
              qq[>> $transcriptome_maker_log.errored 2>&1]),
        256, 'tool: tophat2 - script fails when missing path to bowtie2 indexes');

    # test full execution
    is(system(qq[$transcriptome_maker ].
              qq[--annotation=$tmp_reference_dir/annotation ].
              qq[--genome=$tmp_reference_dir/fasta ].
              qq[--dictionary=$tmp_reference_dir/fasta ].
              qq[--bowtie2=$tmp_reference_dir/bowtie2 ].
              qq[> $transcriptome_maker_log 2>&1]),
        0, 'script runs successfully');

    # verify md5 checksum for all other files
    my %expected_md5 = (
        '10X/star/chrStart.txt' => 'c8c7ad562987f74373c534f6d3e3e132',
        '10X/star/exonInfo.tab' => '2d2accee56dd355f02f9b1d896363b79',
        '10X/star/SAindex' => '7ca7e8e2469741cf540e9c51ace5cd6e',
        '10X/star/transcriptInfo.tab' => '03f7ec35743d3876fcbfd615c9f25558',
        '10X/star/chrName.txt' => '8e994d412c75a98e0233688e4cf6d7b1',
        '10X/star/chrNameLength.txt' => 'c3ffbdc3c0adfe4d961f67ab8fc2208d',
        '10X/star/SA' => '064f5be0a780269539d2e2deef6e89be',
        '10X/star/chrLength.txt' => '41ab06e38c19f53ec725387d106ada0a',
        '10X/star/sjdbList.fromGTF.out.tab' => 'd41d8cd98f00b204e9800998ecf8427e',
        '10X/star/exonGeTrInfo.tab' => '3e56e718cc8b1454bb659aa52d363dcc',
        '10X/star/Genome' => '3ca03428bf591c721f387762b5c01be3',
        '10X/star/geneInfo.tab' => 'c58ab7863c15e8263ed288ce6845b7b3',
        '10X/star/sjdbInfo.txt' => '1082ab459363b3f2f7aabcef0979c1ed',
        '10X/star/sjdbList.out.tab' => 'd41d8cd98f00b204e9800998ecf8427e',
        '10X/genes/genes.gtf' => '5b1f48380946eb7315e9703f1f3fd838',
        '10X/fasta/genome.fa' => '63a3f28b06e34cdfe7c8f990961bfc2a',
        'fasta/FM180568.transcripts.fa' => '6851bdd1c6eb8d2c9af0dff519729ce6',
        'gtf/E_coli_o127_h6_str_e2348_69.gtf' => '5f0cc0bc84e8a9aeefa621c9c605487f',
        'salmon0_10/txpInfo.bin' => '85defcc2df47afd73a03f775693468e3',
        'salmon0_10/sa.bin' => '1b70a94fb36dfaed561c9d39ea48a10e',
        'salmon0_10/rsd.bin' => 'fd0fc0e706765673934db17a784a4dc5',
        'salmon0_10/duplicate_clusters.tsv' => '15f17d12eb9b908605bd6dbb1e9ea5c5',
        'salmon0_10/hash.bin' => '446d94b09b80b83ea08ce2c77484accc',
        'salmon0_8/txpInfo.bin' => '85defcc2df47afd73a03f775693468e3',
        'salmon0_8/sa.bin' => '1b70a94fb36dfaed561c9d39ea48a10e',
        'salmon0_8/rsd.bin' => 'fd0fc0e706765673934db17a784a4dc5',
        'salmon0_8/hash.bin' => '6868268b4b76727708b9b95e4e1b801b',
        'tophat2/FM180568.known.2.bt2' => '8a2a6c9b774d955f01faa08d3118d527',
        'tophat2/FM180568.known.4.bt2' => 'e59f3ab141058a7115b9295cb6b48370',
        'tophat2/FM180568.known.fa' => 'e49ee6859733bfd7b8c243514032a6f8',
        'tophat2/FM180568.known.ver' => 'c7467b56b4f9ab91dcb154041eae5219',
        'tophat2/FM180568.known.gff' => '5f0cc0bc84e8a9aeefa621c9c605487f',
        'tophat2/FM180568.known.1.bt2' => 'b0de45ff97afa0e4d223dbf75e212612',
        'tophat2/FM180568.known.rev.1.bt2' => 'bb774f1e7681a3aeb93adfb87efeffa5',
        'tophat2/FM180568.known.rev.2.bt2' => '2676049b65c4c0f1025a59356559feac',
        'tophat2/FM180568.known.fa.tlst' => '7470f9c9e2370035aec4452e62cc72e9',
        'tophat2/FM180568.known.3.bt2' => '4619fea360398efa5d05806e271fea35',
        'RNA-SeQC/E_coli_o127_h6_str_e2348_69.gtf' => '5f0cc0bc84e8a9aeefa621c9c605487f',
    );

    # JSON files expected structures change at runtime
    my %expected_json = (
        '10X/reference.json' => {
            'fasta_hash' => 'adf82a106f6fbb7b2cd5c3af8f8ac756b099fb64',
            'genomes' => ['10X'], 'version' => undef, 'input_fasta_files' => ['FM180568.fa'],
            'gtf_hash' => '507774471004a118beff5e72106d86dab2351d9c',
            'input_gtf_files' => ['tmp.E_coli_o127_h6_str_e2348_69-filtered.gtf'],
            'mem_gb' => 32, 'mkref_version' => '2.1.1', 'threads' => 4,
        },
        'salmon0_10/refInfo.json' => {
            'ReferenceFiles' => ['../fasta/FM180568.transcripts.fa'],
        },
        'salmon0_10/versionInfo.json' => {
            'auxKmerLength' => 31, 'indexType' => 1, 'indexVersion' => 2,
            'hasAuxIndex' => bless( do{\(my $o = 0)}, 'JSON::PP::Boolean' ),
        },
        'salmon0_10/header.json' => {
            'value0' => {'IndexType' => 1, 'IndexVersion' => 'q5', 'KmerLen' => 31,
                'BigSA' => bless( do{\(my $o = 0)}, 'JSON::PP::Boolean' ),
                'UsesKmers' => bless( do{\(my $o = 1)}, 'JSON::PP::Boolean' ),
                'PerfectHash' => bless( do{\(my $o = 0)}, 'JSON::PP::Boolean' ),
                'NameHash' => 'c5f182ab7e454aabd827e10ef0bfb7b944c06254b93b8bd0fdaf661307ebb16b',
                'NameHash512' => '816739f28356b83dbe9dbf77bcb0a85c0f3e0469a88c263b428057c47f64d1ee6cf3a532623aa72e63c9b90a85007a323bd658a8d0ddbd77b4cb066799b29b38',
                'SeqHash' => '1e75b20b0c86e9d02035d83dca9f793d6dc25a09e6bd9ec38ed791264d4c6406',
                'SeqHash512' => 'b6ed02f8151cbc44ac3e4e179d3cbcc4c6171cdd9bf048dc654d8a7551436aa54750baf60fca41f3cb5aac12418edb1f533e3b557c94cd8c5234e4523c600b9d',
            }
        },
        'salmon0_8/refInfo.json' => {
            'ReferenceFiles' => ['../fasta/FM180568.transcripts.fa'],
        },
        'salmon0_8/versionInfo.json' => {
            'auxKmerLength' => 31, 'indexType' => 1, 'indexVersion' => 2,
             'hasAuxIndex' => bless( do{\(my $o = 0)}, 'JSON::PP::Boolean' ),
        },
        'salmon0_8/header.json' => {
            'value0' => {'IndexType' => 1, 'IndexVersion' => 'q5', 'KmerLen' => 31,
                'BigSA' => bless( do{\(my $o = 0)}, 'JSON::PP::Boolean' ),
                'PerfectHash' => bless( do{\(my $o = 0)}, 'JSON::PP::Boolean' ),
                'UsesKmers' => bless( do{\(my $o = 1)}, 'JSON::PP::Boolean' ),
                'NameHash' => 'c5f182ab7e454aabd827e10ef0bfb7b944c06254b93b8bd0fdaf661307ebb16b',
                'SeqHash' => '1e75b20b0c86e9d02035d83dca9f793d6dc25a09e6bd9ec38ed791264d4c6406',
            }
        },
    );

    # log-type files contain dates, working dirs and other changeable
    # contents, so their checksums are hard to verify but at least
    # they exist. Also included: odd per-run binary files
    my %expected_log = (
        '10X/logs/cellranger_mkref_Log.out' => 1,
        '10X/pickle/genes.pickle' => 1,
        '10X/star/genomeParameters.txt' => 1,
        'salmon0_10/quasi_index.log' => 1,
        'salmon0_10/indexing.log' => 1,
        'salmon0_8/quasi_index.log' => 1,
        'salmon0_8/indexing.log' => 1,
        'tophat2/logs/tophat.log' => 1,
        'tophat2/logs/run.log' => 1,
        'tophat2/logs/g2f.out' => 1,
        'tophat2/logs/g2f.err' => 1,
    );

    # still in transcriptome dir...
    # immutable files exist and are the same
    foreach my $path (keys %expected_md5) {
        my $file = catfile($tmp, q[transcriptome], $path);
        ok(-e $file, qq[file $path exists]);
        my $fh = IO::File->new($file, 'r');
        is(Digest::MD5->new->addfile($fh)->hexdigest, $expected_md5{$path}, qq[match MD5 checksum of file $path]);
        $fh->close();
    }

    # JSON files exist and have the same structure
    foreach my $path (keys %expected_json) {
        my $file = catfile($tmp, q[transcriptome], $path);
        ok(-e $file, qq[file $path exists]);
        my $json = from_json(read_file($file));
        my $json_hash = $expected_json{$path};
        is_deeply($json, $json_hash, qq[match contents of JSON file $path]);
    }

    # log-type and other changing odd files just exist
    foreach my $path (keys %expected_log) {
        my $file = catfile($tmp, q[transcriptome], $path);
        ok(-e $file, qq[file $path exists]);
    }

    chdir($tmp);
}

1;