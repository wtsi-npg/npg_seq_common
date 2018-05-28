
use warnings;
use Test::More tests => 84;
use Cwd qw/getcwd/;
use File::Temp qw/tempdir/;
use File::Slurp;
use Digest::MD5;
use JSON;

# test the Ref_Maker script by building references for E coli
# confirm md5 checksum of expected output files

SKIP: {
  skip 'Third party bioinformatics tools required. Set TOOLS_INSTALLED to true to run.',
    84 unless ($ENV{'TOOLS_INSTALLED'});
  my $startDir = getcwd();
  my $fastaMaster = 't/data/references/E_coli/K12/fasta/E-coli-K12.fa';
  unless (-e $fastaMaster) {
    die "Cannot find FASTA master file $fastaMaster\n";
  }
  my $tmp = tempdir('Ref_Maker_test_XXXXXX', CLEANUP => 1, DIR => '/tmp');
  note "Created temporary directory: $tmp";
  my $tmpFasta = $tmp."/fasta";
  mkdir($tmpFasta);
  system("cp $fastaMaster $tmpFasta");
  local $ENV{'PATH'} = join q[:], join(q[/], $startDir, 'bin'), $ENV{'PATH'};

  chdir($tmp);

  # Run Ref_Maker, redirecting stdout and stderr to .log file
  is(system("$startDir/bin/Ref_Maker > Ref_Maker.log 2>&1"), 0, 'Ref_Maker exit status');

  # can't use checksum on Picard .dict, as it contains full path to fasta file
  my $picard = "picard/E-coli-K12.fa.dict";
  ok(-e $picard, "Picard .dict file exists");
  my $picardsl = "fasta/E-coli-K12.fa.dict";
  ok(-e $picardsl, "Picard .dict symlink exists");

  # now verify md5 checksum for all other files
  my %expectedMD5 = (
    'blat/E-coli-K12.fa.2bit' => 'd40176801d2f23f76f7c575843350923',
    'bowtie/E-coli-K12.fa.1.ebwt' => '3c990c336037da8dcd5b1e7794c3d9de',
    'bowtie/E-coli-K12.fa.2.ebwt' => 'de2a7524129643b72c0b9c12289c0ec2',
    'bowtie/E-coli-K12.fa.3.ebwt' => 'be250db6550b5e06c6d7c36beeb11707',
    'bowtie/E-coli-K12.fa.4.ebwt' => 'b5a28fd5c0e83d467e6eadb971b3a913',
    'bowtie/E-coli-K12.fa.rev.1.ebwt' => '65c083971ad3b8a8c0324b80c4398c3c',
    'bowtie/E-coli-K12.fa.rev.2.ebwt' => 'cead6529b4534fd0e0faf09d69ff8661',
    'bowtie2/E-coli-K12.fa.1.bt2' => '757da19e3e1425b223004881d61efa48',
    'bowtie2/E-coli-K12.fa.2.bt2' => 'aa8c2b1e74071eb0296fc832e33f5094',
    'bowtie2/E-coli-K12.fa.3.bt2' => 'be250db6550b5e06c6d7c36beeb11707',
    'bowtie2/E-coli-K12.fa.4.bt2' => 'b5a28fd5c0e83d467e6eadb971b3a913',
    'bowtie2/E-coli-K12.fa.rev.1.bt2' => '8c9502dfff924d4dac0b33df0d20b07e',
    'bowtie2/E-coli-K12.fa.rev.2.bt2' => '5a3d15836114aa132267808e4b281066',
    'bwa0_6/E-coli-K12.fa.amb' => 'fd2be0b3b8f7e2702450a3c9dc1a5d93',
    'bwa0_6/E-coli-K12.fa.ann' => '84365967cebedbee51467604ae27a1f9',
    'bwa0_6/E-coli-K12.fa.bwt' => '09f551b8f730df82221bcb6ed8eea724',
    'bwa0_6/E-coli-K12.fa.pac' => 'ca740caf5ee4feff8a77d456ad349c23',
    'bwa0_6/E-coli-K12.fa.sa' => '6e5b71027ce8766ce5e2eea08d1da0ec',
    'fasta/E-coli-K12.fa' => '7285062348a4cb07a23fcd3b44ffcf5d',
    'fasta/E-coli-K12.fa.fai' => '3bfb02378761ec6fe2b57e7dc99bd2b5',
    'minimap2/E-coli-K12.fa.mmi' => '57c2841e5d5507e77a71d27d99937be0',
    'samtools/E-coli-K12.fa.fai' => '3bfb02378761ec6fe2b57e7dc99bd2b5',
    'star/chrLength.txt' => '9be57bfb0f37bd0c34ce8c3b58c62f0e',
    'star/chrNameLength.txt' => 'f67cb9999741fd15db9fa7a111c9008f',
    'star/chrName.txt' => '3f941ce91fb7e8eca71405051a1cd7c7',
    'star/chrStart.txt' => 'faf5c55020c99eceeef3e34188ac0d2f',
    'star/Genome' => '4b797f9f5143695d6215711ecad3b609',
    'star/genomeParameters.txt' => 'f27946102070db744f6f8922d1202b59',
    'star/SA' => '7fbb7e96f3037de8f575576c72ecb67c',
    'star/SAindex' => '60c5b0bb5dbd4692facd7b58318f8f98',
  );

  ok (-e 'npgqc/E-coli-K12.fa.json', 'json file exists');

  my $json_hash = {
    'reference_path'=>'fasta/E-coli-K12.fa',
    '_summary'=>{'ref_length'=>4639675,
                 'counts'=>{'A'=>1142228,'T'=>1140970,'C'=>1179554,'G'=>1176923}}};
  my $json = from_json(read_file('npgqc/E-coli-K12.fa.json'));
  delete $json->{__CLASS__};
  is_deeply($json,$json_hash,'Compare the JSON file');

  chdir($startDir);
  foreach my $path (keys %expectedMD5) {
    my $file = join q[/], $tmp, $path;
    ok(-e $file, "file $file exists");
    open my $fh, "<", $file || die "Cannot open $file for reading";
    is(Digest::MD5->new->addfile($fh)->hexdigest, $expectedMD5{$path},
      "$path MD5 checksum");
    close $fh;
  }

  chdir($tmp);

  SKIP: {
    skip '10X longranger pipeline is not installed in travis, skip testing!',
      19 if ($ENV{'TRAVIS'});

    # Run Ref_Maker --longranger, redirecting stdout and stderr to .log file
    is(system("$startDir/bin/Ref_Maker --longranger > Ref_Maker_longranger.log 2>&1"), 0, 'Ref_Maker_longranger exit status');

    # now verify md5 checksum for all other files
    my %expectedLongrangerMD5 = (
      '10X/fasta/genome.fa' => '7285062348a4cb07a23fcd3b44ffcf5d',
      '10X/fasta/genome.fa.amb' => 'fd2be0b3b8f7e2702450a3c9dc1a5d93',
      '10X/fasta/genome.fa.ann' => '84365967cebedbee51467604ae27a1f9',
      '10X/fasta/genome.fa.bwt' => '09f551b8f730df82221bcb6ed8eea724',
      '10X/fasta/genome.fa.fai' => '3bfb02378761ec6fe2b57e7dc99bd2b5',
      '10X/fasta/genome.fa.flat' => '05dc7a37701cdc6bcf154344a227983d',
      '10X/fasta/genome.fa.gdx' => '8d41ec62e1b566f03b3b4a8f240d20e6',
      '10X/fasta/genome.fa.pac' => 'ca740caf5ee4feff8a77d456ad349c23',
      '10X/fasta/genome.fa.sa' => '6e5b71027ce8766ce5e2eea08d1da0ec',
    );
  
    chdir($startDir);
    foreach my $path (keys %expectedLongrangerMD5) {
      my $file = join q[/], $tmp, $path;
      ok(-e $file, "file $file exists");
      open my $fh, "<", $file || die "Cannot open $file for reading";
      is(Digest::MD5->new->addfile($fh)->hexdigest, $expectedLongrangerMD5{$path},
        "$path MD5 checksum");
      close $fh;
    }
  }
} # end SKIP no tool installed
1;
