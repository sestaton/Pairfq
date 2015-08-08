use 5.010;
use strict;
use warnings FATAL => 'all';
use File::Temp;
use File::Spec;
use File::Basename;
use autodie qw(open);
use Test::More tests => 22;

my $cmd     = File::Spec->catfile('blib', 'bin', 'pairfq');
my $fq_data = _build_fq_data();
my $fa_data = _build_fa_data();

my $tmpfqfp_out = File::Temp->new( TEMPLATE => "pairfq_fqfp_XXXX",
				  DIR      => 't',
				  SUFFIX   => ".fastq",
				  UNLINK   => 0 );

my $tmpfqfs_out = File::Temp->new( TEMPLATE => "pairfq_fqfs_XXXX",
				  DIR      => 't',
				  SUFFIX   => ".fastq",
				  UNLINK   => 0 );

my $tmpfqrp_out = File::Temp->new( TEMPLATE => "pairfq_fqrp_XXXX",
				  DIR      => 't',
				  SUFFIX   => ".fastq",
				  UNLINK   => 0 );

my $tmpfqrs_out = File::Temp->new( TEMPLATE => "pairfq_fqrs_XXXX",
				  DIR      => 't',
				  SUFFIX   => ".fastq",
				  UNLINK   => 0 );

my $tmpfafp_out = File::Temp->new( TEMPLATE => "pairfq_fafp_XXXX",
				  DIR      => 't',
				  SUFFIX   => ".fasta",
				  UNLINK   => 0 );

my $tmpfafs_out = File::Temp->new( TEMPLATE => "pairfq_fafs_XXXX",
				  DIR      => 't',
				  SUFFIX   => ".fasta",
				  UNLINK   => 0 );

my $tmpfarp_out = File::Temp->new( TEMPLATE => "pairfq_farp_XXXX",
				  DIR      => 't',
				  SUFFIX   => ".fasta",
				  UNLINK   => 0 );

my $tmpfars_out = File::Temp->new( TEMPLATE => "pairfq_fars_XXXX",
				  DIR      => 't',
				  SUFFIX   => ".fasta",
				  UNLINK   => 0 );
my @pfq_fqout = qx($cmd makepairs -i $fq_data -fp $tmpfqfp_out -rp $tmpfqrp_out -fs $tmpfqfs_out -rs $tmpfqrs_out --stats);
system("$cmd makepairs -i $fa_data -fp $tmpfafp_out -rp $tmpfarp_out -fs $tmpfafs_out -rs $tmpfars_out") == 0 
    or die "system failed: $?";

my $tmpfqfp = $tmpfqfp_out->filename;
my $tmpfqrp = $tmpfqrp_out->filename;
my $tmpfqfs = $tmpfqfs_out->filename;
my $tmpfqrs = $tmpfqrs_out->filename;
my $tmpfafp = $tmpfafp_out->filename;
my $tmpfarp = $tmpfarp_out->filename;
my $tmpfafs = $tmpfafs_out->filename;
my $tmpfars = $tmpfars_out->filename;

open my $fqfp, '<', $tmpfqfp;
open my $fqrp, '<', $tmpfqrp;
open my $fqfs, '<', $tmpfqfs;
open my $fqrs, '<', $tmpfqrs;
open my $fafp, '<', $tmpfafp;
open my $farp, '<', $tmpfarp;
open my $fafs, '<', $tmpfafs;
open my $fars, '<', $tmpfars;

my ($fqfp_ct, $fqrp_ct, $fqfs_ct, $fqrs_ct, 
    $fafp_ct, $farp_ct, $fafs_ct, $fars_ct) = (0, 0, 0, 0, 0, 0, 0, 0);

#TODO: add method to return these counts instead of repeating this code
while (<$fqfp>) {
    s/^\s+|\s+$//g;
    $fqfp_ct++ if /\/1$/;
}
close $fqfp;

while (<$fqrp>) {
    s/^\s+|\s+$//g;
    $fqrp_ct++ if /\/2$/;
}
close $fqrp;

while (<$fqfs>) {
    s/^\s+|\s+$//g;
    $fqfs_ct++ if /\/1$/;
}
close $fqfs;

while (<$fqrs>) {
    s/^\s+|\s+$//g;
    $fqrs_ct++ if /\/2$/;
}
close $fqrs;

is($fqfp_ct,              4, 'Correct number of forward paired fastq reads');
is($fqrp_ct,              4, 'Correct number of reverse paired fastq reads');
is($fqfs_ct,              1, 'Correct number of forward singleton fastq reads');
is($fqrs_ct,              1, 'Correct number of reverse singleton fastq reads');
is($fqfs_ct + $fqrs_ct,   2, 'Correct number of total singleton fastq reads');
is($fqfp_ct + $fqrp_ct,   8, 'Correct number of total paired fastq reads');
is($fqfp_ct + $fqrp_ct + 
   $fqfs_ct + $fqrs_ct,  10, 'Correct number of total fastq reads');

for my $fqo (@pfq_fqout) {
    if ($fqo =~ /Total forward reads\s.*(\d+)/) { 
	is($1, 5, 'Correct number of forward fastq reads');
    }
    elsif ($fqo =~ /Total reverse reads\s.*(\d+)/) {
	is($1, 5, 'Correct number of reverse fastq reads');
    }
    elsif ($fqo =~ /Total forward paired reads\s.*(\d+)/) {
	is($1, 4, 'Correct number of paired forward fastq reads');
    }
    elsif ($fqo =~ /Total reverse paired reads\s.*(\d+)/) {
	is($1, 4, 'Correct number of paired reverse fastq reads');
    }
    elsif ($fqo =~ /Total forward unpaired reads\s.*(\d+)/) {
	is($1, 1, 'Correct number of unpaired forward fastq reads');
    }
    elsif ($fqo =~ /Total reverse unpaired reads\s.*(\d+)/) {
	is($1, 1, 'Correct number of unpaired reverse fastq reads');
    }
    elsif ($fqo =~ /Total paired reads\s.*(\s\d+)/) {
	my $tot = $1; $tot =~ s/^\s//;
	is($tot, 8, 'Correct number of total paired fastq reads');
    }
    elsif ($fqo =~ /Total unpaired reads\s.*(\d+)/) {
	is($1, 2, 'Correct number of total upaired fastq reads');
    }
}

unlink $fq_data, $tmpfqfp, $tmpfqrp, $tmpfqfs, $tmpfqrs;

while (<$fafp>) {
    s/^\s+|\s+$//g;
    $fafp_ct++ if /\/1$/;
} 
close $fafp;

while (<$farp>) {
    s/^\s+|\s+$//g;
    $farp_ct++ if /\/2$/;
}
close $farp;

while (<$fafs>) {
    s/^\s+|\s+$//g;
    $fafs_ct++ if /\/1$/;
} 
close $fafs;

while (<$fars>) {
    s/^\s+|\s+$//g;
    $fars_ct++ if /\/2$/;
}
close $fars;

is($fafp_ct,              4, 'Correct number of forward paired fasta reads');
is($farp_ct,              4, 'Correct number of reverse paired fasta reads');
is($fafs_ct,              1, 'Correct number of forward singleton fasta reads');
is($fars_ct,              1, 'Correct number of reverse singleton fasta reads');
is($fafs_ct + $fars_ct,   2, 'Correct number of total singleton fasta reads');
is($fafp_ct + $farp_ct,   8, 'Correct number of total paired fasta reads');
is($fafp_ct + $farp_ct +
   $fafs_ct + $fars_ct,  10, 'Correct number of total fasta reads');

unlink $fa_data, $tmpfafp, $tmpfarp, $tmpfafs, $tmpfars;
#
# private methods
#
sub _build_fq_data {
    my $tmpfq = File::Temp->new( TEMPLATE => "pairfq_fq_XXXX",
				 DIR      => 't',
				 SUFFIX   => ".fastq",
				 UNLINK   => 0 );
    
    my $tmpfq_name = $tmpfq->filename;
    
    say $tmpfq '@HWI-ST765:123:D0TEDACXX:5:1101:2872:2088/1';
    say $tmpfq 'TTGTCTTCCAGATAATTCCGCTATGTTCAACAAATATGTTAGATTCAAGTTTTTCTTGATAAACCTATTTAAAACCATGAAACTGATTCAATCGATTCAAT';
    say $tmpfq '+';
    say $tmpfq 'CCCFFFFFHHHHGJJJJJJIJJJJJJJJJJJJJJJJIIIIJJJIJJIJJHIJJJJJIJIJJJIJJJJJJJJJIJJIJJHHHHHHFFFFFFEEEDEDDDDED';
    say $tmpfq '@HWI-ST765:123:D0TEDACXX:5:1101:2872:2088/2';
    say $tmpfq 'ATTGTGTTATAAAGTTTATTTCTATTTCCGTTGCAACTTAAATCTGATTTACATTCATTTTACTTAAACAAACACAATCAAAAGAAACTCAGATCCTACAA';
    say $tmpfq '+';
    say $tmpfq 'CCCFDFFFHHHHHIHIJJJJJJJJJJJJJJJJJJJJJJJIJIIIJJJJJJJJJJJJIJIJJJJJJJIJJJJJJJJJGIJJHHHHHFBEDFDDEEDDDDDDD';
    say $tmpfq '@HWI-ST765:123:D0TEDACXX:5:1101:6511:2225/1';
    say $tmpfq 'GGGGTTTGAATTGGAATTGAACAAACTCGTGGGACCCCTAGACTGACCGGCTATATACTCAACCTGCTCTAAAGTAAGTGTGGGACACTCGAGCGTGTCGT';
    say $tmpfq '+';
    say $tmpfq '@BCFDFFFHHHHHJJJJJJJJIJJJJIJJGIJJIIJJJEIJJJJIJJJJJJIJHHHHHHHFFFFFAEEEEDDDCDCCDCCEDDDDDDDDDDBBDDBDBDB<';
    say $tmpfq '@HWI-ST765:123:D0TEDACXX:5:1101:12346:2094/1';
    say $tmpfq 'CGTTCGTTAATTTAGTTAATTTAATTTAATTGCAAATTTTGGATTTTTAGAAACTCTCCCTCTCAAACATAAAAAATAGTTAGTGTCGCATCAGTGTCAGT';
    say $tmpfq '+';
    say $tmpfq '@C@FFFFFHFHHHHGHHIJJJJIIJJJIJGGIGHIJHIJJJJIIJJJJIGIIIJEIJJIIIGHIIGIJJJIIIIIHEFFBCFFC>AECDDDDDDCCCDDDD';
    say $tmpfq '@HWI-ST765:123:D0TEDACXX:5:1101:12346:2094/2';
    say $tmpfq 'TCAATTAAGTCCAAATAAAGTAATCAATGCAATTGCCAAAGAGTCCGCGGCAACGGCGCCAAAAAACTTGATGTGCTAAAAGTAGTTTAATAAAACAACTA';
    say $tmpfq '+';
    say $tmpfq '@CCFFFEFHFFFFIIJIIJIIFHIJJGGHIJJIJIJJIJJIIIFHIJIJIIIGJGHGDDDDDDDDDDDDDDDDDCDDDDDDDDDDCDDEDDEDCCCBDDDD';
    say $tmpfq '@HWI-ST765:123:D0TEDACXX:5:1101:16473:2185/2';
    say $tmpfq 'GTTGATTATGTTCTCATGCATACAGGGGTATGGCGATCCCGGACCCAAGTCAGCGACATGGACTCAAGCTTTTAATCGAAGACTACCCGTACGCTTCTGAC';
    say $tmpfq '+';
    say $tmpfq '@@BFFFFDHHHHDHIJHJIJJJJGHIJIFDCHJCHIGIJJIJJJIIIIIHCEIGHHFFFEECEEDDD>A>ACDDDDDBBABBDDDCDBDDBBDDB@BCDCC';
    say $tmpfq '@HWI-ST765:123:D0TEDACXX:5:1101:16583:2310/1';
    say $tmpfq 'CCATCCCCTCTCATCTATCCAAAGCCAACCGTATAATCATGGAACTTGAGAAACAACGCATTCGAGCAAAATATCTCAACAAGAAGTCTATGTTTATGTTT';
    say $tmpfq '+';
    say $tmpfq 'CCCFFFFFHHHHHIGIIIJJJJJJJJJJJHIGHGJIJIJJIJIJIIGIIJGJGHGIEIJGIHGGHFEFFFEEEFFDEDDDDBDDDD>ACADD@DEDDDDED';
    say $tmpfq '@HWI-ST765:123:D0TEDACXX:5:1101:16583:2310/2';
    say $tmpfq 'TCCCTTTTATTTATTTTGTTTTTATGAACTTTTGTGATATTGTTGATCACTAGCAGTGGTGTAGCATTGGTGCTATTTGGTACGGTTTACCCTGCACGCGG';
    say $tmpfq '+';
    say $tmpfq 'CCCFFFFFGHHHHIJJJJHIJJJGJJJJIIJJJJFEHIIIIIIIJJJJJJIIJJIIGHJFDDFGGHCHIGFFHIJJJJJIEHGHFEFDECCCE(;ACBDDD';
    say $tmpfq '@HWI-ST765:123:D0TEDACXX:5:1101:17034:2404/1';
    say $tmpfq 'CATTTGCGGTACTTCACACTAGCATGATGATGAAGGGTGCAATCGTTGCACAAGGGTGCAGTTCCGTTGTATGGCTTTCTAGCAGGGGGTTGAGTTGGTTG';
    say $tmpfq '+';
    say $tmpfq '@CCFFFFFHHHHHIJJJJJJJJJJJJJJJJJJJJJJJ?FHIJIIJIIIJIJJJJJJ@FHIJEIHHHAHFFFFFEDEEDDDEDDDDDDDD9@DCD@DCD@BB';
    say $tmpfq '@HWI-ST765:123:D0TEDACXX:5:1101:17034:2404/2';
    say $tmpfq 'AAAGGTGACAAGAAACCAATCGAAGAATCAAAACCTAAGGATAAACAGACTGAATCCTCCAAGAAGTCAAAGAAGCGGAAGGCTTCTCAGAACTTCACCGT';
    say $tmpfq '+';
    say $tmpfq 'BCCFFDDFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJJJJJJJJIJJJJJJJJGIJJJJJHIHGHHHHFFFDCDDDDDDDDDEDDCDDDDDDDDB';

    return $tmpfq_name;
}

sub _build_fa_data {
    my $tmpfa = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
				 DIR      => 't',
				 SUFFIX   => ".fasta",
				 UNLINK   => 0 );
    
    my $tmpfa_name = $tmpfa->filename;

    say $tmpfa '>HWI-ST765:123:D0TEDACXX:5:1101:2872:2088/2';
    say $tmpfa 'ATTGTGTTATAAAGTTTATTTCTATTTCCGTTGCAACTTAAATCTGATTTACATTCATTTTACTTAAACAAACACAATCAAAAGAAACTCAGATCCTACAA';
    say $tmpfa '>HWI-ST765:123:D0TEDACXX:5:1101:6511:2225/1';
    say $tmpfa 'GGGGTTTGAATTGGAATTGAACAAACTCGTGGGACCCCTAGACTGACCGGCTATATACTCAACCTGCTCTAAAGTAAGTGTGGGACACTCGAGCGTGTCGT';
    say $tmpfa '>HWI-ST765:123:D0TEDACXX:5:1101:6511:2225/2';
    say $tmpfa 'AAAGAAGACGGTGACTGAGTGCAATGATTTGTTCGAGAGTTTTGCACATTCTGATATGGACTACAGCACTGCCAGCAGGACTTCCATTCCTGTTACTACCA';
    say $tmpfa '>HWI-ST765:123:D0TEDACXX:5:1101:12346:2094/1';
    say $tmpfa 'CGTTCGTTAATTTAGTTAATTTAATTTAATTGCAAATTTTGGATTTTTAGAAACTCTCCCTCTCAAACATAAAAAATAGTTAGTGTCGCATCAGTGTCAGT';
    say $tmpfa '>HWI-ST765:123:D0TEDACXX:5:1101:12346:2094/2';
    say $tmpfa 'TCAATTAAGTCCAAATAAAGTAATCAATGCAATTGCCAAAGAGTCCGCGGCAACGGCGCCAAAAAACTTGATGTGCTAAAAGTAGTTTAATAAAACAACTA';
    say $tmpfa '>HWI-ST765:123:D0TEDACXX:5:1101:16473:2185/1';
    say $tmpfa 'GCGTGTTGGGCTGACGCCAACCAAATGACGGTGGTTAGGATGGCGGTCAGGTCCTCGACGTTAGCCAATGTGGGCCACCATGTCTCATTGCGAAGTTCAGC';
    say $tmpfa '>HWI-ST765:123:D0TEDACXX:5:1101:16583:2310/1';
    say $tmpfa 'CCATCCCCTCTCATCTATCCAAAGCCAACCGTATAATCATGGAACTTGAGAAACAACGCATTCGAGCAAAATATCTCAACAAGAAGTCTATGTTTATGTTT';
    say $tmpfa '>HWI-ST765:123:D0TEDACXX:5:1101:16583:2310/2';
    say $tmpfa 'TCCCTTTTATTTATTTTGTTTTTATGAACTTTTGTGATATTGTTGATCACTAGCAGTGGTGTAGCATTGGTGCTATTTGGTACGGTTTACCCTGCACGCGG';
    say $tmpfa '>HWI-ST765:123:D0TEDACXX:5:1101:17034:2404/1';
    say $tmpfa 'CATTTGCGGTACTTCACACTAGCATGATGATGAAGGGTGCAATCGTTGCACAAGGGTGCAGTTCCGTTGTATGGCTTTCTAGCAGGGGGTTGAGTTGGTTG';
    say $tmpfa '>HWI-ST765:123:D0TEDACXX:5:1101:17034:2404/2';
    say $tmpfa 'AAAGGTGACAAGAAACCAATCGAAGAATCAAAACCTAAGGATAAACAGACTGAATCCTCCAAGAAGTCAAAGAAGCGGAAGGCTTCTCAGAACTTCACCGT';
    
    return $tmpfa_name;
}
