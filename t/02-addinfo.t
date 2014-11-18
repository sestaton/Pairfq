use 5.010;
use strict;
use warnings FATAL => 'all';
use File::Temp;
use File::Spec;
use autodie qw(open);
use Test::More tests => 6;

my $cmd     = File::Spec->catfile('bin', 'pairfq');
my $fq_data = _build_fq_data();
my $fa_data = _build_fa_data();

my $tmpfq_out = File::Temp->new( TEMPLATE => "pairfq_fq_XXXX",
				 DIR      => 't',
				 SUFFIX   => ".fastq",
				 UNLINK   => 0 );

my $tmpfa_out = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
                                 DIR      => 't',
                                 SUFFIX   => ".fasta",
                                 UNLINK   => 0 );

system("$cmd addinfo -i $fq_data -o $tmpfq_out -p 1 2>&1") == 0 
    or die "system failed: $?";
system("$cmd addinfo -i $fa_data -o $tmpfa_out -p 1 2>&1") == 0 
    or die "system failed: $?";

open my $fq, '<', $tmpfq_out;
open my $fa, '<', $tmpfa_out;

while (my $line = <$fq>) {
    if ($line =~ /^\@HWI/) {
	ok($line =~ qr/\/1$/, 'Can properly add pair information to fastq');
    }
}

while (my $line = <$fa>) {
    if ($line =~ /^\>HWI/) {
        ok($line =~ qr/\/1$/, 'Can properly add pair information to fasta');
    }
}

my $fqo = $tmpfq_out->filename;
my $fao = $tmpfa_out->filename;
unlink $fq_data, $fa_data, $fqo, $fao;

#
# private methods
#
sub _build_fq_data {
    my $tmpfq = File::Temp->new( TEMPLATE => "pairfq_fq_XXXX",
				 DIR      => 't',
				 SUFFIX   => ".fastq",
				 UNLINK   => 0 );
    
    my $tmpfq_name = $tmpfq->filename;
    
    say $tmpfq '@HWI-ST765:123:D0TEDACXX:5:1101:2872:2088';
    say $tmpfq 'TTGTCTTCCAGATAATTCCGCTATGTTCAACAAATATGTTAGATTCAAGTTTTTCTTGATAAACCTATTTAAAACCATGAAACTGATTCAATCGATTCAAT';
    say $tmpfq '+';
    say $tmpfq 'CCCFFFFFHHHHGJJJJJJIJJJJJJJJJJJJJJJJIIIIJJJIJJIJJHIJJJJJIJIJJJIJJJJJJJJJIJJIJJHHHHHHFFFFFFEEEDEDDDDED';
    say $tmpfq '@HWI-ST765:123:D0TEDACXX:5:1101:6511:2225';
    say $tmpfq 'GGGGTTTGAATTGGAATTGAACAAACTCGTGGGACCCCTAGACTGACCGGCTATATACTCAACCTGCTCTAAAGTAAGTGTGGGACACTCGAGCGTGTCGT';
    say $tmpfq '+';
    say $tmpfq '@BCFDFFFHHHHHJJJJJJJJIJJJJIJJGIJJIIJJJEIJJJJIJJJJJJIJHHHHHHHFFFFFAEEEEDDDCDCCDCCEDDDDDDDDDDBBDDBDBDB<';
    say $tmpfq '@HWI-ST765:123:D0TEDACXX:5:1101:12346:2094';
    say $tmpfq 'CGTTCGTTAATTTAGTTAATTTAATTTAATTGCAAATTTTGGATTTTTAGAAACTCTCCCTCTCAAACATAAAAAATAGTTAGTGTCGCATCAGTGTCAGT';
    say $tmpfq '+';
    say $tmpfq '@C@FFFFFHFHHHHGHHIJJJJIIJJJIJGGIGHIJHIJJJJIIJJJJIGIIIJEIJJIIIGHIIGIJJJIIIIIHEFFBCFFC>AECDDDDDDCCCDDDD';
    
    return $tmpfq_name;
}

sub _build_fa_data {
    my $tmpfa = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
                                 DIR      => 't',
                                 SUFFIX   => ".fasta",
                                 UNLINK   => 0 );
    
    my $tmpfa_name = $tmpfa->filename;
    
    say $tmpfa '>HWI-ST765:123:D0TEDACXX:5:1101:2872:2088';
    say $tmpfa 'TTGTCTTCCAGATAATTCCGCTATGTTCAACAAATATGTTAGATTCAAGTTTTTCTTGATAAACCTATTTAAAACCATGAAACTGATTCAATCGATTCAAT';
    say $tmpfa '>HWI-ST765:123:D0TEDACXX:5:1101:6511:2225';
    say $tmpfa 'GGGGTTTGAATTGGAATTGAACAAACTCGTGGGACCCCTAGACTGACCGGCTATATACTCAACCTGCTCTAAAGTAAGTGTGGGACACTCGAGCGTGTCGT';
    say $tmpfa '>HWI-ST765:123:D0TEDACXX:5:1101:12346:2094';
    say $tmpfa 'CGTTCGTTAATTTAGTTAATTTAATTTAATTGCAAATTTTGGATTTTTAGAAACTCTCCCTCTCAAACATAAAAAATAGTTAGTGTCGCATCAGTGTCAGT';
    
    return $tmpfa_name;
}
