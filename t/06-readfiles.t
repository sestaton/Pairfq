#!/usr/bin/env perl

# this is a test to fix a bug with reading files with dashes
# https://github.com/sestaton/Pairfq/commit/405fea5e6eaad3fcab3700ab3bf251992708005f
use 5.010;
use strict;
use warnings FATAL => 'all';
use IPC::System::Simple qw(capture system);
use File::Temp;
use File::Basename;
use autodie qw(open);
use Test::More tests => 32;

#TODO: Add tests that sequences and IDs are correct between tests
my $fq_data = _build_fq_data();

makepairs_inmemory($fq_data, $fa_data);

#
# methods
#
sub makepairs_inmemory {
    my ($fq_data, $fa_data) = @_;
    my $fpfq = File::Temp->new( TEMPLATE => "pairfq_fq_XXXX",
				DIR      => 't',
				SUFFIX   => ".fastq",
				UNLINK   => 0 );
    
    my $rpfq = File::Temp->new( TEMPLATE => "pairfq_fq_XXXX",
				DIR      => 't',
				SUFFIX   => ".fastq",
				UNLINK   => 0 );
    
    my $fsfq = File::Temp->new( TEMPLATE => "pairfq_fq_XXXX",
				DIR      => 't',
				SUFFIX   => ".fastq",
				UNLINK   => 0 );
    
    my $rsfq = File::Temp->new( TEMPLATE => "pairfq_fq_XXXX",
				DIR      => 't',
				SUFFIX   => ".fastq",
				UNLINK   => 0 );
    
    my @pfq_fqout = capture([0..5],"bin/pairfq makepairs -f $fq_data->[0] -r $fq_data->[1] -fp $fpfq -rp $rpfq -fs $fsfq -rs $rsfq --stats");
    

    for my $fqo (@pfq_fqout) {
	if ($fqo =~ /Total forward reads\s.*(\d+)/) { 
	    is($1, 8, 'Correct number of forward fastq reads calculated in memory');
	}
	elsif ($fqo =~ /Total reverse reads\s.*(\d+)/) {
	    is($1, 6, 'Correct number of reverse fastq reads calculated in memory');
	}
	elsif ($fqo =~ /Total forward paired reads\s.*(\d+)/) {
	    is($1, 6, 'Correct number of paired forward fastq reads in memory');
	}
	elsif ($fqo =~ /Total reverse paired reads\s.*(\d+)/) {
	    is($1, 6, 'Correct number of paired reverse fastq reads in memory');
	}
	elsif ($fqo =~ /Total forward unpaired reads\s.*(\d+)/) {
	    is($1, 2, 'Correct number of unpaired forward fastq reads in memory');
	}
	elsif ($fqo =~ /Total reverse unpaired reads\s.*(\d+)/) {
	    is($1, 0, 'Correct number of unpaired reverse fastq reads in memory');
	}
	elsif ($fqo =~ /Total paired reads\s.*(\s\d+)/) {
	    my $tot = $1; $tot =~ s/^\s//;
	    is($tot, 12, 'Correct number of total paired fastq reads in memory');
	}
	elsif ($fqo =~ /Total unpaired reads\s.*(\d+)/) {
	    is($1, 2, 'Correct number of total upaired fastq reads in memory');
	}
    }
    
    unlink $fpfq, $rpfq, $fsfq, $rsfq;
    
}

sub _build_fq_data {
    my $tmpfq1 = File::Temp->new( TEMPLATE => "pairfq-fq-XXXX",
				  DIR      => 't',
				  SUFFIX   => ".fastq",
				  UNLINK   => 0 );
    
    my $tmpfq1_name = $tmpfq1->filename;
    
    say $tmpfq1 '@HWI-ST765:123:D0TEDACXX:5:1101:2872:2088/1';
    say $tmpfq1 'TTGTCTTCCAGATAATTCCGCTATGTTCAACAAATATGTTAGATTCAAGTTTTTCTTGATAAACCTATTTAAAACCATGAAACTGATTCAATCGATTCAAT';
    say $tmpfq1 '+';
    say $tmpfq1 'CCCFFFFFHHHHGJJJJJJIJJJJJJJJJJJJJJJJIIIIJJJIJJIJJHIJJJJJIJIJJJIJJJJJJJJJIJJIJJHHHHHHFFFFFFEEEDEDDDDED';
    say $tmpfq1 '@HWI-ST765:123:D0TEDACXX:5:1101:6511:2225/1';
    say $tmpfq1 'GGGGTTTGAATTGGAATTGAACAAACTCGTGGGACCCCTAGACTGACCGGCTATATACTCAACCTGCTCTAAAGTAAGTGTGGGACACTCGAGCGTGTCGT';
    say $tmpfq1 '+';
    say $tmpfq1 '@BCFDFFFHHHHHJJJJJJJJIJJJJIJJGIJJIIJJJEIJJJJIJJJJJJIJHHHHHHHFFFFFAEEEEDDDCDCCDCCEDDDDDDDDDDBBDDBDBDB<';
    say $tmpfq1 '@HWI-ST765:123:D0TEDACXX:5:1101:12346:2094/1';
    say $tmpfq1 'CGTTCGTTAATTTAGTTAATTTAATTTAATTGCAAATTTTGGATTTTTAGAAACTCTCCCTCTCAAACATAAAAAATAGTTAGTGTCGCATCAGTGTCAGT';
    say $tmpfq1 '+';
    say $tmpfq1 '@C@FFFFFHFHHHHGHHIJJJJIIJJJIJGGIGHIJHIJJJJIIJJJJIGIIIJEIJJIIIGHIIGIJJJIIIIIHEFFBCFFC>AECDDDDDDCCCDDDD';
    say $tmpfq1 '@HWI-ST765:123:D0TEDACXX:5:1101:16473:2185/1';
    say $tmpfq1 'GCGTGTTGGGCTGACGCCAACCAAATGACGGTGGTTAGGATGGCGGTCAGGTCCTCGACGTTAGCCAATGTGGGCCACCATGTCTCATTGCGAAGTTCAGC';
    say $tmpfq1 '+';
    say $tmpfq1 '@CCDDFDFFHGHFJIGGHGEHGGHGGIEGHI6BF?BFHECBHHGGHDEFCE>;>@CB@BBBBDBDC@CDCDD@B22<?@BD>:C:AA>CCDDB@@8@DCCC';
    say $tmpfq1 '@HWI-ST765:123:D0TEDACXX:5:1101:11717:2411/1';
    say $tmpfq1 'ACACAATGTGCAAGCCAATTAGAAGCCAACTGGACAGCACTGAAGGCTTGGAAAAGTGGCTATAAAAGTTACATAAATAAAGAAGATGTTTTATTTCAAAT';
    say $tmpfq1 '+';
    say $tmpfq1 '@@@DD;DA=C>FB?GFFGECHBHIEFBGBCFGF<BFEI;BGFF4?>BGEIFCCF;;BCC;FFFEE@)=7A;==CA:7?BBD>@CCCCB;;@BBBBDEA>A;';
    say $tmpfq1 '@HWI-ST765:123:D0TEDACXX:5:1101:16191:2473/1';
    say $tmpfq1 'TTTTATCATCTTCCTCTTAGTTTGTTCTCTCTATTTATTCGTGTCCCTTTTTTTTATTTATTGTATTAGCAAACTAAATATCTATATCTAAAATATGGTTA';
    say $tmpfq1 '+';
    say $tmpfq1 'CCCFFFFFHHHHHJJJJJJJHIJJIJIJIGGJJJJJIJJJJGHHIIJJJJJJIJJGIJJJJJIIJHHHBEEEFFFFFEEEFEDEEFEEFEDCDEDEEDACC';
    say $tmpfq1 '@HWI-ST765:123:D0TEDACXX:5:1101:16583:2310/1';
    say $tmpfq1 'CCATCCCCTCTCATCTATCCAAAGCCAACCGTATAATCATGGAACTTGAGAAACAACGCATTCGAGCAAAATATCTCAACAAGAAGTCTATGTTTATGTTT';
    say $tmpfq1 '+';
    say $tmpfq1 'CCCFFFFFHHHHHIGIIIJJJJJJJJJJJHIGHGJIJIJJIJIJIIGIIJGJGHGIEIJGIHGGHFEFFFEEEFFDEDDDDBDDDD>ACADD@DEDDDDED';
    say $tmpfq1 '@HWI-ST765:123:D0TEDACXX:5:1101:17034:2404/1';
    say $tmpfq1 'CATTTGCGGTACTTCACACTAGCATGATGATGAAGGGTGCAATCGTTGCACAAGGGTGCAGTTCCGTTGTATGGCTTTCTAGCAGGGGGTTGAGTTGGTTG';
    say $tmpfq1 '+';
    say $tmpfq1 '@CCFFFFFHHHHHIJJJJJJJJJJJJJJJJJJJJJJJ?FHIJIIJIIIJIJJJJJJ@FHIJEIHHHAHFFFFFEDEEDDDEDDDDDDDD9@DCD@DCD@BB';
 
    my $tmpfq2 = File::Temp->new( TEMPLATE => "pairfq-fq-XXXX",
				  DIR      => 't',
				  SUFFIX   => ".fastq",
				  UNLINK   => 0 );
    
    my $tmpfq2_name = $tmpfq2->filename;

    say $tmpfq2 '@HWI-ST765:123:D0TEDACXX:5:1101:2872:2088/2';
    say $tmpfq2 'ATTGTGTTATAAAGTTTATTTCTATTTCCGTTGCAACTTAAATCTGATTTACATTCATTTTACTTAAACAAACACAATCAAAAGAAACTCAGATCCTACAA';
    say $tmpfq2 '+';
    say $tmpfq2 'CCCFDFFFHHHHHIHIJJJJJJJJJJJJJJJJJJJJJJJIJIIIJJJJJJJJJJJJIJIJJJJJJJIJJJJJJJJJGIJJHHHHHFBEDFDDEEDDDDDDD';
    say $tmpfq2 '@HWI-ST765:123:D0TEDACXX:5:1101:6511:2225/2';
    say $tmpfq2 'AAAGAAGACGGTGACTGAGTGCAATGATTTGTTCGAGAGTTTTGCACATTCTGATATGGACTACAGCACTGCCAGCAGGACTTCCATTCCTGTTACTACCA';
    say $tmpfq2 '+';
    say $tmpfq2 'CCCFFFFFHHHDHIIIGIJHGIJJJJJJJJIIJJIGIIJGIIJJJJJJJJIJIIIJJIJIJJIJJJIJJJIHHHHFFFFDDECEEDDEEFDDDDDDDEDD>';
    say $tmpfq2 '@HWI-ST765:123:D0TEDACXX:5:1101:12346:2094/2';
    say $tmpfq2 'TCAATTAAGTCCAAATAAAGTAATCAATGCAATTGCCAAAGAGTCCGCGGCAACGGCGCCAAAAAACTTGATGTGCTAAAAGTAGTTTAATAAAACAACTA';
    say $tmpfq2 '+';
    say $tmpfq2 '@CCFFFEFHFFFFIIJIIJIIFHIJJGGHIJJIJIJJIJJIIIFHIJIJIIIGJGHGDDDDDDDDDDDDDDDDDCDDDDDDDDDDCDDEDDEDCCCBDDDD';
    say $tmpfq2 '@HWI-ST765:123:D0TEDACXX:5:1101:16473:2185/2';
    say $tmpfq2 'GTTGATTATGTTCTCATGCATACAGGGGTATGGCGATCCCGGACCCAAGTCAGCGACATGGACTCAAGCTTTTAATCGAAGACTACCCGTACGCTTCTGAC';
    say $tmpfq2 '+';
    say $tmpfq2 '@@BFFFFDHHHHDHIJHJIJJJJGHIJIFDCHJCHIGIJJIJJJIIIIIHCEIGHHFFFEECEEDDD>A>ACDDDDDBBABBDDDCDBDDBBDDB@BCDCC';
    say $tmpfq2 '@HWI-ST765:123:D0TEDACXX:5:1101:16583:2310/2';
    say $tmpfq2 'TCCCTTTTATTTATTTTGTTTTTATGAACTTTTGTGATATTGTTGATCACTAGCAGTGGTGTAGCATTGGTGCTATTTGGTACGGTTTACCCTGCACGCGG';
    say $tmpfq2 '+';
    say $tmpfq2 'CCCFFFFFGHHHHIJJJJHIJJJGJJJJIIJJJJFEHIIIIIIIJJJJJJIIJJIIGHJFDDFGGHCHIGFFHIJJJJJIEHGHFEFDECCCE(;ACBDDD';
    say $tmpfq2 '@HWI-ST765:123:D0TEDACXX:5:1101:17034:2404/2';
    say $tmpfq2 'AAAGGTGACAAGAAACCAATCGAAGAATCAAAACCTAAGGATAAACAGACTGAATCCTCCAAGAAGTCAAAGAAGCGGAAGGCTTCTCAGAACTTCACCGT';
    say $tmpfq2 '+';
    say $tmpfq2 'BCCFFDDFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJJJJJJJJIJJJJJJJJGIJJJJJHIHGHHHHFFFDCDDDDDDDDDEDDCDDDDDDDDB';

    return [$tmpfq1_name, $tmpfq2_name];
}


