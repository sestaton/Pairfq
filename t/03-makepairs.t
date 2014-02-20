#!/usr/bin/env perl

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
my $fa_data = _build_fa_data();

makepairs_inmemory($fq_data, $fa_data);
makepairs_ondisk($fq_data, $fa_data);

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
    
    my @pfq_fqout = capture([0..5],"bin/pairfq makepairs -f $fq_data->[0] -r $fq_data->[1] -fp $fpfq -rp $rpfq -fs $fsfq -rs $rsfq");
    
    my $fpfa = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
				DIR      => 't',
				SUFFIX   => ".fasta",
				UNLINK   => 0 );
    
    my $rpfa = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
				DIR      => 't',
				SUFFIX   => ".fasta",
				UNLINK   => 0 );
    
    my $fsfa = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
				DIR      => 't',
				SUFFIX   => ".fasta",
				UNLINK   => 0 );
    
    my $rsfa = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
				DIR      => 't',
				SUFFIX   => ".fasta",
				UNLINK   => 0 );
    
    my @pfq_faout = capture([0..5],"bin/pairfq makepairs -f $fa_data->[0] -r $fa_data->[1] -fp $fpfa -rp $rpfa -fs $fsfa -rs $rsfa");
    
    for my $fqo (@pfq_fqout) {
	if ($fqo =~ /Total forward reads\s.*(\d+)/) { 
	    is($1, 8, 'Correct number of forward fastq reads calculated');
	}
	elsif ($fqo =~ /Total reverse reads\s.*(\d+)/) {
	    is($1, 6, 'Correct number of reverse fastq reads calculated');
	}
	elsif ($fqo =~ /Total forward paired reads\s.*(\d+)/) {
	    is($1, 6, 'Correct number of paired forward fastq reads');
	}
	elsif ($fqo =~ /Total reverse paired reads\s.*(\d+)/) {
	    is($1, 6, 'Correct number of paired reverse fastq reads');
	}
	elsif ($fqo =~ /Total forward unpaired reads\s.*(\d+)/) {
	    is($1, 2, 'Correct number of unpaired forward fastq reads');
	}
	elsif ($fqo =~ /Total reverse unpaired reads\s.*(\d+)/) {
	    is($1, 0, 'Correct number of unpaired reverse fastq reads');
	}
	elsif ($fqo =~ /Total paired reads\s.*(\s\d+)/) {
	    my $tot = $1; $tot =~ s/^\s//;
	    is($tot, 12, 'Correct number of total paired fastq reads');
	}
	elsif ($fqo =~ /Total upaired reads\s.*(\d+)/) {
	    is($1, 2, 'Correct number of total upaired fastq reads');
	}
    }
    
    unlink $fpfq, $rpfq, $fsfq, $rsfq;
    
    for my $fao (@pfq_faout) {
	if ($fao =~ /Total forward reads\s.*(\d+)/) { 
	    is($1, 8, 'Correct number of forward fasta reads calculated');
	}
	elsif ($fao =~ /Total reverse reads\s.*(\d+)/) {
	    is($1, 6, 'Correct number of reverse fasta reads calculated');
	}
	elsif ($fao =~ /Total forward paired reads\s.*(\d+)/) {
	    is($1, 6, 'Correct number of paired forward fasta reads');
	}
	elsif ($fao =~ /Total reverse paired reads\s.*(\d+)/) {
	    is($1, 6, 'Correct number of paired reverse fasta reads');
	}
	elsif ($fao =~ /Total forward unpaired reads\s.*(\d+)/) {
	    is($1, 2, 'Correct number of unpaired forward fasta reads');
	}
	elsif ($fao =~ /Total reverse unpaired reads\s.*(\d+)/) {
	    is($1, 0, 'Correct number of unpaired reverse fasta reads');
	}
	elsif ($fao =~ /Total paired reads\s.*(\s\d+)/) {
	    my $tot = $1; $tot =~ s/^\s//;
	    is($tot, 12, 'Correct number of total paired fasta reads');
	}
	elsif ($fao =~ /Total upaired reads\s.*(\d+)/) {
	    is($1, 2, 'Correct number of total upaired fasta reads');
	}
    }
    
    unlink $fpfa, $rpfa, $fsfa, $rsfa;   
}

sub makepairs_ondisk {
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
    
    my @pfq_fqout = capture([0..5],"bin/pairfq makepairs -f $fq_data->[0] -r $fq_data->[1] -fp $fpfq -rp $rpfq -fs $fsfq -rs $rsfq -idx");

    my $fpfa = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
                                DIR      => 't',
                                SUFFIX   => ".fasta",
				UNLINK   => 0 );
    
    my $rpfa = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
                                DIR      => 't',
                                SUFFIX   => ".fasta",
                                UNLINK   => 0 );
    
    my $fsfa = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
                                DIR      => 't',
                                SUFFIX   => ".fasta",
                                UNLINK   => 0 );
    
    my $rsfa = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
                                DIR      => 't',
                                SUFFIX   => ".fasta",
                                UNLINK   => 0 );
    
    my @pfq_faout = capture([0..5],"bin/pairfq makepairs -f $fa_data->[0] -r $fa_data->[1] -fp $fpfa -rp $rpfa -fs $fsfa -rs $rsfa -idx");

    for my $fqo (@pfq_fqout) {
        if ($fqo =~ /Total forward reads\s.*(\d+)/) { 
            is($1, 8, 'Correct number of forward fastq reads calculated');
        }
        elsif ($fqo =~ /Total reverse reads\s.*(\d+)/) {
            is($1, 6, 'Correct number of reverse fastq reads calculated');
        }
        elsif ($fqo =~ /Total forward paired reads\s.*(\d+)/) {
            is($1, 6, 'Correct number of paired forward fastq reads');
        }
        elsif ($fqo =~ /Total reverse paired reads\s.*(\d+)/) {
            is($1, 6, 'Correct number of paired reverse fastq reads');
        }
        elsif ($fqo =~ /Total forward unpaired reads\s.*(\d+)/) {
            is($1, 2, 'Correct number of unpaired forward fastq reads');
        }
        elsif ($fqo =~ /Total reverse unpaired reads\s.*(\d+)/) {
            is($1, 0, 'Correct number of unpaired reverse fastq reads');
        }
        elsif ($fqo =~ /Total paired reads\s.*(\s\d+)/) {
            my $tot = $1; $tot =~ s/^\s//;
            is($tot, 12, 'Correct number of total paired fastq reads');
        }
        elsif ($fqo =~ /Total upaired reads\s.*(\d+)/) {
            is($1, 2, 'Correct number of total upaired fastq reads');
        }
    }
    
    unlink $fq_data->[0], $fq_data->[1], $fpfq, $rpfq, $fsfq, $rsfq;

    for my $fao (@pfq_faout) {
        if ($fao =~ /Total forward reads\s.*(\d+)/) { 
            is($1, 8, 'Correct number of forward fasta reads calculated');
        }
        elsif ($fao =~ /Total reverse reads\s.*(\d+)/) {
            is($1, 6, 'Correct number of reverse fasta reads calculated');
        }
        elsif ($fao =~ /Total forward paired reads\s.*(\d+)/) {
            is($1, 6, 'Correct number of paired forward fasta reads');
        }
        elsif ($fao =~ /Total reverse paired reads\s.*(\d+)/) {
            is($1, 6, 'Correct number of paired reverse fasta reads');
        }
        elsif ($fao =~ /Total forward unpaired reads\s.*(\d+)/) {
            is($1, 2, 'Correct number of unpaired forward fasta reads');
        }
        elsif ($fao =~ /Total reverse unpaired reads\s.*(\d+)/) {
            is($1, 0, 'Correct number of unpaired reverse fasta reads');
        }
        elsif ($fao =~ /Total paired reads\s.*(\s\d+)/) {
            my $tot = $1; $tot =~ s/^\s//;
            is($tot, 12, 'Correct number of total paired fasta reads');
        }
        elsif ($fao =~ /Total upaired reads\s.*(\d+)/) {
            is($1, 2, 'Correct number of total upaired fasta reads');
        }
    }
    
    unlink $fa_data->[0], $fa_data->[1], $fpfa, $rpfa, $fsfa, $rsfa;
}

sub _build_fq_data {
    my $tmpfq1 = File::Temp->new( TEMPLATE => "pairfq_fq_XXXX",
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
 
    my $tmpfq2 = File::Temp->new( TEMPLATE => "pairfq_fq_XXXX",
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

sub _build_fa_data {
    my $tmpfa1 = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
				  DIR      => 't',
				  SUFFIX   => ".fasta",
				  UNLINK   => 0 );
    
    my $tmpfa1_name = $tmpfa1->filename;
    
    say $tmpfa1 '>HWI-ST765:123:D0TEDACXX:5:1101:2872:2088/1';
    say $tmpfa1 'TTGTCTTCCAGATAATTCCGCTATGTTCAACAAATATGTTAGATTCAAGTTTTTCTTGATAAACCTATTTAAAACCATGAAACTGATTCAATCGATTCAAT';
    say $tmpfa1 '>HWI-ST765:123:D0TEDACXX:5:1101:6511:2225/1';
    say $tmpfa1 'GGGGTTTGAATTGGAATTGAACAAACTCGTGGGACCCCTAGACTGACCGGCTATATACTCAACCTGCTCTAAAGTAAGTGTGGGACACTCGAGCGTGTCGT';
    say $tmpfa1 '>HWI-ST765:123:D0TEDACXX:5:1101:12346:2094/1';
    say $tmpfa1 'CGTTCGTTAATTTAGTTAATTTAATTTAATTGCAAATTTTGGATTTTTAGAAACTCTCCCTCTCAAACATAAAAAATAGTTAGTGTCGCATCAGTGTCAGT';
    say $tmpfa1 '>HWI-ST765:123:D0TEDACXX:5:1101:16473:2185/1';
    say $tmpfa1 'GCGTGTTGGGCTGACGCCAACCAAATGACGGTGGTTAGGATGGCGGTCAGGTCCTCGACGTTAGCCAATGTGGGCCACCATGTCTCATTGCGAAGTTCAGC';
    say $tmpfa1 '>HWI-ST765:123:D0TEDACXX:5:1101:11717:2411/1';
    say $tmpfa1 'ACACAATGTGCAAGCCAATTAGAAGCCAACTGGACAGCACTGAAGGCTTGGAAAAGTGGCTATAAAAGTTACATAAATAAAGAAGATGTTTTATTTCAAAT';
    say $tmpfa1 '>HWI-ST765:123:D0TEDACXX:5:1101:16191:2473/1';
    say $tmpfa1 'TTTTATCATCTTCCTCTTAGTTTGTTCTCTCTATTTATTCGTGTCCCTTTTTTTTATTTATTGTATTAGCAAACTAAATATCTATATCTAAAATATGGTTA';
    say $tmpfa1 '>HWI-ST765:123:D0TEDACXX:5:1101:16583:2310/1';
    say $tmpfa1 'CCATCCCCTCTCATCTATCCAAAGCCAACCGTATAATCATGGAACTTGAGAAACAACGCATTCGAGCAAAATATCTCAACAAGAAGTCTATGTTTATGTTT';
    say $tmpfa1 '>HWI-ST765:123:D0TEDACXX:5:1101:17034:2404/1';
    say $tmpfa1 'CATTTGCGGTACTTCACACTAGCATGATGATGAAGGGTGCAATCGTTGCACAAGGGTGCAGTTCCGTTGTATGGCTTTCTAGCAGGGGGTTGAGTTGGTTG';

    my $tmpfa2 = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
                                  DIR      => 't',
                                  SUFFIX   => ".fasta",
				  UNLINK   => 0 );

    my $tmpfa2_name = $tmpfa2->filename;

    say $tmpfa2 '>HWI-ST765:123:D0TEDACXX:5:1101:2872:2088/2';
    say $tmpfa2 'ATTGTGTTATAAAGTTTATTTCTATTTCCGTTGCAACTTAAATCTGATTTACATTCATTTTACTTAAACAAACACAATCAAAAGAAACTCAGATCCTACAA';
    say $tmpfa2 '>HWI-ST765:123:D0TEDACXX:5:1101:6511:2225/2';
    say $tmpfa2 'AAAGAAGACGGTGACTGAGTGCAATGATTTGTTCGAGAGTTTTGCACATTCTGATATGGACTACAGCACTGCCAGCAGGACTTCCATTCCTGTTACTACCA';
    say $tmpfa2 '>HWI-ST765:123:D0TEDACXX:5:1101:12346:2094/2';
    say $tmpfa2 'TCAATTAAGTCCAAATAAAGTAATCAATGCAATTGCCAAAGAGTCCGCGGCAACGGCGCCAAAAAACTTGATGTGCTAAAAGTAGTTTAATAAAACAACTA';
    say $tmpfa2 '>HWI-ST765:123:D0TEDACXX:5:1101:16473:2185/2';
    say $tmpfa2 'GTTGATTATGTTCTCATGCATACAGGGGTATGGCGATCCCGGACCCAAGTCAGCGACATGGACTCAAGCTTTTAATCGAAGACTACCCGTACGCTTCTGAC';
    say $tmpfa2 '>HWI-ST765:123:D0TEDACXX:5:1101:16583:2310/2';
    say $tmpfa2 'TCCCTTTTATTTATTTTGTTTTTATGAACTTTTGTGATATTGTTGATCACTAGCAGTGGTGTAGCATTGGTGCTATTTGGTACGGTTTACCCTGCACGCGG';
    say $tmpfa2 '>HWI-ST765:123:D0TEDACXX:5:1101:17034:2404/2';
    say $tmpfa2 'AAAGGTGACAAGAAACCAATCGAAGAATCAAAACCTAAGGATAAACAGACTGAATCCTCCAAGAAGTCAAAGAAGCGGAAGGCTTCTCAGAACTTCACCGT';
    
    return [$tmpfa1_name, $tmpfa2_name];
}
