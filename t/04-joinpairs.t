#!perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use File::Temp;
use File::Spec;
use File::Basename;
use autodie qw(open);
use Test::More tests => 30;

#TODO: Add tests that sequences and IDs are correct between tests 
my $cmd     = File::Spec->catfile('blib', 'bin', 'pairfq');
my $fq_data = _build_fq_data();
my $fa_data = _build_fa_data();

joinpairs_inmemory($cmd, $fq_data, $fa_data);
joinpairs_stdout($cmd, $fq_data, $fa_data);
joinpairs_ondisk($cmd, $fq_data, $fa_data);

#
# methods
#
sub joinpairs_inmemory {
    my ($cmd, $fq_data, $fa_data) = @_;
    my $tmpfq_out = File::Temp->new( TEMPLATE => "pairfq_fq_XXXX",
				     DIR      => 't',
				     SUFFIX   => ".fastq",
				     UNLINK   => 0 );
    
    my $tmpfa_out = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
				     DIR      => 't',
				     SUFFIX   => ".fasta",
				     UNLINK   => 0 );
    
    system("$cmd joinpairs -f $fq_data->[0] -r $fq_data->[1] -o $tmpfq_out") == 0 
	or die "system failed: $?";
    system("$cmd joinpairs -f $fa_data->[0] -r $fa_data->[1] -o $tmpfa_out") == 0 
	or die "system failed: $?";

    open my $fqo, '<', $tmpfq_out;
    open my $fao, '<', $tmpfa_out;
    
    my ($fq_fct, $fq_rct, $fa_fct, $fa_rct) = (0, 0, 0, 0);
    my (%fqpairs, %fapairs);
    while (<$fqo>) {
	if (/^\@(\S+)\/1$/) {
	    $fqpairs{$1}++;
	    $fq_fct++;
	}
	elsif (/^\@(\S+)\/2$/) {
	    $fqpairs{$1}++;
	    $fq_rct++;
	}
    }
    close $fqo;
    
    my $fqpaired = (keys %fqpairs);
    my @fqpairnum = (values %fqpairs);
    my $allfqpaired = all_the_same(\@fqpairnum);
    ok($allfqpaired == 1, 'All fastq reads were paired');

    is($fqpaired,          6, 'Correct number of fastq pairs');
    is($fq_fct,            6, 'Correct number of forward fastq reads paired');
    is($fq_rct,            6, 'Correct number of reverse fastq reads paired');
    is($fq_fct + $fq_fct, 12, 'Correct number of total fastq reads paired');
    unlink $tmpfq_out;
    
    while (<$fao>) {
	if (/^\>(\S+)\/1$/) {
	    $fapairs{$1}++;
	    $fa_fct++;
	}
	elsif (/^\>(\S+)\/2$/) {
	    $fapairs{$1}++;
	    $fa_rct++;
	}
    }
    close $fao;
    
    my $fapaired = (keys %fapairs);
    my @fapairnum = (values %fapairs);
    my $allfapaired = all_the_same(\@fapairnum);
    ok($allfapaired == 1, 'All fasta reads were paired');

    is($fapaired,          6, 'Correct number of fasta pairs');
    is($fa_fct,            6, 'Correct number of forward fasta reads paired');
    is($fa_rct,            6, 'Correct number of reverse fasta reads paired');
    is($fa_fct + $fa_fct, 12, 'Correct number of total fasta reads paired');
    unlink $tmpfa_out;
}

sub joinpairs_stdout {
    my ($cmd, $fq_data, $fa_data) = @_;
    my $tmpfq_out = File::Temp->new( TEMPLATE => "pairfq_fq_XXXX",
				     DIR      => 't',
				     SUFFIX   => ".fastq",
				     UNLINK   => 0 );
    
    my $tmpfa_out = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
				     DIR      => 't',
				     SUFFIX   => ".fasta",
				     UNLINK   => 0 );

    {
        my $fqstdo = \*STDOUT;
        open $fqstdo, '>', $tmpfq_out;
        system("$cmd joinpairs -f $fq_data->[0] -r $fq_data->[1] -o -") == 0 
	    or die "system failed: $?";
	#close $fqstdo;
    }
    {
        my $fastdo = \*STDOUT;
        open $fastdo, '>', $tmpfa_out;
 	system("$cmd joinpairs -f $fa_data->[0] -r $fa_data->[1] -o -") == 0 
	    or die "system failed: $?";
	#close $fastdo;
    }

    open my $fqo, '<', $tmpfq_out;
    open my $fao, '<', $tmpfa_out;
    
    my ($fq_fct, $fq_rct, $fa_fct, $fa_rct) = (0, 0, 0, 0);
    my (%fqpairs, %fapairs);
    while (<$fqo>) {
	if (/^\@(\S+)\/1$/) {
	    $fqpairs{$1}++;
	    $fq_fct++;
	}
	elsif (/^\@(\S+)\/2$/) {
	    $fqpairs{$1}++;
	    $fq_rct++;
	}
    }
    close $fqo;
    
    my $fqpaired = (keys %fqpairs);
    my @fqpairnum = (values %fqpairs);
    my $allfqpaired = all_the_same(\@fqpairnum);
    ok($allfqpaired == 1, 'All fastq reads were paired');

    is($fqpaired,          6, 'Correct number of fastq pairs');
    is($fq_fct,            6, 'Correct number of forward fastq reads paired');
    is($fq_rct,            6, 'Correct number of reverse fastq reads paired');
    is($fq_fct + $fq_fct, 12, 'Correct number of total fastq reads paired');
    unlink $tmpfq_out;
    
    while (<$fao>) {
	if (/^\>(\S+)\/1$/) {
	    $fapairs{$1}++;
	    $fa_fct++;
	}
	elsif (/^\>(\S+)\/2$/) {
	    $fapairs{$1}++;
	    $fa_rct++;
	}
    }
    close $fao;
    
    my $fapaired = (keys %fapairs);
    my @fapairnum = (values %fapairs);
    my $allfapaired = all_the_same(\@fapairnum);
    ok($allfapaired == 1, 'All fasta reads were paired');

    is($fapaired,          6, 'Correct number of fasta pairs');
    is($fa_fct,            6, 'Correct number of forward fasta reads paired');
    is($fa_rct,            6, 'Correct number of reverse fasta reads paired');
    is($fa_fct + $fa_fct, 12, 'Correct number of total fasta reads paired');
    unlink $tmpfa_out;
}

sub joinpairs_ondisk {
    my ($cmd, $fq_data, $fa_data) = @_;
    my $tmpfq_out = File::Temp->new( TEMPLATE => "pairfq_fq_XXXX",
                                     DIR      => 't',
                                     SUFFIX   => ".fastq",
                                     UNLINK   => 0 );
    
    my $tmpfa_out = File::Temp->new( TEMPLATE => "pairfq_fa_XXXX",
                                     DIR      => 't',
                                     SUFFIX   => ".fasta",
                                     UNLINK   => 0 );
    
    system("$cmd joinpairs -f $fq_data->[0] -r $fq_data->[1] -o $tmpfq_out -idx") == 0 
	or die "system failed: $?";
    system("$cmd joinpairs -f $fa_data->[0] -r $fa_data->[1] -o $tmpfa_out -idx") == 0 
	or die "system failed: $?";

    open my $fqo, '<', $tmpfq_out;
    open my $fao, '<', $tmpfa_out;
    
    my ($fq_fct, $fq_rct, $fa_fct, $fa_rct) = (0, 0, 0, 0);
    my (%fqpairs, %fapairs);
    while (<$fqo>) {
        if (/^\@(\S+)\/1$/) {
            $fqpairs{$1}++;
            $fq_fct++;
        }
        elsif (/^\@(\S+)\/2$/) {
            $fqpairs{$1}++;
            $fq_rct++;
        }
    }
    close $fqo;
    
    my $fqpaired = (keys %fqpairs);
    my @fqpairnum = (values %fqpairs);
    my $allfqpaired = all_the_same(\@fqpairnum);
    ok($allfqpaired == 1, 'All fastq reads were paired');

    is($fqpaired,          6, 'Correct number of fastq pairs');
    is($fq_fct,            6, 'Correct number of forward fastq reads paired');
    is($fq_rct,            6, 'Correct number of reverse fastq reads paired');
    is($fq_fct + $fq_fct, 12, 'Correct number of total fastq reads paired');
    unlink $fq_data->[0], $fq_data->[1], $tmpfq_out;

    while (<$fao>) {
        if (/^\>(\S+)\/1$/) {
            $fapairs{$1}++;
            $fa_fct++;
        }
        elsif (/^\>(\S+)\/2$/) {
            $fapairs{$1}++;
            $fa_rct++;
        }
    }
    close $fao;
    
    my $fapaired = (keys %fapairs);
    my @fapairnum = (values %fapairs);
    my $allfapaired = all_the_same(\@fapairnum);
    ok($allfapaired == 1, 'All fasta reads were paired');

    is($fapaired,          6, 'Correct number of fasta pairs');
    is($fa_fct,            6, 'Correct number of forward fasta reads paired');
    is($fa_rct,            6, 'Correct number of reverse fasta reads paired');
    is($fa_fct + $fa_fct, 12, 'Correct number of total fasta reads paired');
    unlink $fa_data->[0], $fa_data->[1], $tmpfa_out;
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

# This solution is provided by Sinan Unur
# http://stackoverflow.com/a/2305879/1543853
sub all_the_same {
    my ($ref) = @_;
    return 1 unless @$ref;
    my $cmpv = \ $ref->[-1];
    for my $i (0 .. $#$ref - 1)  {
        my $this = \ $ref->[$i];
        return unless defined $$cmpv == defined $$this;
        return if defined $$this
            and ( $$cmpv ne $$this );
    }
    return 1;
}
