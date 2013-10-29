#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Getopt::Long;
use autodie qw(open);

my $usage = "\n$0 -i interleaved.fasta -f reads_1.fas -r reads_2.fas\n\n";
my $forward;
my $reverse;
my $infile; 

GetOptions(
	   'i|infile=s'  => \$infile,
	   'f|forward=s' => \$forward,
	   'r|reverse=s' => \$reverse,
	   );

die $usage if !$infile or !$forward or !$reverse;

open my $in, '<', $infile;
open my $f, '>', $forward;
open my $r, '>', $reverse; 

{
    local $/ = '>';

    while (my $line = <$in>) {
        chomp $line;
        my ($seqid, @seqparts) = split /\n/, $line;
        my $seq = join '', @seqparts;
        next unless defined $seqid && defined $seq;
	say $f join "\n", ">".$seqid, $seq if $seqid =~ m/1$/;
	say $r join "\n", ">".$seqid, $seq if $seqid =~ m/2$/;
    }
}

close $i;
close $f;
close $r;
