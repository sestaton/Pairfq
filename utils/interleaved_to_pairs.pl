#!/usr/bin/env perl

=head1 NAME 
                                                                       
pairfq.pl - Match paired-end sequences from separate FastA/Q files

=head1 SYNOPSIS    
 
pairfq.pl -f s_1_1_trim.fq -r s_1_2_trim.fq -fp s_1_1_trim_paired.fq -rp s_1_2_trim_paired.fq -fs s_1_1_trim_unpaired.fq -rs s_1_2_trim_unpaired.fq

=head1 DESCRIPTION
     
Re-pair paired-end sequences that may have been separated by quality trimming.
This script also writes the unpaired forward and reverse sequences to separate 
files so that they may be used for assembly or mapping. The input may be FastA
or FastQ format in either Illumina 1.3+ or Illumina 1.8 format.

=head1 DEPENDENCIES

Only core Perl is required, no external dependencies. See below for information
on which Perls have been tested.

=head1 LICENSE
 
The MIT License should included with the project. If not, it can be found at: http://opensource.org/licenses/mit-license.php

Copyright (C) 2013 S. Evan Staton
 
=head1 TESTED WITH:

=over

=item *
Perl 5.14.1 (Red Hat Enterprise Linux Server release 5.7 (Tikanga))

=item *
Perl 5.14.2 (Red Hat Enterprise Linux Desktop release 6.2 (Santiago); Fedora 17)

=item *
Perl 5.18.0 (Red Hat Enterprise Linux Server release 5.9 (Tikanga))

=back

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -f, --forward

The file of forward sequences from an Illumina paired-end sequencing run.

=item -r, --reverse                                                                                                                                                       
The file of reverse sequences from an Illumina paired-end sequencing run.

=item -fp, --forw_paired

The output file to place the paired forward reads.

=item -rp, --rev_paired                                                                                                                                                  
The output file to place the paired reverse reads. 

=item -fs, --forw_unpaired                                                                                                                                                  
The output file to place the unpaired forward reads. 

=item -rs, --rev_unpaired                                                                                                                                                  
The output file to place the unpaired reverse reads. 

=back

=head1 OPTIONS

=over 2

=item -im, --in_memory

Construct the database in memory. May be faster, but will obviously use more memory.

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut

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
	say $f join "\n", ">".$seqid, $seq if $seqid =~ m/1$|\s+1/;
	say $r join "\n", ">".$seqid, $seq if $seqid =~ m/2$|\s+2/;
    }
}

close $in;
close $f;
close $r;

sub usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script [-f] [-r] [-fp] [-rp] [-fs] [-rs] [-im] [-h] [-im]

Required:
    -f|forward        :       File of foward reads (usually with "/1" or " 1" in the header).
    -r|reverse        :       File of reverse reads (usually with "/2" or " 2" in the header).
    -fp|forw_paired   :       Name for the file of paired forward reads.
    -rp|rev_paired    :       Name for the file of paired reverse reads.
    -fs|forw_unpaired :       Name for the file of singleton forward reads.
    -rs|rev_unpaired  :       Name for the file of singleton reverse reads.

Options:
    -im|in_memory     :       Construct a database in memory for faster execution.
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}
