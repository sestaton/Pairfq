#!/usr/bin/env perl

=head1 NAME 
                                                                       
interleaved_to_pairs.pl - Split forward and reverse reads into separate files

=head1 SYNOPSIS    
 
interleaved_to_pairs.pl -i seq_interl.fq -f seq_1.fq -r seq_2.fq

=head1 DESCRIPTION
     
Interleaving paired-end files is necessary for some assembly or mapping programs, but
that format is not advantageous for all operations. This script will take the interleaved
file and split the forward and reverse reads in to separate files. The input may be FastA or
FastQ, and the output will be the same as the input format.

=head1 DEPENDENCIES

Only core Perl is required, no external dependencies. See below for information
on which Perls have been tested.

=head1 LICENSE
 
The MIT License should included with the project. If not, it can be found at: http://opensource.org/licenses/mit-license.php

Copyright (C) 2013 S. Evan Staton
 
=head1 TESTED WITH:

=over

=item *
Perl 5.18.0 (Red Hat Enterprise Linux Server release 5.9 (Tikanga))

=back

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

The file of interleaved FastA/Q forward and reverse reads.

=item -f, --forward

The file to place the forward reads.

=item -r, --reverse

The file to place the reverse reads.

=back

=head1 OPTIONS

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut

use 5.012;
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $forward;
my $reverse;
my $infile;
my $help;
my $man;

GetOptions(
	   'i|infile=s'  => \$infile,
	   'f|forward=s' => \$forward,
	   'r|reverse=s' => \$reverse,
	   'h|help'      => \$help,
	   'm|man'       => \$man,
	   ) || pod2usage( "Try '$0 --man' for more information." );

#
# Check @ARGV
#
usage() and exit(0) if $help;

pod2usage( -verbose => 2 ) if $man;

if (!$infile || !$forward || !$reverse) {
    say "\nERROR: Command line not parsed correctly. Check input.\n";
    usage();
    exit(1);
}

open my $in, '<', $infile or die "\nERROR: Could not open file: $!\n";
open my $f, '>', $forward or die "\nERROR: Could not open file: $!\n";
open my $r, '>', $reverse or die "\nERROR: Could not open file: $!\n"; 

my @aux = undef;
my ($name, $comm, $seq, $qual);

while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
    if (defined $comm && $comm =~ /^1/ || $name =~ /1$/) {
	say $f join "\n", ">".$name, $seq if !defined $qual && !defined $comm;
	say $f join "\n", ">".$name.q{ }.$comm, $seq, if !defined $qual && defined $comm;
	say $f join "\n", "@".$name, $seq, '+', $qual if defined $qual && !defined $comm;
	say $f join "\n", "@".$name.q{ }.$comm, $seq, '+', $qual if defined $qual && defined $comm;
    }
    elsif (defined $comm && $comm =~ /^2/ || $name =~ /2$/) {
	say $r join "\n", ">".$name, $seq if !defined $qual && !defined $comm;
	say $r join "\n", ">".$name.q{ }.$comm, $seq if !defined $qual && defined $comm;
	say $r join "\n", "@".$name, $seq, "+", $qual if defined $qual && !defined $comm;
	say $r join "\n", "@".$name.q{ }.$comm, $seq, "+", $qual if defined $qual && defined $comm;
    }
}

close $in;
close $f;
close $r;

exit;
#
# subroutines
#
sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!@$aux);
    return if ($aux->[1]);
    if (!defined($aux->[0])) {
        while (<$fh>) {
            chomp;
            if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                $aux->[0] = $_;
                last;
            }
        }
        if (!defined($aux->[0])) {
            $aux->[1] = 1;
            return;
        }
    }
    my ($name, $comm);
    defined $_ && do {
        ($name, $comm) = /^.(\S+)(?:\s+)(\S+)/ ? ($1, $2) : 
                         /^.(\S+)/ ? ($1, '') : ('', '');
    };
    my $seq = '';
    my $c;
    $aux->[0] = undef;
    while (<$fh>) {
        chomp;
        $c = substr($_, 0, 1);
        last if ($c eq '>' || $c eq '@' || $c eq '+');
        $seq .= $_;
    }
    $aux->[0] = $_;
    $aux->[1] = 1 if (!defined($aux->[0]));
    return ($name, $comm, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $comm, $seq, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq);
}

sub usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script [-i] [-f] [-r] [-h] [-m]

Required:
    -i|infile         :       File of interleaved forward and reverse reads.
    -f|forward        :       File to place the foward reads.
    -r|reverse        :       File to place the reverse reads.

Options:
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}
