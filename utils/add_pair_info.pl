#!/usr/bin/env perl

=head1 NAME 
                                                                       
add_pair_info.pl - Append pair information to the sequence name

=head1 SYNOPSIS    
 
add_pair_info.pl -i seq_1.fq -o seq_1_pair.fq -p 1

=head1 DESCRIPTION
     
With Casava 1.8+ and other data, the pair information is in the comment, separated
from the sequence name by a space. This information is usually lost if any quality trimming
or sequence sampling is done. This script will restore the pair information in the Casava 1.4 style.
The input may be FastA or FastQ.

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

The file of sequences without the pair information in the sequence name.

=item -o, --outfile
              
The file of sequences that will contain the sequence names with the pair information.

=item -p, pairnum

The number to append to the sequence name. Integer. Must be 1 or 2.

=back

=head1 OPTIONS

=over 2

=item -c, --compress

The output files should be compressed. If given, this option must be given the arguments 'gzip' to compress with gzip,
or 'bzip2' to compress with bzip2.

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
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
use Pod::Usage;

my $infile;
my $outfile;
my $pairnum;
my $pair;
my $compress;
my $help;
my $man;

GetOptions(
	   'i|infile=s'    => \$infile,
	   'o|outfile=s'   => \$outfile,
	   'p|pairnum=i'   => \$pairnum,
	   'c|compress=s'  => \$compress,
	   'h|help'        => \$help,
	   'm|man'         => \$man,
	   ) || pod2usage( "Try '$0 --man' for more information." );

#
# Check @ARGV
#
usage() and exit(0) if $help;

pod2usage( -verbose => 2 ) if $man;

if (!$infile || !$outfile || !$pairnum) {
    say "\nERROR: Command line not parsed correctly. Check input.\n";
    usage();
    exit(1);
}

if ($compress) {
    unless ($compress =~ /gzip/i || $compress =~ /bzip2/i) {
        say "\nERROR: $compress is not recognized as an argument to the --compress option. Must be 'gzip' or 'bzip2'. Exiting";
        exit(1);
    }
}

if ($pairnum == 1) {
    $pair = "/1";
}
elsif ($pairnum == 2) {
    $pair = "/2";
}
else {
    say "\nERROR: $pairnum is not correct. Must be 1 or 2. Exiting.";
    exit(1);
}

my $fh = get_fh($infile);
open my $out, '>', $outfile or die "\nERROR: Could not open file: $!\n";

my @aux = undef;
my ($name, $comm, $seq, $qual);

while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
    if (defined $qual) {
	say $out join "\n", "@".$name.$pair, $seq, '+', $qual;
    }
    else {
	say $out join "\n", ">".$name.$pair, $seq;
    }
}
close $fh;
close $out;

compress($outfile) if $compress;
exit;
#
# subs
#
sub get_fh {
    my ($file) = @_;

    my $fh;
    if ($file =~ /\.gz$/) {
        open $fh, '-|', 'zcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /\.bz2$/) {
        open $fh, '-|', 'bzcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    else {
        open $fh, '<', $file or die "\nERROR: Could not open file: $file\n";
    }

    return $fh;
}

sub compress {
    my ($outfile) = @_;
    if ($compress =~ /gzip/i) {
        my $outfilec = $outfile.".gz";
        gzip $outfile => $outfilec or die "gzip failed: $GzipError\n";
        unlink $outfile;
    }
    elsif ($compress =~ /bzip2/i) {
	my $outfilec = $outfile.".bz2";
        bzip2 $outfile => $outfilec or die "bzip2 failed: $Bzip2Error\n";
        unlink $outfile;
    }
}

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

sub mk_key { join "\N{INVISIBLE SEPARATOR}", @_ }

sub mk_vec { split "\N{INVISIBLE SEPARATOR}", shift }

sub usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script [-i] [-o] [-p] [-h] [-m]

Required:
    -i|infile        :       The file of sequences without the pair information in the sequence name.
    -o|outfile       :       The file of sequences that will contain the sequence names with the pair information.
    -p|pairnum       :       The number to append to the sequence name. Integer (Must be 1 or 2).

Options:
    -c|compress       :       Compress the output files. Options are 'gzip' or 'bzip2' (Default: No).
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}
