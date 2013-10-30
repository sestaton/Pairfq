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

=item -i, --infile

The file of sequences without the pair information in the sequence name.

=item -o, --outfile
              
The file of sequences that will contain the sequence names with the pair information.

=item -p, pairnum

The number to append to the sequence name. Integer. Must be 1 or 2.

=back

=head1 OPTIONS

=over 2

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;

my $infile;
my $outfile;
my $pairnum;
my $pair;

GetOptions(
	   'i|infile=s'    => \$infile,
	   'o|outfile=s'   => \$outfile,
	   'p|pairnum=i'   => \$pairnum,
	   );

if (!$infile || !$outfile || !$pairnum) {
    say "\nERROR: Command line not parsed correctly. Exiting.";
    usage();
    exit(1);
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

open my $f, '<', $infile or die "\nERROR: Could not open file: $!\n";
open my $out, '>', $outfile or die "\nERROR: Could not open file: $!\n";

my @aux = undef;
my ($name, $comm, $seq, $qual);

while (($name, $comm, $seq, $qual) = readfq(\*$f, \@aux)) {
    if (defined $qual) {
	say $out join "\n", "@".$name.$pair, $seq, '+', $qual;
    }
    else {
	say $out join "\n", ">".$name.$pair, $seq;
    }
}
close $f;
close $out;

exit;
#
# subs
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
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}
