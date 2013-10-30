#!/usr/bin/env perl

=head1 NAME 
                                                                       
interleaved_to_pairs.pl - Split forward and reverse reads into separate files

=head1 SYNOPSIS    
 
interleaved_to_pairs.pl -i seq_interl.fq -f seq_1.fq -r seq_2.fq

=head1 DESCRIPTION
     
Interleaving paired-end files is necessary for some assembly or mapping programs, but
that format is not advantageous for all operations. This script will take the interleaved
file and split the forward and reverse reads in to separate files.

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

use 5.010;
use strict;
use warnings;
use Getopt::Long;
use autodie qw(open);

my $forward;
my $reverse;
my $infile; 

GetOptions(
	   'i|infile=s'  => \$infile,
	   'f|forward=s' => \$forward,
	   'r|reverse=s' => \$reverse,
	   );

if (!$infile || !$forward || !$reverse) {
    say "\nERROR: Command line not parsed correctly. Exiting.";
    usage();
    exit(1);
}

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
