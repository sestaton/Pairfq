#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use Getopt::Long;

my $usage = "$0 -i file -o fixed -p num";
my $infile;
my $outfile;
my $pairnum;
my $pair;

GetOptions(
	   'i|infile=s'    => \$infile,
	   'o|outfile=s'   => \$outfile,
	   'p|pairnum=i'   => \$pairnum,
	   );

say $usage and exit(1) if !$infile or !$pairnum or !$outfile;

if ($pairnum == 1) {
    $pair = "/1";
}
elsif ($pairnum == 2) {
    $pair = "/2";
}
else {
    say "ERROR: $pairnum is not correct. Must be 1 or 2. Exiting.";
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
