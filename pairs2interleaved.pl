#!/usr/bin/env perl

use 5.010;
use utf8;
use strict;
use warnings;
use warnings FATAL => "utf8";
use charnames qw(:full :short);
use autodie qw(open);
use DB_File;
use Getopt::Long;
use Data::Dump qw(dd);
use IO::Handle;

my $usage = "$0 -f forward -r reverse -o interleaved.fas [-im]";
my $forward;
my $reverse;
my $outfile;
my $memory;

GetOptions(
	   'f|forward=s'   => \$forward,
	   'r|reverse=s'   => \$reverse,
	   'o|outfile=s'   => \$outfile,
	   'im|memory'     => \$memory,
	   );

say $usage and exit(1) if !$forward or !$reverse or !$outfile;

my $pairs = store_pair($forward);
open my $r, '<', $reverse or die "\nERROR: Could not open file: $!\n";
open my $out, '>', $outfile or die "\nERROR: Could not open file: $!\n";
binmode $out, ":utf8";

my (@faux, @raux) = (undef, undef);
my ($fname, $fcomm, $fseq, $fqual, $rname, $rcomm, $rseq, $rqual);
my ($fct, $rct) = (0, 0);

#dd $pairs and exit;

while (($rname, $rcomm, $rseq, $rqual) = readfq(\*$r, \@raux)) {
    $rname =~ s/\d$// if defined $rname && $rname =~ /\/\d$/;
    if (defined $rcomm && $rcomm =~ /^\d/) {
	$rcomm =~ s/^\d//;
	$rname = mk_key($rname, $rcomm);
    }
    if (defined $rname && exists $pairs->{$rname}) {
	if (defined $rcomm) {
	    my ($name, $comm) = mk_vec($rname);
	    if (defined $rqual) {
		my ($seqf, $qualf) = mk_vec($pairs->{$rname});
		say $out join "\n", "@".$rname.q{ 1}.$comm, $seqf, '+', $qualf;
		say $out join "\n", "@".$rname.q{ 2}.$comm, $rseq, '+', $rqual;
	    }
	    else {
		say $out join "\n", ">".$name.q{ 1}.$comm, $pairs->{$rname};
		say $out join "\n", ">".$name.q{ 2}.$comm, $rseq;
	    }
	}
	else {
	    if (defined $rqual) {
		my ($seqf, $qualf) = mk_vec($pairs->{$rname});
		say $out join "\n", "@".$rname.q{/1}, $seqf, '+', $qualf;
		say $out join "\n", "@".$rname.q{/2}, $rseq, '+', $rqual;
	    }
	    else {
		say $out join "\n", ">".$rname.q{/1}, $pairs->{$rname};
		say $out join "\n", ">".$rname.q{/2}, $rseq;                                               
	    }
	}
    }
}
close $r;
close $out;

exit;
#
# subs
#
sub store_pair {
    my ($file) = @_;

    my %pairs;
    $DB_BTREE->{cachesize} = 100000;
    $DB_BTREE->{flags} = R_DUP;
    my $db_file = "pairs_to_interleaved.bdb";
    unlink $db_file if -e $db_file;

    unless (defined $memory) {
	tie %pairs, 'DB_File', $db_file, O_RDWR|O_CREAT, 0666, $DB_BTREE
	    or die "\nERROR: Could not open DBM file $db_file: $!\n";
    }

    my @faux = undef;
    my ($fname, $fcomm, $fseq, $fqual);
    my ($fct, $rct) = (0, 0);
    open my $f, '<', $file or die "\nERROR: Could not open file: $!\n";

    while (($fname, $fcomm, $fseq, $fqual) = readfq(\*$f, \@faux)) {
	$fname =~ s/\d$// if $fname =~ /\/\d$/;
	if (defined $fcomm && $fcomm =~ /^\d/) {
	    $fcomm =~ s/^\d//;
	    $fname = mk_key($fname, $fcomm);
	}
	$pairs{$fname} = mk_key($fseq, $fqual) if defined $fqual;
	$pairs{$fname} = $fseq if !defined $fqual;
    }
    close $f;
    
    return \%pairs;
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
