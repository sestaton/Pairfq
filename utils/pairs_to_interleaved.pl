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
