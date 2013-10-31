#!/usr/bin/env perl

=head1 NAME 
                                                                       
pairs_to_interleaved.pl - Create an interleaved FastA/Q file for assembly or mapping

=head1 SYNOPSIS    
 
pairs_to_interleaved.pl -f seq_1_p.fq -r seq_2_p.fq -o seqs_interl.fq

=head1 DESCRIPTION
     
For some assembly programs, such as Velvet, paired-end sequence files must be interleaved.
This program creates a single interleaved file from two paired files, which have been created
by pairfq.pl (or some other method). The input may be FastA or FastQ.

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

The file of paired forward sequences from an Illumina paired-end sequencing run.

=item -r, --reverse                                                                                                                                                       
The file of paired reverse sequences from an Illumina paired-end sequencing run.

=item -o, --outfile

The interleaved file to produce from the forward and reverse files.

=back

=head1 OPTIONS

=over 2

=item -im, --memory

The computation should be done in memory instead of on the disk. This will be faster, but may use a large amount
of RAM if there are many millions of sequences in each input file.

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=cut

use 5.012;
use utf8;
use strict;
use warnings;
use warnings FATAL => "utf8";
use charnames qw(:full :short);
use Encode qw(encode);
use File::Basename;
use DB_File;
use DBM_Filter;
use Getopt::Long;
use Pod::Usage;

my $forward;
my $reverse;
my $outfile;
my $memory;
my $help;
my $man;

GetOptions(
	   'f|forward=s'   => \$forward,
	   'r|reverse=s'   => \$reverse,
	   'o|outfile=s'   => \$outfile,
	   'im|memory'     => \$memory,
	   'h|help'        => \$help,
	   'm|man'         => \$man,
	   ) || pod2usage( "Try '$0 --man' for more information." );;

#
# Check @ARGV
#
usage() and exit(0) if $help;

pod2usage( -verbose => 2 ) if $man;

if (!$forward || !$reverse || !$outfile) {
    say "\nERROR: Command line not parsed correctly. Exiting.";
    usage();
    exit(1);
}

my ($pairs, $db_file) = store_pair($forward);
open my $r, '<', $reverse or die "\nERROR: Could not open file: $!\n";
open my $out, '>', $outfile or die "\nERROR: Could not open file: $!\n";
binmode $out, ":utf8";

my @raux = undef;
my ($rname, $rcomm, $rseq, $rqual, $forw_id, $rev_id);

while (($rname, $rcomm, $rseq, $rqual) = readfq(\*$r, \@raux)) {
    $rname =~ s/\d$// if defined $rname && $rname =~ /\/\d$/;
    if (defined $rcomm && $rcomm =~ /^\d/) {
	$rcomm =~ s/^\d//;
	$rname = mk_key($rname, $rcomm);
    }

    if ($fname =~ /\N{INVISIBLE SEPARATOR}/) {
        my ($name, $comm) = mk_vec($fname);
        $forw_id = $name.q{ 1}.$comm;
        $rev_id  = $name.q{ 2}.$comm;
    }

    $rname = encode('UTF-8', $rname);
    if (exists $pairs->{$rname}) {
	if (defined $rqual) {
	    my ($seqf, $qualf) = mk_vec($pairs->{$rname});
	    if ($rname =~ /\N{INVISIBLE SEPARATOR}/) {
		say $out join "\n", "@".$forw_id, $seqf, "+", $qualf;
		say $out join "\n", "@".$rev_id, $rseq, "+", $rqual;
	    }
	    else {
		# problem with initialized vars at 141 and 142
		say $out join "\n", "@".$rname.q{/1}, $seqf, "+", $qualf;
                say $out join "\n", "@".$rname.q{/2}, $rseq, "+", $rqual;
	    }
	}
	else {
	    if ($rname =~ /\N{INVISIBLE SEPARATOR}/) {
		say $out join "\n", ">".$forw_id, $seqf;
		say $out join "\n", ">".$rev_id, $rseq;
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

untie %$pairs if defined $memory;
unlink $db_file if -e $db_file;

exit;
#
# subroutines
#
sub store_pair {
    my ($file) = @_;

    my %pairs;
    $DB_BTREE->{cachesize} = 100000;
    $DB_BTREE->{flags} = R_DUP;
    my $db_file = "pairs_to_interleaved.bdb";
    unlink $db_file if -e $db_file;

    unless (defined $memory) {
	my $db = tie %pairs, 'DB_File', $db_file, O_RDWR|O_CREAT, 0666, $DB_BTREE
	    or die "\nERROR: Could not open DBM file $db_file: $!\n";
	$db->Filter_Value_Push("utf8");
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
	$fname = encode('UTF-8', $fname);
	$pairs{$fname} = mk_key($fseq, $fqual) if defined $fqual;
	$pairs{$fname} = $fseq if !defined $fqual;
    }
    close $f;
    
    return (\%pairs, $db_file);
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
USAGE: $script [-f] [-r] [-fp] [-o] [-im] [-h] [-m]

Required:
    -f|forward        :       File of foward reads (usually with "/1" or " 1" in the header).
    -r|reverse        :       File of reverse reads (usually with "/2" or " 2" in the header).
    -o|outfile        :       File of interleaved reads.

Options:
    -im|in_memory     :       Construct a database in memory for faster execution.
                              NB: This may result in large RAM usage for a large number of sequences. 
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}
