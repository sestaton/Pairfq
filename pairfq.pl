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
use Encode qw(encode decode);
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use DB_File;
use Pod::Usage;

my ($fread, $rread, $fpread, $rpread, $fsread, $rsread, $memory, $help, $man);

GetOptions(
           'f|forward=s'        => \$fread,
           'r|reverse=s'        => \$rread,
           'fp|forw_paired=s'   => \$fpread,
           'rp|rev_paired=s'    => \$rpread,
           'fs|forw_unpaired=s' => \$fsread,
           'rs|rev_unpaired=s'  => \$rsread,
           'im|in_memory'       => \$memory,
           'h|help'             => \$help,
           'm|man'              => \$man,
          ) || pod2usage( "Try '$0 --man' for more information." );

#
# Check @ARGV
#
usage() and exit(0) if $help;

pod2usage( -verbose => 2 ) if $man;

if (!$fread  || !$rread  ||
    !$fpread || !$rpread ||
    !$fsread || !$rsread) {
    say "\nERROR: Command line not parsed correctly. Check input.\n";
    usage();
    exit(1);
}

my ($rseqpairs, $db_file, $rct) = store_pair($rread);

my @faux = undef;
my ($fname, $fcomm, $fseq, $fqual, $forw_id, $rev_id, $fname_enc);
my ($fct, $fpct, $rpct, $pct, $fsct, $rsct, $sct) = (0, 0, 0, 0, 0, 0, 0);

open my $f, '<', $fread or die "\nERROR: Could not open file: $fread\n";
open my $fp, '>', $fpread or die "\nERROR: Could not open file: $fpread\n";
open my $rp, '>', $rpread or die "\nERROR: Could not open file: $rpread\n";
open my $fs, '>', $fsread or die "\nERROR: Could not open file: $fsread\n";

binmode $fp, ":utf8";
binmode $rp, ":utf8";
binmode $fs, ":utf8";

while (($fname, $fcomm, $fseq, $fqual) = readfq(\*$f, \@faux)) {
    $fct++;
    if ($fname =~ /(\/\d)$/) {
	$fname =~ s/\/\d//;
    }
    elsif (defined $fcomm && $fcomm =~ /^\d/) {
        $fcomm =~ s/^\d//;
	$fname = mk_key($fname, $fcomm);
    }
    else {
	say "\nERROR: Could not determine FastA/Q format. ".
	    "Please see https://github.com/sestaton/Pairfq or the README for supported formats. Exiting.\n";
	exit(1);
    }
    
    if ($fname =~ /\N{INVISIBLE SEPARATOR}/) {
	my ($name, $comm) = mk_vec($fname);
	$forw_id = $name.q{ 1}.$comm;
	$rev_id  = $name.q{ 2}.$comm;
    }
    $fname_enc = encode('UTF-8', $fname);
    if (exists $rseqpairs->{$fname_enc}) {
	$fpct++; $rpct++;
	if (defined $fqual) {
	    my ($rread, $rqual) = mk_vec($rseqpairs->{$fname});
	    if ($fname =~ /\N{INVISIBLE SEPARATOR}/) {
		say $fp join "\n","@".$forw_id, $fseq, "+", $fqual;
		say $rp join "\n","@".$rev_id, $rread, "+", $rqual;
	    } 
	    else {
		say $fp join "\n","@".$fname.q{/1}, $fseq, "+", $fqual;
                say $rp join "\n","@".$fname.q{/2}, $rread, "+", $rqual;
	    }
	} 
	else {
	    if ($fname =~ /\N{INVISIBLE SEPARATOR}/) {
		say $fp join "\n",">".$forw_id, $fseq;
		say $rp join "\n",">".$rev_id, $rseqpairs->{$fname_enc};
	    } 
	    else {
                say $fp join "\n",">".$fname.q{/1}, $fseq;
                say $rp join "\n",">".$fname.q{/2}, $rseqpairs->{$fname_enc};
            }
	}
        delete $rseqpairs->{$fname_enc};
    } 
    else {
	$fsct++;
	if (defined $fqual) {
	    if ($fname =~ /\N{INVISIBLE SEPARATOR}/) {
		say $fs join "\n","@".$forw_id, $fseq, "+", $fqual;
	    } 
	    else {
		say $fs join "\n","@".$fname.q{/1}, $fseq, "+", $fqual;
	    }
	} 
	else {
	    if ($fname =~ /\N{INVISIBLE SEPARATOR}/) {
		say $fs join "\n",">".$forw_id, $fseq;
	    } 
	    else {
		say $fs join "\n",">".$fname.q{/1}, $fseq;
	    }
	}
    }
}
close $f;
close $fp;
close $rp;
close $fs;

open my $rs, '>', $rsread or die "\nERROR: Could not open file: $rsread\n";
binmode $rs, ":utf8";

my $rev_id_up;
while (my ($rname_up, $rseq_up) = each %$rseqpairs) {
    $rsct++;
    if ($rname_up =~ /\N{INVISIBLE SEPARATOR}/) {
	my ($uname, $uid) = mk_vec($rname_up);
	$rev_id_up .= $uname.q{ 2}.$uid;
    }
    if ($rseq_up =~ /\N{INVISIBLE SEPARATOR}/) {
	my ($rread_up, $rqual_up) = mk_vec($rseq_up);
	if ($rname_up =~ /\|/) {
	    say $rs join "\n","@".$rev_id_up, $rread_up, "+", $rqual_up;
	} 
	else {
	    say $rs join "\n","@".$rname_up."/2", $rread_up, "+", $rqual_up;
	}
    } 
    else {
	if ($rname_up =~ /\N{INVISIBLE SEPARATOR}/) {
	    say $rs join "\n",">".$rev_id_up, $rseq_up;
	} 
	else {
	    say $rs join "\n",">".$rname_up."/2", $rseq_up;
	}
    }
    $rev_id_up = undef;
}
close $rs;

$pct = $fpct + $rpct;
$sct = $fsct + $rsct;
untie %$rseqpairs if defined $memory;
unlink $db_file if -e $db_file;

say "Total forward reads in $fread:              $fct";
say "Total reverse reads in $rread:              $rct";
say "Total paired reads in $fpread and $rpread:  $pct";
say "Total forward paired reads in $fpread:      $fpct";
say "Total reverse paired reads in $rpread:      $rpct";
say "Total upaired reads in $fsread and $rsread: $sct"; 
say "Total forward unpaired reads in $fsread:    $fsct";
say "Total reverse unpaired reads in $rsread:    $rsct";

exit;
#
# Subs
#
sub store_pair {
    my ($file) = @_;

    my $rct = 0;
    my %rseqpairs;
    $DB_BTREE->{cachesize} = 100000;
    $DB_BTREE->{flags} = R_DUP;
    my $db_file = "pairfq.bdb";
    unlink $db_file if -e $db_file;

    unless (defined $memory) { 
	tie %rseqpairs, 'DB_File', $db_file, O_RDWR|O_CREAT, 0666, $DB_BTREE
	    or die "\nERROR: Could not open DBM file $db_file: $!\n";
    }

    my @raux = undef;
    my ($rname, $rcomm, $rseq, $rqual, $rname_k, $rname_enc);

    open my $r, '<', $file or die "\nERROR: Could not open file: $file\n";

    {
	local @SIG{qw(INT TERM HUP)} = sub { if (defined $memory && -e $db_file) { untie %rseqpairs; unlink $db_file if -e $db_file; } };

	while (($rname, $rcomm, $rseq, $rqual) = readfq(\*$r, \@raux)) {
	    $rct++;
	    if ($rname =~ /(\/\d)$/) {
		$rname =~ s/$1//;
	    }
	    elsif (defined $rcomm && $rcomm =~ /^\d/) {
		$rcomm =~ s/^\d//;
		$rname_k = mk_key($rname, $rcomm);
		$rname_enc = encode('UTF-8', $rname_k);
	    }
	    else {
		say "\nERROR: Could not determine FastA/Q format. ".
		    "Please see https://github.com/sestaton/Pairfq or the README for supported formats. Exiting.\n";
		exit(1);
	    }
	    $rseqpairs{$rname_enc} = mk_key($rseq, $rqual) if defined $rqual && defined $rcomm;
	    $rseqpairs{$rname} = mk_key($rseq, $rqual) if defined $rqual;
	    $rseqpairs{$rname_enc} = $rseq if !defined $rqual && defined $rcomm;
	    $rseqpairs{$rname} = $rseq if !defined $rqual && !defined $rcomm;
	}
	close $r;
    }
    return(\%rseqpairs, $db_file, $rct);
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
