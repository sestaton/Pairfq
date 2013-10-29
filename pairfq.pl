#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
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

my %rseqhash;
$DB_BTREE->{cachesize} = 100000;
$DB_BTREE->{flags} = R_DUP;
my $db_file = "pairfq.bdb";
unlink $db_file if -e $db_file;

unless (defined $memory) { 
   tie %rseqhash, 'DB_File', $db_file, O_RDWR|O_CREAT, 0666, $DB_BTREE
       or die "\nERROR: Could not open DBM file $db_file: $!\n";
}

my @raux = undef;
my ($fname, $fcomm, $fseq, $fqual, $fid, $rname, $rcomm, $rseq, $rqual, $rid);
my ($fct, $rct, $fpct, $rpct, $pct, $fsct, $rsct, $sct) = (0, 0, 0, 0, 0, 0, 0, 0);
open my $r, '<', $rread or die "\nERROR: Could not open file: $rread\n";

while (($rname, $rcomm, $rseq, $rqual) = readfq(\*$r, \@raux)) {
    $rct++;
    if ($rname =~ /(\/\d)$/) {
	$rname =~ s/$1//;
    }
    elsif (defined $rcomm && $rcomm =~ /^\d/) {
	$rcomm =~ s/^\d//;
	$rname = mk_key($rname, $rcomm);
    }
    else {
	die "\nERROR: Could not determine FastA/Q format. Please see ... Exiting.\n";
    }

    $rseqhash{$rname} = mk_key($rseq, $rqual) if defined $rqual;
    $rseqhash{$rname} = $rseq if !defined $rqual;
}
close $r;

open my $f, '<', $fread or die "\nERROR: Could not open file: $fread\n";
open my $fp, '>', $fpread or die "\nERROR: Could not open file: $fpread\n";
open my $rp, '>', $rpread or die "\nERROR: Could not open file: $rpread\n";
open my $fs, '>', $fsread or die "\nERROR: Could not open file: $fsread\n";

binmode $fp, ":utf8";
binmode $rp, ":utf8";
binmode $fs, ":utf8";

my ($forw_id, $rev_id);
my @faux = undef;
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
        die "\nERROR: Could not determine FastA/Q format. Please see ... Exiting.\n";
    }

    if ($fname =~ /\N{INVISIBLE SEPARATOR}/) {
	my ($name, $comm) = mk_vec($fname);
	$forw_id = $name.q{ 1}.$comm;
	$rev_id  = $name.q{ 2}.$comm;
    }
    if (exists $rseqhash{$fname}) {
	$fpct++; $rpct++;
	if (defined $fqual) {
	    my ($rread, $rqual) = mk_vec($rseqhash{$fname});
	    if ($fname =~ /\N{INVISIBLE SEPARATOR}/) {
		say $fp join "\n","@".$forw_id, $fseq, "+", $fqual;
		say $rp join "\n","@".$rev_id, $rread, "+", $rqual;
	    } 
	    else {
		say $fp join "\n","@".$fname."/1", $fseq, "+", $fqual;
                say $rp join "\n","@".$fname."/2", $rread, "+", $rqual;
	    }
	} 
	else {
	    if ($fname =~ /\N{INVISIBLE SEPARATOR}/) {
		say $fp join "\n",">".$forw_id, $fseq;
		say $rp join "\n",">".$rev_id, $rseqhash{$fname};
	    } 
	    else {
                say $fp join "\n",">".$fname."/1", $fseq;
                say $rp join "\n",">".$fname."/2", $rseqhash{$fname};
            }
	}
        delete $rseqhash{$fname};
    } 
    else {
	$fsct++;
	if (defined $fqual) {
	    if ($fname =~ /\N{INVISIBLE SEPARATOR}/) {
		say $fs join "\n","@".$forw_id, $fseq, "+", $fqual;
	    } 
	    else {
		say $fs join "\n","@".$fname."/1", $fseq, "+", $fqual;
	    }
	} 
	else {
	    if ($fname =~ /\N{INVISIBLE SEPARATOR}/) {
		say $fs join "\n",">".$forw_id, $fseq;
	    } 
	    else {
		say $fs join "\n",">".$fname."/1", $fseq;
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
while (my ($rname_up, $rseq_up) = each %rseqhash) {
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
untie %rseqhash if defined $memory;
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
    my ($name, $comm) = /^.(\S+)(?:\s+)(\S+)/ ? ($1, $2) : 
                        /^.(\S+)/ ? ($1, '') : ('', '');
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
