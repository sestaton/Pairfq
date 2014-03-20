#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use File::Basename;
use File::Temp;
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
use List::Util qw(max);
use Getopt::Long;
use Pod::Usage;

our $VERSION = '0.09.2';

my $infile;     # input file for 'addinfo' and 'splitpairs' methods
my $outfile;    # output file for 'addinfo' and 'splitpairs' methods
my $fread;      # file of forward reads for 'splitpairs', 'makepairs' and 'joinpairs' methods
my $rread;      # file of forward reads for 'splitpairs', 'makepairs' and 'joinpairs' methods
my $fpread;     # file of paired forward reads for 'makepairs' method
my $rpread;     # file of paired reverse reads for 'makepairs' method
my $fsread;     # file of unpaired forward reads for 'makepairs' method
my $rsread;     # file of unpaired reverse reads for 'makepairs' method 
my $pairnum;    # for the 'addinfo' method
my $compress;   # for any/all the methods
my $uppercase;  # for 'addinfo' method
my $stats;      # currently, for 'makepairs' option only

my $version; 
my $help;
my $man;
my $script = basename($0);

GetOptions(
	   'i|infile=s'         => \$infile,
	   'o|outfile=s'        => \$outfile,
	   'p|pairnum=i'        => \$pairnum, 
	   'f|forward=s'        => \$fread,
	   'r|reverse=s'        => \$rread,
	   'fp|forw_paired=s'   => \$fpread,
	   'rp|rev_paired=s'    => \$rpread,
	   'fs|forw_unpaired=s' => \$fsread,
	   'rs|rev_unpaired=s'  => \$rsread,
	   'c|compress=s'       => \$compress,
	   'uc|uppercase'       => \$uppercase,
	   's|stats'            => \$stats,
	   'version'            => \$version,
	   'h|help'             => \$help,
	   'm|man'              => \$man,
	  ) or pod2usage( "Try '$0 --man' for more information." );

#
# Check @ARGV
#
usage() and exit(0) if $help;

pod2usage( -verbose => 2 ) if $man;

print $VERSION and exit(0) if $version;

my $method = shift;
if (!defined $method) {
    print "\nERROR: Command line not parsed correctly. Check input.\n\n";
    usage();
    exit(1);
}

if ($compress) {
    unless ($compress eq 'gzip' || $compress eq 'bzip2') {
	print "\nERROR: '$compress' is not recognized as an argument to the --compress option.".
	    " Must be 'gzip' or 'bzip2'. Exiting.\n\n";
	exit(1);
    }
}

if ($method eq 'addinfo') {
    if (!$pairnum || !$infile || !$outfile) {
	print "\nERROR: Command line not parsed correctly. Check input.\n\n";
	addinfo_usage();
	exit(1);
    }
    add_pair_info($pairnum, $infile, $outfile, $compress, $uppercase);
}
elsif ($method eq 'makepairs') {
    if (!$fread || !$rread || !$fpread || !$rpread || !$fsread || !$rsread) {
	print "\nERROR: Command line not parsed correctly. Check input.\n\n";
	makepairs_usage();
	exit(1);
    }
    make_pairs_and_singles($fread, $rread, $fpread, $rpread, $fsread, $rsread, $compress, $stats);
}
elsif ($method eq 'joinpairs') {
    if (!$fread || !$rread || !$outfile) {
	print "\nERROR: Command line not parsed correctly. Check input.\n\n";
	joinpairs_usage();
	exit(1);
    }
    pairs_to_interleaved($fread, $rread, $outfile, $compress);
}
elsif ($method eq 'splitpairs') {
    if (!$infile || !$fread || !$rread) {
	print "\nERROR: Command line not parsed correctly. Check input.\n\n";
	splitpairs_usage();
	exit(1);
    }
    interleaved_to_pairs($infile, $fread, $rread, $compress);
}
else {
    print "\nERROR: '$method' is not recognized. See the manual by typing 'perl $script -m',".
	" or see https://github.com/sestaton/Pairfq.\n\n";
    exit(1);
}

exit;
#
# Methods
#
sub add_pair_info {
    my ($pairnum, $infile, $outfile, $compress, $uppercase) = @_;

    my $pair;
    if ($pairnum == 1) {
	$pair = "/1";
    }
    elsif ($pairnum == 2) {
	$pair = "/2";
    }
    else {
	print "\nERROR: $pairnum is not correct. Must be 1 or 2. Exiting.\n";
	exit(1);
    }

    my $fh = get_fh($infile);
    open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

    my @aux = undef;
    my ($name, $comm, $seq, $qual);

    while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
	$seq = uc($seq) if $uppercase;
	if (defined $qual) {
	    print $out join "\n", "@".$name.$pair, $seq, '+', $qual, "\n";
	}
	else {
	    print $out join "\n", ">".$name.$pair, $seq, "\n";
	}
    }
    close $fh;
    close $out;

    compress($compress, $outfile) if $compress;
    exit;
}

sub make_pairs_and_singles {
    my ($fread, $rread, $fpread, $rpread, $fsread, $rsread, $index, $compress, $stats) = @_;

    my ($rseqpairs, $rct) = store_pair($rread);

    my @faux = undef;
    my ($fname, $fcomm, $fseq, $fqual, $forw_id, $rev_id, $fname_enc);
    my ($fct, $fpct, $rpct, $pct, $fsct, $rsct, $sct) = (0, 0, 0, 0, 0, 0, 0);

    my $fh = get_fh($fread);
    open my $fp, '>', $fpread or die "\nERROR: Could not open file: $fpread\n";
    open my $rp, '>', $rpread or die "\nERROR: Could not open file: $rpread\n";
    open my $fs, '>', $fsread or die "\nERROR: Could not open file: $fsread\n";

    while (($fname, $fcomm, $fseq, $fqual) = readfq(\*$fh, \@faux)) {
	$fct++;
	if ($fname =~ /(\/\d)$/) {
	    $fname =~ s/$1//;
	}
	elsif (defined $fcomm && $fcomm =~ /^\d/) {
	    $fcomm =~ s/^\d//;
	    $fname = mk_key($fname, $fcomm);
	}
	else {
	    print "\nERROR: Could not determine FASTA/Q format. ".
		"Please see https://github.com/sestaton/Pairfq or the README for supported formats. Exiting.\n\n";
	    exit(1);
	}
	
	if ($fname =~ /\|\|/) {
	    my ($name, $comm);
	    ($name, $comm) = mk_vec($fname);
	    $forw_id = $name.q{ 1}.$comm;
	    $rev_id  = $name.q{ 2}.$comm;
	}

	if (exists $rseqpairs->{$fname}) {
	    $fpct++; $rpct++;
	    if (defined $fqual) {
		my ($rread, $rqual) = mk_vec($rseqpairs->{$fname});
		if ($fname =~ /\|\|/) {
		    print $fp join "\n","@".$forw_id, $fseq, "+", $fqual, "\n";
		    print $rp join "\n","@".$rev_id, $rread, "+", $rqual, "\n";
		} 
		else {
		    print $fp join "\n","@".$fname.q{/1}, $fseq, "+", $fqual, "\n";
		    print $rp join "\n","@".$fname.q{/2}, $rread, "+", $rqual, "\n";
		}
	    } 
	    else {
		if ($fname =~ /\|\|/) {
		    print $fp join "\n",">".$forw_id, $fseq, "\n";
		    print $rp join "\n",">".$rev_id, $rseqpairs->{$fname}, "\n";
		} 
		else {
		    print $fp join "\n",">".$fname.q{/1}, $fseq, "\n";
		    print $rp join "\n",">".$fname.q{/2}, $rseqpairs->{$fname}, "\n";
		}
	    }
	    delete $rseqpairs->{$fname};
	} 
	else {
	    $fsct++;
	    if (defined $fqual) {
		if ($fname =~ /\|\|/) {
		    print $fs join "\n","@".$forw_id, $fseq, "+", $fqual, "\n";
		} 
		else {
		    print $fs join "\n","@".$fname.q{/1}, $fseq, "+", $fqual, "\n";
		}
	    } 
	    else {
		if ($fname =~ /\|\|/) {
		    print $fs join "\n",">".$forw_id, $fseq, "\n";
		} 
		else {
		    print $fs join "\n",">".$fname.q{/1}, $fseq, "\n";
		}
	    }
	}
	undef $forw_id;
	undef $rev_id;
    }
    close $fh;
    close $fp;
    close $rp;
    close $fs;

    open my $rs, '>', $rsread or die "\nERROR: Could not open file: $rsread\n";

    while (my ($rname_up_unenc, $rseq_up_unenc) = each %$rseqpairs) {
	$rsct++;
	my ($rname_up, $rcomm_up) = mk_vec($rname_up_unenc);
	my ($rseq_up, $rqual_up) = mk_vec($rseq_up_unenc);

	my $rev_id_up .= $rname_up.q{ 2}.$rcomm_up if defined $rcomm_up;
    
	if (defined $rcomm_up && defined $rqual_up) {
	    print $rs join "\n","@".$rev_id_up, $rseq_up, "+", $rqual_up, "\n";
	} 
	elsif (defined $rcomm_up && !defined $rqual_up) {
	    print $rs join "\n",">".$rev_id_up, $rseq_up_unenc, "\n";
	} 
	elsif (!defined $rcomm_up && defined $rqual_up) {
	    print $rs join "\n", "@".$rname_up.q{/2}, $rseq_up, "+", $rqual_up, "\n";
	}
	else {
	    print $rs join "\n",">".$rname_up.q{/2}, $rseq_up_unenc, "\n";
	}
    }
    close $rs;

    compress($compress, $fpread, $rpread, $fsread, $rsread) if $compress;
    $pct = $fpct + $rpct;
    $sct = $fsct + $rsct;

    if (defined $stats) {
	my $maxfn = max(length($fread), length($rread), length($fpread), length($rpread), length($fsread), length($rsread));
	my $offset = $maxfn + 38;
	my $date = qx(date); chomp $date;
	my $prog = basename($0, ());
	print "========= $prog version : $VERSION (completion time: $date)";
	printf "%-${offset}s %s %10d\n", "Total forward reads ($fread)", ":",$fct;
	printf "%-${offset}s %s %10d\n", "Total reverse reads ($rread)", ":", $rct;
	printf "%-${offset}s %s %10d\n", "Total forward paired reads ($fpread)", ":", $fpct;
	printf "%-${offset}s %s %10d\n", "Total reverse paired reads ($rpread)", ":", $rpct;
	printf "%-${offset}s %s %10d\n", "Total forward unpaired reads ($fsread)", ":", $fsct;
	printf "%-${offset}s %s %10d\n\n", "Total reverse unpaired reads ($rsread)", ":", $rsct;
	printf "%-${offset}s %s %10d\n", "Total paired reads", ":",  $pct;
	printf "%-${offset}s %s %10d\n", "Total unpaired reads", ":", $sct;
    }
    exit;
}

sub pairs_to_interleaved {
    my ($forward, $reverse, $outfile, $compress) = @_;

    my ($pairs, $ct) = store_pair($forward);

    my $fh = get_fh($reverse);
    open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

    my @raux = undef;
    my ($rname, $rcomm, $rseq, $rqual, $forw_id, $rev_id, $rname_enc);

    while (($rname, $rcomm, $rseq, $rqual) = readfq(\*$fh, \@raux)) {
	if ($rname =~ /(\/\d)$/) {
	    $rname =~ s/$1//;
	}
	elsif (defined $rcomm && $rcomm =~ /^\d/) {
	    $rcomm =~ s/^\d//;
	    $rname = mk_key($rname, $rcomm);
	}
	else {
	    print "\nERROR: Could not determine FastA/Q format. ".
		"Please see https://github.com/sestaton/Pairfq or the README for supported formats. Exiting.\n\n";
	    exit(1);
	}

	if ($rname =~ /\|\|/) {
	    my ($name, $comm) = mk_vec($rname);
	    $forw_id = $name.q{ 1}.$comm;
	    $rev_id  = $name.q{ 2}.$comm;
	}

	if (exists $pairs->{$rname}) {
	    if (defined $rqual) {
		my ($seqf, $qualf) = mk_vec($pairs->{$rname});
		if ($rname =~ /\|\|/) {
		    print $out join "\n", "@".$forw_id, $seqf, "+", $qualf, "\n";
		    print $out join "\n", "@".$rev_id, $rseq, "+", $rqual, "\n";
		}
		else {
		    print $out join "\n", "@".$rname.q{/1}, $seqf, "+", $qualf, "\n";
		    print $out join "\n", "@".$rname.q{/2}, $rseq, "+", $rqual, "\n";
		}
	    }
	    else {
		if ($rname =~ /\|\|/) {
		    print $out join "\n", ">".$forw_id, $pairs->{$rname}, "\n";
		    print $out join "\n", ">".$rev_id, $rseq, "\n";
		}
		else {
		    print $out join "\n", ">".$rname.q{/1}, $pairs->{$rname}, "\n";
		    print $out join "\n", ">".$rname.q{/2}, $rseq, "\n";                                               
		}
	    }
	}
    }
    close $fh;
    close $out;

    compress($compress, $outfile) if $compress;
    exit;
}

sub interleaved_to_pairs {
    my ($infile, $forward, $reverse, $compress) = @_;

    my $fh = get_fh($infile);
    open my $f, '>', $forward or die "\nERROR: Could not open file: $forward\n";
    open my $r, '>', $reverse or die "\nERROR: Could not open file: $reverse\n"; 

    my @aux = undef;
    my ($name, $comm, $seq, $qual);

    while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
	if (defined $comm && $comm =~ /^1/ || $name =~ /\/1$/) {
	    print $f join "\n", ">".$name, $seq, "\n" if !defined $qual && !defined $comm;
	    print $f join "\n", ">".$name.q{ }.$comm, $seq, "\n"if !defined $qual && defined $comm;
	    print $f join "\n", "@".$name, $seq, '+', $qual, "\n" if defined $qual && !defined $comm;
	    print $f join "\n", "@".$name.q{ }.$comm, $seq, '+', $qual, "\n" if defined $qual && defined $comm;
	}
	elsif (defined $comm && $comm =~ /^2/ || $name =~ /\/2$/) {
	    print $r join "\n", ">".$name, $seq, "\n" if !defined $qual && !defined $comm;
	    print $r join "\n", ">".$name.q{ }.$comm, $seq, "\n" if !defined $qual && defined $comm;
	    print $r join "\n", "@".$name, $seq, "+", $qual, "\n" if defined $qual && !defined $comm;
	    print $r join "\n", "@".$name.q{ }.$comm, $seq, "+", $qual, "\n" if defined $qual && defined $comm;
	}
    }
    close $fh;
    close $f;
    close $r;

    compress($compress, $forward, $reverse) if $compress;
    exit;
}

sub get_fh {
    my ($file) = @_;

    my $fh;
    if ($file =~ /\.gz$/) {
        open $fh, '-|', 'zcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /\.bz2$/) {
        open $fh, '-|', 'bzcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /-|STDIN/) {
	open $fh, '< -' or die "\nERROR: Could not open STDIN\n";
    }
    else {
        open $fh, '<', $file or die "\nERROR: Could not open file: $file\n";
    }

    return $fh;
}

sub compress {
    my ($compress, $fp, $rp, $fs, $rs) = @_;

    if ($compress eq 'gzip') {
	my $fpo = $fp.".gz";
	my $rpo = $rp.".gz" if defined $rp;
	my $fso = $fs.".gz" if defined $fs;
	my $rso = $rs.".gz" if defined $rs;
	
	gzip $fp => $fpo or die "gzip failed: $GzipError\n";
	gzip $rp => $rpo or die "gzip failed: $GzipError\n" if defined $rp;
	gzip $fs => $fso or die "gzip failed: $GzipError\n" if defined $fs;
        gzip $rs => $rso or die "gzip failed: $GzipError\n" if defined $rs;

	unlink $fp;
	unlink $rp if defined $rp;
	unlink $fs if defined $fs;
	unlink $rs if defined $rs;
    }
    elsif ($compress eq 'bzip2') {
	my $fpo = $fp.".bz2";
	my $rpo = $rp.".bz2" if defined $rp;
	my $fso = $fs.".bz2" if defined $fs;
	my $rso = $rs.".bz2" if defined $rs;

	bzip2 $fp => $fpo or die "bzip2 failed: $Bzip2Error\n";
	bzip2 $rp => $rpo or die "bzip2 failed: $Bzip2Error\n" if defined $rp;
	bzip2 $fs => $fso or die "bzip2 failed: $Bzip2Error\n" if defined $fs;
	bzip2 $rs => $rso or die "bzip2 failed: $Bzip2Error\n" if defined $rs;

	unlink $fp; 
	unlink $rp if defined $rp;
	unlink $fs if defined $fs;
	unlink $rs if defined $rs;
    }
    return;
}

sub store_pair {
    my ($file) = @_;

    my $rct = 0;
    my %rseqpairs;
    my $cwd = getcwd();
    my $db_file = File::Temp->new( TEMPLATE => "pairfq_XXXX",
				   DIR      => $cwd,
				   SUFFIX   => ".bdb",
				   UNLINK   => 0 );
    
    my @raux = undef;
    my ($rname, $rcomm, $rseq, $rqual, $rname_k, $rname_enc);

    my $fh = get_fh($file);

    while (($rname, $rcomm, $rseq, $rqual) = readfq(\*$fh, \@raux)) {
	$rct++;
	if ($rname =~ /(\/\d)$/) {
	    $rname =~ s/$1//;
	}
	elsif (defined $rcomm && $rcomm =~ /^\d/) {
	    $rcomm =~ s/^\d//;
	    $rname = mk_key($rname, $rcomm);
	}
	else {
	    print "\nERROR: Could not determine FASTA/Q format. ".
		"Please see https://github.com/sestaton/Pairfq or the README for supported formats. Exiting.\n\n";
	    exit(1);
	}
	
	$rseqpairs{$rname} = mk_key($rseq, $rqual) if defined $rqual;
	$rseqpairs{$rname} = $rseq if !defined $rqual;
    }
    close $fh;
    return (\%rseqpairs, $rct);
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

sub mk_key { return join "||", @_ }

sub mk_vec { return split /\|\|/, shift }

sub usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script [-h] [-m] [--version]

Required:
    addinfo           :      Add the pair info back to the FASTA/Q header.
    makepairs         :      Pair the forward and reverse reads and write singletons 
                             for both forward and reverse reads to separate files.
    joinpairs         :      Interleave the paired forward and reverse files.
    splitpairs        :      Split the interleaved file into separate files for the 
                             forward and reverse reads.

Options:
    --version         :       Print the program version and exit.
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}

sub makepairs_usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script makepairs [-f] [-r] [-fp] [-rp] [-fs] [-rs] [-im] [-h] [-m]

Required:
    -f|forward        :       File of foward reads (usually with "/1" or " 1" in the header).
    -r|reverse        :       File of reverse reads (usually with "/2" or " 2" in the header).
    -fp|forw_paired   :       Name for the file of paired forward reads.
    -rp|rev_paired    :       Name for the file of paired reverse reads.
    -fs|forw_unpaired :       Name for the file of singleton forward reads.
    -rs|rev_unpaired  :       Name for the file of singleton reverse reads.

Options:
    -idx|index        :       Construct an index for limiting memory usage.
                              NB: This may result in long run times for a large number of sequences. 
    -c|compress       :       Compress the output files. Options are 'gzip' or 'bzip2' (Default: No).
    -s|stats          :       Print statistics on the pairing results to STDOUT (Default: No).
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}

sub addinfo_usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script addinfo [-i] [-o] [-p] [-uc] [-h] [-m]

Required:
    -i|infile        :       The file of sequences without the pair information in the sequence name.
    -o|outfile       :       The file of sequences that will contain the sequence names with the pair information.
    -p|pairnum       :       The number to append to the sequence name. Integer (Must be 1 or 2).

Options:
    -c|compress      :       Compress the output files. Options are 'gzip' or 'bzip2' (Default: No).
    -uc|uppercase    :       Convert the sequence to uppercase.
    -h|help          :       Print a usage statement.
    -m|man           :       Print the full documentation.

EOF
}

sub splitpairs_usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script splitpairs [-i] [-f] [-r] [-c] [-h] [-m]

Required:
    -i|infile         :       File of interleaved forward and reverse reads.
    -f|forward        :       File to place the foward reads.
    -r|reverse        :       File to place the reverse reads.

Options:
    -c|compress       :       Compress the output files. Options are 'gzip' or 'bzip2' (Default: No).
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}

sub joinpairs_usage {
    my $script = basename($0);
    print STDERR<<EOF
USAGE: $script joinpairs [-f] [-r] [-o] [-im] [-h] [-m]

Required:
    -f|forward        :       File of foward reads (usually with "/1" or " 1" in the header).
    -r|reverse        :       File of reverse reads (usually with "/2" or " 2" in the header).
    -o|outfile        :       File of interleaved reads.

Options:
    -idx|index        :       Construct an index for limiting memory usage.
                              NB: This may result in long run times for a large number of sequences.
    -c|compress       :       Compress the output files. Options are 'gzip' or 'bzip2' (Default: No).
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}

__END__

=head1 NAME 
                                                                       
pairfq_lite.pl - Sync paired-end sequences from separate FASTA/Q files

=head1 SYNOPSIS    

## Add pair information back to the reads

pairfq_lite.pl addinfo -i s_1_1_trim.fq -o s_1_1_trim_info.fq -p 1

pairfq_lite.pl addinfo -i s_1_2_trim.fq -o s_1_2_trim_info.fq -p 2

## Sync paired-end reads and write singleton reads to separate files

pairfq_lite.pl makepairs -f s_1_1_trim_info.fq -r s_1_2_trim_info.fq -fp s_1_1_trim_paired.fq -rp s_1_2_trim_paired.fq -fs s_1_1_trim_unpaired.fq -rs s_1_2_trim_unpaired.fq --stats

## Interleave the paired-end reads

pairfq_lite.pl joinpairs -f s_1_1_trim_paired.fq -r s_1_2_trim_paired.fq -o s_1_interl.fq -im

## Split the interleaved reads into separate forward and reverse files

pairfq_lite.pl splitpairs -i s_1_interl.fq -f s_1_1_trim_p.fq -r s_1_2_trim_p.fq

=head1 DESCRIPTION
     
Re-pair paired-end sequences that may have been separated by quality trimming.
This script also writes the unpaired forward and reverse sequences to separate 
files so that they may be used for assembly or mapping. The input may be FastA
or FastQ format in either Illumina 1.3+ or Illumina 1.8 format. The input files
may be compressed with gzip or bzip2. Optionally, the script can interleave paired
files, separate interleaved files into separate forward and reverse files, and 
fix paired-end files which have lost the pair information. 

=head1 DEPENDENCIES

There are no external dependencies with the 'pairfq_lite.pl' script. See below for 
information on which Perls have been tested.

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
Perl 5.16.0 (Red Hat Enterprise Linux Server release 5.9 (Tikanga)); BerkeleyDB 0.54

=item *
Perl 5.18.0 (Red Hat Enterprise Linux Server release 5.9 (Tikanga)); BerkeleyDB 0.54

=back

=head1 AUTHOR 

S. Evan Staton                                                

=head1 CONTACT
 
statonse at gmail dot com

=head1 REQUIRED ARGUMENTS

=over 2

=item The method to perform. Must be one of: 'addinfo', 'makepairs', 'joinpairs', or 'splitpairs'.

  addinfo    | Add the pair info back to the FASTA/Q header.
  makepairs  | Pair the forward and reverse reads and write singletons for both forward and reverse reads to separate files.
  joinpairs  | Interleave the paired forward and reverse files.
  splitpairs | Split the interleaved file into separate files for the forward and reverse reads.

=back

=head1 OPTIONS

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

=item -idx, --index

  The computation should be not be in memory but on disk. This will be much slower, but it will not use hardly any
  RAM, even if there are many millions of sequences in each input file.

=item -c, --compress

  The output files should be compressed. If given, this option must be given the arguments 'gzip' to compress with gzip,
  or 'bzip2' to compress with bzip2.

=item -uc, --uppercase

 For the 'addinfo' method, uppercase the sequence.

=item -s, --stats

 For the 'makepairs' method, print (to STDOUT) statistics for paired/unpaired forward and reverse reads. This is useful for
 record keeping and debugging. The reason this is not the default is that people may want to run multiple instances of this
 command and redirect the output to the same file or to another program.

=item --version

 Get the program version and exit.

=item -h, --help

  Print a usage statement. 

=item -m, --man
  
  Print the full documentation.

=back

=cut
