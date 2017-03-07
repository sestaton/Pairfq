#!/usr/bin/env perl

use strict;
use warnings;
use Cwd qw(getcwd abs_path);
use File::Basename;
use Getopt::Long;
use Pod::Usage;

our $VERSION = '0.17.0';

my $infile;     # input file for 'addinfo', 'splitpairs' and 'makepairs' methods
my $outfile;    # output file for 'addinfo' method
my $fread;      # file of forward reads for 'splitpairs', 'makepairs' and 'joinpairs' methods
my $rread;      # file of forward reads for 'splitpairs', 'makepairs' and 'joinpairs' methods
my $fpread;     # file of paired forward reads for 'makepairs' method
my $rpread;     # file of paired reverse reads for 'makepairs' method
my $fsread;     # file of unpaired forward reads for 'makepairs' method
my $rsread;     # file of unpaired reverse reads for 'makepairs' method 
my $pairnum;    # for the 'addinfo' method
my $uppercase;  # for 'addinfo' method
my $stats;      # currently, for 'makepairs' option only

my $version; 
my $help;
my $man;
my $script = basename($0, ());
$script = "pairfq_lite" if $script =~ /^-$|stdin/i;

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
	   'uc|uppercase'       => \$uppercase,
	   's|stats'            => \$stats,
	   'version'            => \$version,
	   'h|help'             => \$help,
	   'm|man'              => \$man,
	  ) or pod2usage( "Try '$0 --man' for more information." );

#
# Check @ARGV
#
usage($script) and exit(0) if $help;
pod2usage( -verbose => 2 ) if $man;
print $VERSION and exit(0) if $version;

my $method = shift;
if (!defined $method) {
    print "\nERROR: Command line not parsed correctly. Check input.\n\n";
    usage($script);
    exit(1);
}

if ($method eq 'addinfo') {
    if (!$pairnum || !$infile || !$outfile) {
	print "\nERROR: Command line not parsed correctly. Check input.\n\n";
	addinfo_usage($script);
	exit(1);
    }
    add_pair_info($pairnum, $infile, $outfile, $uppercase);
}
elsif ($method eq 'makepairs') {
    if ($infile && $fpread && $rpread && $fsread && $rsread) {
        interleaved_to_pairs_and_singles($script, $infile, $fpread, $rpread, $fsread, $rsread, $stats);
    }
    elsif (!$infile && $fread && $rread && $fpread && $rpread && $fsread && $rsread) {
        make_pairs_and_singles($script, $fread, $rread, $fpread, $rpread, $fsread, $rsread, $stats);
    }
    else {
        print "\nERROR: Command line not parsed correctly. Check input.\n\n";
        makepairs_usage($script);
        exit(1);
    }
}
elsif ($method eq 'joinpairs') {
    if (!$fread || !$rread || !$outfile) {
	print "\nERROR: Command line not parsed correctly. Check input.\n\n";
	joinpairs_usage($script);
	exit(1);
    }
    pairs_to_interleaved($fread, $rread, $outfile);
}
elsif ($method eq 'splitpairs') {
    if (!$infile || !$fread || !$rread) {
	print "\nERROR: Command line not parsed correctly. Check input.\n\n";
	splitpairs_usage($script);
	exit(1);
    }
    interleaved_to_pairs($infile, $fread, $rread);
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
    my ($pairnum, $infile, $outfile, $uppercase) = @_;

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

    my $fh  = get_fh($infile);
    my $out = get_outfh($outfile);

    my @aux = undef;
    my ($name, $comm, $seq, $qual);

    while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
	$seq = uc($seq) if $uppercase;
	print $out join "\n", "@".$name.$pair, $seq, "+", "$qual\n" if defined $qual;
	print $out join "\n", ">".$name.$pair, "$seq\n" if !defined $qual;
    }
    close $fh;
    close $out;

    exit;
}

sub make_pairs_and_singles {
    my ($script, $fread, $rread, $fpread, $rpread, $fsread, $rsread, $stats) = @_;

    my ($rseqpairs, $rct) = store_pair($rread);

    my @faux = undef;
    my ($fname, $fcomm, $fseq, $fqual, $forw_id, $rev_id);
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
	    $forw_id = $name." 1".$comm if defined $comm;
	    $forw_id = $name."/1" if !defined $comm;
	    $rev_id = $name." 2".$comm if defined $comm;
	    $rev_id = $name."/2" if !defined $comm;
	}

	if (exists $rseqpairs->{$fname}) {
	    $fpct++; 
	    $rpct++;
	    if (defined $fqual) {
		my ($rread, $rqual) = mk_vec($rseqpairs->{$fname});
		if ($fname =~ /\|\|/) {
		    print $fp join "\n", "@".$forw_id, $fseq, "+", "$fqual\n";
		    print $rp join "\n", "@".$rev_id, $rread, "+", "$rqual\n";
		} 
		else {
		    print $fp join "\n", "@".$fname."/1", $fseq, "+", "$fqual\n";
		    print $rp join "\n", "@".$fname."/2", $rread, "+", "$rqual\n";
		}
	    } 
	    else {
		if ($fname =~ /\|\|/) {
		    print $fp join "\n", ">".$forw_id, "$fseq\n";
		    print $rp join "\n", ">".$rev_id, "$rseqpairs->{$fname}\n";
		} 
		else {
		    print $fp join "\n", ">".$fname."/1", "$fseq\n";
		    print $rp join "\n", ">".$fname."/2", "$rseqpairs->{$fname}\n";
		}
	    }
	    delete $rseqpairs->{$fname};
	} 
	else {
	    $fsct++;
	    if (defined $fqual) {
		if ($fname =~ /\|\|/) {
		    print $fs join "\n", "@".$forw_id, $fseq, "+", "$fqual\n";
		} 
		else {
		    print $fs join "\n", "@".$fname."/1", $fseq, "+", "$fqual\n";
		}
	    } 
	    else {
		if ($fname =~ /\|\|/) {
		    print $fs join "\n", ">".$forw_id, "$fseq\n";
		} 
		else {
		    print $fs join "\n", ">".$fname."/1", "$fseq\n";
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
	my ($rseq_up, $rqual_up)  = mk_vec($rseq_up_unenc);

	my $rev_id_up .= $rname_up." 2".$rcomm_up if defined $rcomm_up;
    
	if (defined $rcomm_up && defined $rqual_up) {
	    print $rs join "\n", "@".$rev_id_up, $rseq_up, "+", "$rqual_up\n";
	} 
	elsif (defined $rcomm_up && !defined $rqual_up) {
	    print $rs join "\n", ">".$rev_id_up, "$rseq_up_unenc\n";
	} 
	elsif (!defined $rcomm_up && defined $rqual_up) {
	    print $rs join "\n", "@".$rname_up."/2", $rseq_up, "+", "$rqual_up\n";
	}
	else {
	    print $rs join "\n", ">".$rname_up."/2", "$rseq_up_unenc\n";
	}
    }
    close $rs;

    $pct = $fpct + $rpct;
    $sct = $fsct + $rsct;

    if (defined $stats) {
	my $maxfn = max(length($fread), length($rread), length($fpread), length($rpread), length($fsread), length($rsread));
	my $offset = $maxfn + 38;
	my $date = qx(date); chomp $date;
	print "========= $script version : $VERSION (completion time: $date)\n";
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

sub interleaved_to_pairs_and_singles {
    my ($script, $infile, $fpread, $rpread, $fsread, $rsread, $stats) = @_;

    my $fh = get_fh($infile);
    open my $fp, '>', $fpread or die "\nERROR: Could not open file: $fpread\n";
    open my $rp, '>', $rpread or die "\nERROR: Could not open file: $rpread\n";
    open my $fs, '>', $fsread or die "\nERROR: Could not open file: $fsread\n";
    open my $rs, '>', $rsread or die "\nERROR: Could not open file: $rsread\n";

    my @aux = undef;
    my ($fct, $rct, $fpct, $rpct, $pct, $fsct, $rsct, $sct, $pair) = (0, 0, 0, 0, 0, 0, 0, 0, 0);
    my %singles;
    my ($fpairname, $rpairname);
    my ($name, $comm, $seq, $qual);
    my ($fname, $fcomm, $fseq, $fqual, $rname, $rcomm, $rseq, $rqual);

    while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
        if ($name =~ /\/1$/) {
            $fct++;
            $pair = 1;
            ($fname, $fpairname, $fseq, $fqual) = ($name, $name, $seq, $qual);
            $fpairname =~ s/\/1//;
            $singles{$fname} = { name => $fname, seq => $fseq, qual => $fqual, pair => $pair };
        }
        elsif (defined $comm && $comm =~ /^1/) {
            $fct++;
            $pair = 1;
            ($fname, $fpairname, $fseq, $fcomm, $fqual) = ($name, $name, $seq, $comm, $qual);
            $singles{$fname} = { name => $fname, seq => $fseq, comm => $fcomm, qual => $fqual, pair => $pair };
        }
        elsif ($name =~ /\/2$/) {
            $rct++;
            $pair = 2;
            ($rname, $rpairname, $rseq, $rqual) = ($name, $name, $seq, $qual);
            $rpairname =~ s/\/2//;
            $singles{$rname} = { name => $rname, seq => $rseq, qual => $rqual, pair => $pair };
        }
        elsif (defined $comm && $comm =~ /^2/) {
            $rct++;
            $pair = 2;
            ($rname, $rpairname, $rseq, $rcomm, $rqual) = ($name, $name, $seq, $comm, $qual);
            $singles{$rname} = { name => $rname, seq => $rseq, comm => $rcomm, qual => $rqual, pair => $pair };
        }

	next unless defined $fpairname && defined $rpairname;
        if ($fpairname eq $rpairname) {
            $fpct++;
            $rpct++;
            say $fp join "\n", ">".$fname, $fseq 
                  if !defined $fqual && !defined $fcomm;
            say $fp join "\n", ">".$fname." ".$fcomm, $fseq 
                if !defined $fqual && defined $fcomm;
            say $fp join "\n", "@".$fname, $fseq, "+", $fqual 
                if defined $fqual && !defined $fcomm;
            say $fp join "\n", "@".$fname." ".$fcomm, $fseq, "+", $fqual 
                if defined $fqual && defined $fcomm;
            
            say $rp join "\n", ">".$rname, $rseq 
                if !defined $rqual && !defined $rcomm;
            say $rp join "\n", ">".$rname." ".$rcomm, $rseq 
                if !defined $rqual && defined $rcomm;
            say $rp join "\n", "@".$rname, $rseq, "+", $rqual 
                if defined $rqual && !defined $rcomm;
            say $rp join "\n", "@".$rname." ".$rcomm, $rseq, "+", $rqual 
                if defined $rqual && defined $rcomm;
            delete $singles{$fname};
            delete $singles{$rname};
        }
        $pair = 0;
    }
    close $fh;
    close $fp;
    close $rp;

    for my $id (keys %singles) {
        my $sfh;
        if ($singles{$id}->{'pair'} == 1) {
            $fsct++;
            $sfh = $fs;
        }
        else {
            $rsct++;
            $sfh = $rs;
        }

	say $sfh join "\n", ">".$singles{$id}->{'name'}, $singles{$id}->{'seq'} 
            if !defined $singles{$id}->{'qual'} && !defined $singles{$id}->{'comm'};
        say $sfh join "\n", ">".$singles{$id}->{'name'}." ".$singles{$id}->{'comm'}, $singles{$id}->{'seq'} 
            if !defined $singles{$id}->{'qual'} && defined $singles{$id}->{'comm'};
        say $sfh join "\n", "@".$singles{$id}->{'name'}, $singles{$id}->{'seq'}, "+", $singles{$id}->{'qual'} 
            if defined $singles{$id}->{'qual'} && !defined $singles{$id}->{'comm'};
        say $sfh join "\n", "@".$singles{$id}->{'name'}." ".$singles{$id}->{'comm'}, $singles{$id}->{'seq'}, "+", $singles{$id}->{'qual'} 
            if defined $singles{$id}->{'qual'} && defined $singles{$id}->{'comm'};
    }
    close $fs;
    close $rs;

    $pct = $fpct + $rpct;
    $sct = $fsct + $rsct;

    if (defined $stats) {
        my $maxfn = max(length($infile), length($fpread), length($rpread), length($fsread), length($rsread));
        my $offset = $maxfn + 38;
        my $date = qx(date); chomp $date;
        print "========= $script version : $VERSION (completion time: $date)\n";
        printf "%-${offset}s %s %10d\n", "Total forward reads ($infile)", ":",$fct;
        printf "%-${offset}s %s %10d\n", "Total reverse reads ($infile)", ":", $rct;
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
    my ($forward, $reverse, $outfile) = @_;

    my ($pairs, $ct) = store_pair($forward);

    my $fh  = get_fh($reverse);
    my $out = get_outfh($outfile);

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
	    $forw_id = $name." 1".$comm if defined $comm;
	    $forw_id = $name."/1" if !defined $comm;
	    $rev_id  = $name." 2".$comm if defined $comm;
	    $rev_id  = $name."/2" if !defined $comm;
	}

	if (exists $pairs->{$rname}) {
	    if (defined $rqual) {
		my ($seqf, $qualf) = mk_vec($pairs->{$rname});
		if ($rname =~ /\|\|/) {
		    print $out join "\n", "@".$forw_id, $seqf, "+", "$qualf\n";
		    print $out join "\n", "@".$rev_id, $rseq, "+", "$rqual\n";
		}
		else {
		    print $out join "\n", "@".$rname."/1", $seqf, "+", "$qualf\n";
		    print $out join "\n", "@".$rname."/2", $rseq, "+", "$rqual\n";
		}
	    }
	    else {
		if ($rname =~ /\|\|/) {
		    print $out join "\n", ">".$forw_id, "$pairs->{$rname}\n";
		    print $out join "\n", ">".$rev_id, "$rseq\n";
		}
		else {
		    print $out join "\n", ">".$rname."/1", "$pairs->{$rname}\n";
		    print $out join "\n", ">".$rname."/2", "$rseq\n";                                               
		}
	    }
	}
    }
    close $fh;
    close $out;

    exit;
}

sub interleaved_to_pairs {
    my ($infile, $forward, $reverse) = @_;

    my $fh = get_fh($infile);
    open my $f, '>', $forward or die "\nERROR: Could not open file: $forward\n";
    open my $r, '>', $reverse or die "\nERROR: Could not open file: $reverse\n"; 

    my @aux = undef;
    my ($name, $comm, $seq, $qual);

    while (($name, $comm, $seq, $qual) = readfq(\*$fh, \@aux)) {
	if (defined $comm && $comm =~ /^1/ || $name =~ /\/1$/) {
	    print $f join "\n", ">".$name, "$seq\n" if !defined $qual && !defined $comm;
	    print $f join "\n", ">".$name." ".$comm, "$seq\n" if !defined $qual && defined $comm;
	    print $f join "\n", "@".$name, $seq, "+", "$qual\n" if defined $qual && !defined $comm;
	    print $f join "\n", "@".$name." ".$comm, $seq, '+', "$qual\n" if defined $qual && defined $comm;
	}
	elsif (defined $comm && $comm =~ /^2/ || $name =~ /\/2$/) {
	    print $r join "\n", ">".$name, "$seq\n" if !defined $qual && !defined $comm;
	    print $r join "\n", ">".$name." ".$comm, "$seq\n" if !defined $qual && defined $comm;
	    print $r join "\n", "@".$name, $seq, "+", "$qual\n" if defined $qual && !defined $comm;
	    print $r join "\n", "@".$name." ".$comm, $seq, "+", "$qual\n" if defined $qual && defined $comm;
	}
    }
    close $fh;
    close $f;
    close $r;

    exit;
}

sub get_fh {
    my ($file) = @_;

    unless ($file =~ /^-$|STDIN/i) {
        $file = abs_path($file);
    }

    my $fh;
    if ($file =~ /\.gz$/) {
        open $fh, '-|', 'zcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /\.bz2$/) {
        open $fh, '-|', 'bzcat', $file or die "\nERROR: Could not open file: $file\n";
    }
    elsif ($file =~ /^-$|STDIN/i) {
	open $fh, '<&', \*STDIN or die "\nERROR: Could not open STDIN\n";
	    }
    else {
        open $fh, '<', $file or die "\nERROR: Could not open file: $file\n";
    }

    return $fh;
}

sub get_outfh {
    my ($file) = @_;

    unless ($file =~ /^-$|STDOUT/i) {
        $file = abs_path($file);
    }

    my $fh;
    if ($file =~ /^-$|STDOUT/i) {
	open $fh, '>&', \*STDOUT or die "\nERROR: Could not open STDOUT\n";
    }
    else {
	open $fh, '>', $file or die "\nERROR: Could not open file: $file\n";
    }

    return $fh;
}

sub store_pair {
    my ($file) = @_;

    my $rct = 0;
    my %rseqpairs;
    my $cwd = getcwd();
    
    my @raux = undef;
    my ($rname, $rcomm, $rseq, $rqual);

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

sub max {
    my $max = shift;
    for (@_) { $max = $_ if $_ > $max }
    return $max;
}

sub usage {
    my ($script) = @_;
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
    my ($script) = @_;
    print STDERR<<EOF
USAGE: $script makepairs [-i] [-f] [-r] [-fp] [-rp] [-fs] [-rs] [-s] [-h] [-m]

Required:
    -i|infile         :       File of interleaved forward and reverse reads that has been trimmed.
                              As below, the forward and reverse reads must be labeled in the name
                              or comment.
    -f|forward        :       File of foward reads (usually with "/1" or " 1" in the header).
    -r|reverse        :       File of reverse reads (usually with "/2" or " 2" in the header).
    -fp|forw_paired   :       Name for the file of paired forward reads.
    -rp|rev_paired    :       Name for the file of paired reverse reads.
    -fs|forw_unpaired :       Name for the file of singleton forward reads.
    -rs|rev_unpaired  :       Name for the file of singleton reverse reads.

Options:
    -s|stats          :       Print statistics on the pairing results to STDOUT (Default: No).
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}

sub addinfo_usage {
    my ($script) = @_;
    print STDERR<<EOF
USAGE: $script addinfo [-i] [-o] [-p] [-uc] [-h] [-m]

Required:
    -i|infile        :       The file of sequences without the pair information in the sequence name.
    -o|outfile       :       The file of sequences that will contain the sequence names with the pair information.
    -p|pairnum       :       The number to append to the sequence name. Integer (Must be 1 or 2).

Options:
    -uc|uppercase    :       Convert the sequence to uppercase.
    -h|help          :       Print a usage statement.
    -m|man           :       Print the full documentation.

EOF
}

sub splitpairs_usage {
    my ($script) = @_;
    print STDERR<<EOF
USAGE: $script splitpairs [-i] [-f] [-r] [-h] [-m]

Required:
    -i|infile         :       File of interleaved forward and reverse reads.
    -f|forward        :       File to place the foward reads.
    -r|reverse        :       File to place the reverse reads.

Options:
    -h|help           :       Print a usage statement.
    -m|man            :       Print the full documentation.

EOF
}

sub joinpairs_usage {
    my ($script) = @_;
    print STDERR<<EOF
USAGE: $script joinpairs [-f] [-r] [-o] [-h] [-m]

Required:
    -f|forward        :       File of foward reads (usually with "/1" or " 1" in the header).
    -r|reverse        :       File of reverse reads (usually with "/2" or " 2" in the header).
    -o|outfile        :       File of interleaved reads.

Options:
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

OR

pairfq_lite.pl makepairs -i s_interl_trimmed.fq -fp s_1_1_trim_paired.fq -rp s_1_2_trim_paired.fq -fs s_1_1_trim_unpaired.fq -rs s_1_2_trim_unpaired.fq --stats
## Interleave the paired-end reads

pairfq_lite.pl joinpairs -f s_1_1_trim_paired.fq -r s_1_2_trim_paired.fq -o s_1_interl.fq

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

Copyright (C) 2013-2016 S. Evan Staton
 
=head1 TESTED WITH:

=over

=item *
Perl 5.6.2 (Ubuntu 12.04.3 LTS)

=item *
Perl 5.8.9 (Ubuntu 12.04.3 LTS)

=item *
Perl 5.14.1 (Red Hat Enterprise Linux Server release 5.7 (Tikanga))

=item *
Perl 5.14.2 (Red Hat Enterprise Linux Desktop release 6.2 (Santiago); Fedora 17)

=item *
Perl 5.16.0 (Red Hat Enterprise Linux Server release 5.9 (Tikanga))

=item *
Perl 5.18.0 (Red Hat Enterprise Linux Server release 5.9 (Tikanga))

=item *
Perl 5.20.1 (Red Hat Enterprise Linux Server release 5.9 (Tikanga))

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

=item -i, --infile

  For the 'addinfo' method, this would be any FASTA/Q file (or STDIN). For the 'splitpairs' method,
  this would be either the forward or reverse file from a paired-end sequencing run. For the
  'makepairs' method, this would be the interleaved file of forward and reverse reads that
  has been trimmed.

=item -o, --outfile

  The outfile for the 'addinfo' or 'joinpairs' methods (may be STDOUT instead of a file).

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

=item -p, --pairnum

  The pair number to add to the file with the 'addinfo' method. Should be either '1' or '2' and other arguments
  with generate an exception.

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
