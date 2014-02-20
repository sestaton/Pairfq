#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use IPC::System::Simple qw(capture);
use Test::More tests => 8;

my @addinfo    = capture([0..5],"bin/pairfq addinfo 2>&1");
my @makepairs  = capture([0..5],"bin/pairfq makepairs 2>&1");
my @splitpairs = capture([0..5],"bin/pairfq splitpairs 2>&1");
my @joinpairs  = capture([0..5],"bin/pairfq joinpairs 2>&1");
my @wrongtask  = capture([0..5],"bin/pairfq jonpairs 2>&1");

## just look to see if the task was handled correctly to generate a usage statement
for my $out (@addinfo) {
    if ($out =~ /pairnum/i) {
	ok($out =~ qr/pairnum/i, 'Help menu generated correctly for pairnum task');
    }
}

for my $out (@makepairs) {
    if ($out =~ /forward/i) {
        ok($out =~ qr/forward/i, 'Help menu generated correctly for makepairs task');
    }
}

for my $out (@splitpairs) {
    if ($out =~ /forward/i) {
        ok($out =~ qr/forward/i, 'Help menu generated correctly for splitpairs task');
    }
}

for my $out (@joinpairs) {
    if ($out =~ /forward/i) {
        ok($out =~ qr/forward/i, 'Help menu generated correctly for joinpairs task');
    }
}

## look to see if a wrong argument to task is handled correctly
for my $out (@wrongtask) {
    if ($out =~ /^error/i) {
	ok($out =~ qr/\'jonpairs\' is not recognized/, 'Wrong task properly handled');
    }
}
