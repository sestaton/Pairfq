use 5.010;
use strict;
use warnings FATAL => 'all';
use File::Spec;
use Cwd qw(getcwd);
use Test::More tests => 10;

my $cmd        = File::Spec->catfile('blib', 'bin', 'pairfq');
my @addinfo    = qx($cmd addinfo 2>&1);
my @makepairs  = qx($cmd makepairs 2>&1);
my @splitpairs = qx($cmd splitpairs 2>&1);
my @joinpairs  = qx($cmd joinpairs 2>&1);
my @wrongtask  = qx($cmd jonpairs 2>&1);

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
