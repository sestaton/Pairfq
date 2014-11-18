use 5.010;
use strict;
use warnings FATAL => 'all';
use File::Spec;
use File::Find;
use File::Basename;
use Test::More tests => 1;

my $cmd = File::Spec->catfile('bin', 'pairfq');
ok(-x $cmd, 'Can execute pairfq');
