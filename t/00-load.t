#!/usr/bin/env perl

use 5.012;
use strict;
use warnings FATAL => 'all';
use IPC::System::Simple qw(capture);
use Test::More tests => 1;

my $pairfq = "bin/pairfq";
ok(-x $pairfq, 'Can execute pairfq');

my $ver = capture([0..5], "bin/pairfq --version");
diag( "Testing Pairfq $ver, Perl $], $^X" );
