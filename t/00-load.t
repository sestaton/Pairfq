#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use File::Spec;
use IPC::System::Simple qw(capture);
use Test::More tests => 1;

my $cmd = File::Spec->catfile('bin', 'pairfq');
ok(-x $cmd, 'Can execute pairfq');

my $ver = capture([0..5], "$cmd --version");
chomp $ver;
diag( "Testing Pairfq v$ver, Perl $], $^X" );
