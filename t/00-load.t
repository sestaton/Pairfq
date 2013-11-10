#!/usr/bin/env perl

use 5.012;
use strict;
use warnings FATAL => 'all';
use IPC::System::Simple qw(capture);
use Test::More tests => 2;

my $pairfq = capture('which pairfq');
chomp $pairfq;
ok($pairfq =~ qr/pairfq$/, 'Can locate pairfq');
ok(-x $pairfq, 'Can execute pairfq');

diag( "Testing Pairfq::VERSION, Perl $], $^X" );
