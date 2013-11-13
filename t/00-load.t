#!/usr/bin/env perl

use 5.012;
use strict;
use warnings FATAL => 'all';
use Test::More tests => 1;

my $pairfq = "bin/pairfq";
ok(-x $pairfq, 'Can execute pairfq');

diag( "Testing Pairfq::VERSION, Perl $], $^X" );
