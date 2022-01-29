#!perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use Config;
use File::Spec;
use File::Find;
use File::Basename;

use Test::More tests => 2;

my $cmd = File::Spec->catfile('blib', 'bin', 'pairfq');
if ($Config{osname} =~ /Win32/i) {
    # https://perldoc.perl.org/perlport#-X
    # On Windows, -x looks at the suffix
    ok(-e $cmd, 'Can execute pairfq');  
}
else {
    ok(-x $cmd, 'Can execute pairfq');
}

my $vers = qx($cmd --version);
chomp $vers;
ok( $vers =~ /\d+\.\d+\.\d+/, 'Can get pairfq version' );
diag( "Testing Pairfq $vers, Perl $], $^X" );
