use strict;
use warnings FATAL => 'all';
use ExtUtils::MakeMaker;

BEGIN {
    system("chmod +x bin/pairfq");
}

WriteMakefile(
    NAME             => 'bin/pairfq',
    AUTHOR           => q{S. Evan Staton <statonse@gmail.com>},
    VERSION_FROM     => 'bin/pairfq',
    LICENSE          => 'MIT',
    PL_FILES         => {},
    MIN_PERL_VERSION => 5.010,
    CONFIGURE_REQUIRES => {
        'ExtUtils::MakeMaker' => 0,
    },
    BUILD_REQUIRES => {
        'Test::More' => 0,
    },
    PREREQ_PM => {
	'BerkeleyDB'          => 0.54,       
        'IPC::System::Simple' => 1.21,
    },
    INST_SCRIPT => 'bin',
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    test  => { TESTS => 't/*.t' },
);

