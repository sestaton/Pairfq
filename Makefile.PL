use strict;
use warnings FATAL => 'all';
use ExtUtils::MakeMaker;

if (eval ($ExtUtils::MakeMaker::VERSION) <= 6.55) {
    print "\n[ERROR]: ExtUtils::MakeMaker >= 6.5503 is required to build the package correctly.\n";
    print "The installed version is $ExtUtils::MakeMaker::VERSION.\n";
    print "Install the latest ExtUtils::MakeMaker and try again. Exiting.\n";
    exit(1);
}

BEGIN {
    system("chmod +x bin/pairfq");
}

WriteMakefile(
    NAME             => 'Pairfq',
    AUTHOR           => q{S. Evan Staton <statonse@gmail.com>},
    VERSION_FROM     => 'bin/pairfq',
    LICENSE          => 'MIT',
    PL_FILES         => {},
    MIN_PERL_VERSION => 5.010,
    CONFIGURE_REQUIRES => {
        'ExtUtils::MakeMaker' => 6.5503,
    },
    BUILD_REQUIRES => {
        'Test::More' => 0,
    },
    PREREQ_PM => {
	'DBD::SQLite'         => 1.44,       
	'Tie::Hash::DBD'      => 0.13,
        'IPC::System::Simple' => 1.21,
	'List::MoreUtils'     => 0.33,
    },
    INST_SCRIPT => 'bin',
    dist  => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    test  => { TESTS => 't/*.t' },
);

