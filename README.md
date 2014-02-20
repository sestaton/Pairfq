Pairfq
======

Sync paired-end FASTA/Q files and keep singleton reads

[![Build Status](https://travis-ci.org/sestaton/Pairfq.png?branch=master)](https://travis-ci.org/sestaton/Pairfq)

**DEPENDENCIES**

There is an external package used by Pairfq (called BerkeleyDB) and it can be installed through your package manager, though it may be installed already. See the [Installing dependencies](https://github.com/sestaton/Pairfq/wiki/Installing-dependencies) page for suggestions on how to install for common operating systems.

**INSTALLATION**

Perl version 5.12 (or greater) must be installed to use Pairfq, and there are a couple of external modules required (which may be installed already depending on your version of Perl). If you have [cpanminus](http://search.cpan.org/~miyagawa/App-cpanminus-1.6935/lib/App/cpanminus.pm), installation can be done with a single command:

    cpanm git://github.com/sestaton/Pairfq.git

Alternatively, download the latest [release](https://github.com/sestaton/Pairfq/releases) and run the following commands in the top directory:

    perl Makefile.PL
    make
    make test
    make install

**TYPICAL USAGE CASES**

See the [Pairfq wiki](https://github.com/sestaton/Pairfq/wiki) for examples with each method.

**SUPPORT AND DOCUMENTATION**

After installation, you can find documentation for Pairfq with the `perldoc` command.

    perldoc pairfq

The documentation can also be accessed by specifying the manual option with `pairfq -m` or `pairfq --man`. The `pairfq` program will also print a diagnostic help message when executed with no arguments.

**ISSUES**

Report any issues at the Pairfq issue tracker: https://github.com/sestaton/Pairfq/issues

Be aware that Pairfq will not work for every data set given the wide range of FASTA/Q formats. Feel free to fork the project and modify the code to your needs, or submit code changes if you find any bugs. 

**ATTRIBUTION**

This project uses the [readfq](https://github.com/lh3/readfq) library written by Heng Li. The readfq code has been modified for error handling and to parse the comment line in the Casava header.

**LICENSE**

The MIT License should included with the project. If not, it can be found at: http://opensource.org/licenses/mit-license.php

Copyright (C) 2013 S. Evan Staton

[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/sestaton/pairfq/trend.png)](https://bitdeli.com/free "Bitdeli Badge")
