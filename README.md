Pairfq
======

Sync paired-end FASTA/Q files and keep singleton reads

**INSTALLATION**

Perl version 5.12 (or greater) must be installed to use Pairfq, and there are a couple of external modules required (which may be installed already depending on your version of Perl). If you have [cpanminus](http://search.cpan.org/~miyagawa/App-cpanminus-1.6935/lib/App/cpanminus.pm), installation can be done with a single command:

    cpanm git://github.com/sestaton/Pairfq.git

**USAGE**

Type `pairfq` at the command line and you will see a menu describing the usage. 

    $ pairfq

    ERROR: Command line not parsed correctly. Check input.

    USAGE: pairfq [-h] [-m] [--version]

    Required:
        addinfo           :      Add the pair info back to the FASTA/Q header.
        makepairs         :      Pair the forward and reverse reads and write singletons 
                                 for both forward and reverse reads to separate files.
        joinpairs         :      Interleave the paired forward and reverse files.
        splitpairs        :      Split the interleaved file into separate files for the 
                                 forward and reverse reads.
    
    Options:
        --version         :       Print the program version and exit.
        -h|help           :       Print a usage statement.
        -m|man            :       Print the full documentation.

Specifying the method with no arguments will print the usage for that method. For example, 

    $ pairfq makepairs

    ERROR: Command line not parsed correctly. Check input.

    USAGE: pairfq makepairs [-f] [-r] [-fp] [-rp] [-fs] [-rs] [-im] [-h] [-m]

    Required:
        -f|forward        :       File of foward reads (usually with "/1" or " 1" in the header).
        -r|reverse        :       File of reverse reads (usually with "/2" or " 2" in the header).
        -fp|forw_paired   :       Name for the file of paired forward reads.
        -rp|rev_paired    :       Name for the file of paired reverse reads.
        -fs|forw_unpaired :       Name for the file of singleton forward reads.
        -rs|rev_unpaired  :       Name for the file of singleton reverse reads.

    Options:
        -idx|index        :       Construct an index for limiting memory usage.
                                  NB: This may result in long run times for a large number of sequences. 
        -c|compress       :       Compress the output files. Options are 'gzip' or 'bzip2' (Default: No).
        -h|help           :       Print a usage statement.
        -m|man            :       Print the full documentation.

Running the command `pairfq -m` will print the full documentation.

**EXPECTED FORMATS**

The input should be in [FASTA](http://en.wikipedia.org/wiki/FASTA_format) or [FASTQ](http://en.wikipedia.org/wiki/FASTQ_format) format. It is fine if the input files are compressed (with either gzip or bzip2).

Currently, data from the Casava pipeline version 1.4 are supported. For example,

    @HWUSI-EAS100R:6:73:941:1973#0/1

As well as Casava 1.8+ format,

    @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG

The overall format of the sequence name and comment may vary, but there must be an integer (1 or 2) at the end of the sequence name or as the first character in the comment (following a space after the sequence name). If your data is missing this pair information it will be necessary to fix them first (with the `addinfo` method, see below).

**METHODS**

Pairfq has several different methods which can be executed. Below is a brief description of each.

* **makepairs**

  * Pair the forward and reverse reads and write the singletons to separate files.

* **joinpairs**

  * Interleave the paired reads for assembly or mapping.

* **splitpairs**

  * Separate the interleaved FASTA/Q file into separate files for the forward and reverse reads.

* **addinfo**

  * Add the pair information back to the data. After filtering or sampling Casava 1.8+ data, the pair information is often lost, making downstream analyses difficult. For example, `@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG` usually becomes `@EAS139:136:FC706VJ:2:2104:15343:197393`. This script will add the pair information back (to become `@EAS139:136:FC706VJ:2:2104:15343:197393/1`). There is no way to know what was in the comment, so it will not be restored. 

**TYPICAL USAGE CASES**

See the [Pairfq wiki](https://github.com/sestaton/Pairfq/wiki) for examples with each method.

**ISSUES**

Report any issues at: https://github.com/sestaton/Pairfq/issues

Be aware that Pairfq will not work for every data set given the wide range of FASTA/Q formats. Feel free to fork the project and modify the code to your needs, or submit code changes if you find any bugs. 

**ATTRIBUTION**

This project uses the [readfq](https://github.com/lh3/readfq) library written by Heng Li. The readfq code has been modified for error handling and to parse the comment line in the Casava header.

**LICENSE**

The MIT License should included with the project. If not, it can be found at: http://opensource.org/licenses/mit-license.php

Copyright (C) 2013 S. Evan Staton

[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/sestaton/pairfq/trend.png)](https://bitdeli.com/free "Bitdeli Badge")
