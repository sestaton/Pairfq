Pairfq
======

Sync paired-end FastA/Q files and keep orphaned reads

**INSTALLATION**

You can specify the full path to the script, or move to a working directory and use it in place. Perl must be installed to use Pairfq, but there are no external modules required.

**USAGE**

Just type `perl pairfq.pl` and you will see a menu describing the usage. 

    $ perl pairfq.pl

    ERROR: Command line not parsed correctly. Check input.

    USAGE: pairfq.pl [-f] [-r] [-fp] [-rp] [-fs] [-rs] [-im] [-h] [-im]

    Required:
        -f|forward        :       File of foward reads (usually with "/1" or " 1" in the header).
        -r|reverse        :       File of reverse reads (usually with "/2" or " 2" in the header).
        -fp|forw_paired   :       Name for the file of paired forward reads.
        -rp|rev_paired    :       Name for the file of paired reverse reads.
        -fs|forw_unpaired :       Name for the file of singleton forward reads.
        -rs|rev_unpaired  :       Name for the file of singleton reverse reads.

    Options:
        -im|in_memory     :       Construct a database in memory for faster execution.
         -h|help           :       Print a usage statement.
         -m|man            :       Print the full documentation.

Running the command `perl pairfq.pl -m` will print the full documentation.

**EXPECTED FORMATS**

The input will be two files (i.e., forward and reverse) in [FASTA](http://en.wikipedia.org/wiki/FASTA_format) or [FASTQ](http://en.wikipedia.org/wiki/FASTQ_format) format that are expected to have reads out of order due to quality trimming.

Currently, data from the Casava pipeline version 1.4 are supported. For example,

    @HWUSI-EAS100R:6:73:941:1973#0/1

As well Casava 1.8+ format,

    @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG

The overall format of the sequence name and comment may vary, but there must be an integer (1 or 2) at the end of the sequence name or as the first character in the comment (following a space after the sequence name). If your data is missing this pair information it will be necessary to fix them first (see below). 

**UTILITIES**

In the Pairfq/utils subdirectory are several stand-alone scripts for working with paired-end FastA/Q files. Briefly, the scripts (and their functions) included are:

* **pairs_to_interleaved.pl**

  * Interleave the paired reads for assembly.

* **interleaved_to_pairs.pl**

  * Separate the interleaved FastA/Q file into separate files for the forward and reverse reads.

* **add_pair_info.pl**

  * Add the pair information back to the data. After filtering or sampling Casava 1.8+ data, the pair information is often lost, making downstream analyses difficult. For example, `@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG` usually becomes `@EAS139:136:FC706VJ:2:2104:15343:197393`. This script will add the pair information back (to become `@EAS139:136:FC706VJ:2:2104:15343:197393/`). There is no way to know what was in the comment, so will not be restored. 

As with `pairfq.pl`, typing the name of the script will print the usage statement, while typing the name of the script followed by `-m` or `--man` will print the full documentation.

**ISSUES**

Report any issues at: https://github.com/sestaton/Pairfq/issues

Be aware that Pairfq will not work for every data set given the wide range of FastA/Q formats. Feel free to fork the project and modify the code to your needs, or submit code changes if you find any bugs. 

**ATTRIBUTION**

This project would not be possible without the excellent [readfq](https://github.com/lh3/readfq) library written by Heng Li. The code used by Pairfq has been modified for error handling and to parse the comment line in the Casava header.

**LICENSE**

The MIT License should included with the project. If not, it can be found at: http://opensource.org/licenses/mit-license.php

Copyright (C) 2013 S. Evan Staton