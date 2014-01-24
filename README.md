Pairfq
======

Sync paired-end FASTA/Q files and keep singleton reads

**INSTALLATION**

Perl version 5.12 (or greater) must be installed to use Pairfq, but there are no external modules required. If you have [cpanminus](http://search.cpan.org/~miyagawa/App-cpanminus-1.6935/lib/App/cpanminus.pm), installation can be done in one shot:

    cpanm -n git://github.com/sestaton/Pairfq.git

If you don't have cpanminus, it is advised to run these commands in the top directory (requires `make`):

    perl Makefile.PL
    make
    make install

Updating is then as easy as repeating those commands (after a git pull or clone). 

As a third option, if you don't have `make`, or you want to move `pairfq` to a custom location, just make the program (found in the "bin" directory) executable and add it to your PATH. The problem with this third approach is that it is harder to keep updated.

**USAGE**

Type `pairfq` at the command line and you will see a menu describing the usage. 

    $ pairfq

    ERROR: Command line not parsed correctly. Check input.

    USAGE: pairfq [-h] [-m] [--version]

    Required:
        addinfo           :      Add the pair info back to the FastA/Q header.
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

**TASKS**

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

* **makepairs**

You have quality/adapter trimmed two paired-end sequence files and now they are out of sync. In this case, it is necessary to re-pair them, and then interleave the pairs for assembly.

    $ pairfq makepairs -f s_1_1_trimmed.fq -r s_1_2_trimmed.fq -fp s_1_1_trimmed_p.fq -rp s_1_2_trimmed_p.fq -fs s_1_1_trimmed_s.fq -rs s_1_2_trimmed_s.fq --index

In the above command, we specify the `makepairs` positional argument for pairing reads. The short arguements are `-f` for the file of forward reads, `-r` for the reverse reads, `-fp` for the file of paired forward reads, `-rp` for the file of reverse paired reads, `-fs` for the file of forward singleton/unpaired reads, and `-rs` for the singleton/unpaired reverse reads. 

The last argument, `--index`, is optional and specifies that an index will be constructed and all computation will be done on disk. This is a good thing if you are trying to pair two files of around 50 million reads on a machine with, for example, 8 GB of RAM. The computation will be much slower but much less memory will be used. If you have a moderate amount of memory and not so many reads, omit this last option, as the processing will go much faster. Below are some rough benchmarks (with and without the `--index`) for `pairfq makepairs` using a file of 10.7 million forward reads and 10.8 million reverse reads.

    Command                                         Time (utime)    RAM (RSS)
    pairfq makepairs ...                            17min04s        4.98G
    pairfq makepairs ... --index w/BerkeleyDB       4hr20min        4.62G
    pairfq makepairs ... --index 

Note that there are two different benchmarks with the `--index` command. The first one uses the BerkeleyDB package, which is much faster. You can install it with:

    $ cpanm BerkeleyDB

Installing BerkeleyDB isn't required, but it is clear that the indexing method will be much slower if this package is not available. These figures should be taken with caution, as they will vary depending on the machine and obviously, the amount of data being processed.
 
* **joinpairs**

With this command we can interleave the files for assembly or mapping.

    $ pairfq joinpairs -f s_1_1_trimmed_p.fq -r s_1_2_trimmed_p.fq -o s_1_interl.fq

In the above command, we are doing all computation without an index for speed. Now we can use our interleaved pairs, along with the unpaired reads for added coverage. 

* **addinfo**

Some toolkits discard the FASTQ comment which results in losing the pair information (e.g., `seqtk sample`). Therefore, it is necessary to add this information back before pairing or assembly. 

    $ pairfq addinfo -i s_1_1_sample_500k.fq -o s_1_1_sample_500k_pair.fq -p 1

This is a rather simple command, we specify the forward file that has been modified with `-i` and the corrected file with `-o`. We want to add the forward pair information with `-p 1` and if this were the reverse pair we would simply say `-p 2`.

* **splitpairs**

Sometimes it is necessay to split your interleaved file of forward and reverse reads into separate files.

    $ pairfq splitpairs -i s_1_interl.fq -f s_1_1_p.fq -r s_1_2_p.fq

The conventions used here are the same as with all the commands, `-f` specifies the file of forward reads to create and `-r` the reverse reads.  

**ISSUES**

Report any issues at: https://github.com/sestaton/Pairfq/issues

Be aware that Pairfq will not work for every data set given the wide range of FASTA/Q formats. Feel free to fork the project and modify the code to your needs, or submit code changes if you find any bugs. 

**ATTRIBUTION**

This project uses the [readfq](https://github.com/lh3/readfq) library written by Heng Li. The readfq code has been modified for error handling and to parse the comment line in the Casava header.

**LICENSE**

The MIT License should included with the project. If not, it can be found at: http://opensource.org/licenses/mit-license.php

Copyright (C) 2013 S. Evan Staton

[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/sestaton/pairfq/trend.png)](https://bitdeli.com/free "Bitdeli Badge")