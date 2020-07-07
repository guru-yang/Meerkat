Meerkat is a software package designed to detect structural variations (SVs) from sequencing data. The detailed method description can be found at:
https://www.cell.com/fulltext/S0092-8674(13)00451-0

The current version is 0.189. Detailed instructions of installation and usage can be found in Manual_0.189.pdf.

Frequently Asked Questions:

0. I have a problem, is there documentation available?
Yes! A manual is distributed with Meerkat. Before contacting us for help, please be sure to read the manual first.

1. I can't compile the binaries.
If you encounter error like this:
/usr/bin/ld: gzstream.o: undefined reference to symbol 'gzclose'
//lib/x86_64-linux-gnu/libz.so.1: error adding symbols: DSO missing from command line
collect2: error: ld returned 1 exit status
make: *** [bamreader] Error 1
Please modify the Makefile of each binary as follow:
from: ... -lbamtools -lbamtools-utils
to: ... -lbamtools -lbamtools-utils -lz

2. What if I can't get Meerkat to work on my server?
Please use the latest version of Meerkat. Meerkat requires a number of programs, modules, and references, make sure you have all of them in place and the parameters are given properly. The most likely error is some parameters are not specified correctly. If you have difficulty installing any programs, modules or references, please first contact your system admin. If you need further assistance, please contact us at ylixing@gmail.com	
When you contact us, please provide the following: Are you able to run ./bin/bamreader from command?
Are you able to run example.bam?
The output of "ls -l" for run folder.
pre.log file if it's generated.
isinfo file if it's generated.
dre.log file if it's generated.
Error message if there is any.

3. The breakpoints are not properly annotated with fusions.pl script
The refGene.txt downloaded from UCSC needs to be sorted by chromosome and coordinate by following command:
sort refGene.txt -k 3,3 -k 5,5n > refGene_sorted.txt
