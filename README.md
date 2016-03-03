# MismappingAnalysis

This project is intended to help researchers identify single nucleotide variants that could have been falsely called as variants due to mapping errors. It utilizes BLAT to identify homologous regions of the genome where the called variant is the reference. 

## Running the Program

### Start BLAT Server

To run you first need to start a blat server (gfServer). Here is an example of how to do this on Stanford's scg3 server:
/srv/gs1/software/blat/3.5/bin/x86_64/gfServer start localhost 8888 /srv/gs1/projects/snyder/collinmelton/bundle/2.3/b37/d5/hs37d5.2bit&

### Prepare an Input File

The default input file format is as follows (in tsv format): 

| chrom | pos | ref | var |
| --- | --- | --- | --- |
| 1 | 18866399 | A | G |
| 1 | 18866407 | G | A |
| 1 | 18866403 | G | A |

### Run CheckForMappingErrors.py

Run the program as follows:
python CheckForMappingErrors.py --I input file --O outputFile --L size of region to left of mut --R size of region to right of mut --P blat server port number --M fileType (T for merged file, F for not, default is F)

Since Illumina reads are usually 100ish in length I usually run with --L and --R = 50
