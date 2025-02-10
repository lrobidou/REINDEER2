# REINDEER 2

## Installation

### Compilation

The installations steps are as follows:

* Clone the repository :

    `git clone --recursive https://github.com/Yohan-HernandezCourbevoie/REINDEER2.git `

* Then build :

    `cd REINDEER2`
    
    `cargo build`


## Usage

The generale use of REINDEER 2 is divided in to steps : index building and abundance query.

### Index

For the **index** mode, the mandatory parameters are the file of files (a plain text file where each line represented a unitigs file) and the size of the k-mers to be indexed.

`PACA --mode index -fof file_of_files.txt --kmer 31`


General parameters:
- **(-o, --output)** -- an output dir for the index
- **(-a, --abundance)** -- the abundance granularity (number of levels or discretized abundance values)
- **(-A, --abundance-max)** -- the maximal abundance to take into account

Advanced parameters: 
- **(-b, --bloomfilter)** -- the Bloom filter size in log2 scale
- **(-m, --minimizer)** -- the minimizer size
- **(-p, --partitions)** -- the number of partitions



### Query

For **query** mode, the parameters are the FASTA file containing the sequence(s) to be queried and the index directory.

`PACA --mode query --fasta seqeunces_query.fa --index ~/index_directory`

The **(-c, --color) parameter** <true/false> is used to define the query output format.

#### CSV file with --color false (default)

This option outputs a CSV file (header included) with the following structure : \<sequence header\>,\<abundance\>,\<indexed file number\>

#### FASTA file with --color false 

This option outputs the FASTA file given in query annotated with the abundances of all indexed files.


