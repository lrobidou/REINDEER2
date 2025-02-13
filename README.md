# REINDEER 2

## Installation

### Compilation

Clone the repository :

```
git clone https://github.com/Yohan-HernandezCourbevoie/REINDEER2.git 
```

Then build :

```
cd REINDEER2 && cargo build
```

Alternatively, the tool can be installed globally in the system using :

```
cd REINDEER2 && cargo install --path .
```

In the following examples, the tool's command will use its installation name `Reindeer2` but it can be replaced by `cargo run --` when used from the REINDEER2 folder.


## Usage

The generale use of REINDEER 2 is divided in to steps : index building and abundance query.

### Index

For the **index** mode, the mandatory parameters are the file of files (a plain text file where each line represented a unitigs file) and the size of the k-mers to be indexed.

`Reindeer2 --mode index -fof file_of_files.txt --kmer 31`


General parameters:
- `-o, --output` an output directory for the index
- `-a, --abundance` the abundance granularity (number of levels or discretized abundance values)
- `-A, --abundance-max` the maximal abundance to take into account

Advanced parameters: 
- `-b, --bloomfilter` the Bloom filter size in log2 scale
- `-m, --minimizer` the minimizer size
- `-p, --partitions` the number of partitions



### Query

For **query** mode, the parameters are the FASTA file containing the sequence(s) to be queried and the index directory.

`Reindeer2 --mode query --fasta seqeunces_query.fa --index ~/index_directory`

The `-c, --color` parameter (*true*/*false*) is used to define the query output format.

#### CSV file with --color false (default)

This option outputs a CSV file (header included) with the following structure : `<Sequence_header>,<Color>,<Median_abundance>`

#### FASTA file with --color true

This option outputs the FASTA file given in query annotated with the abundances of all indexed files.

## Example

To illustrate how the tool works, a small example is available. The commands are launched from the REINDEER2 main directory.

#### INDEX
How to build the index:
```
Reindeer2 --mode index --fof test_files/fof.txt --kmer 31 -output ../index_test
```

#### QUERY (results: CSV)
With the command:
```
Reindeer2 --mode query --fasta test_files/file1Q.fa --index ../index_test
cat ../index_test/query_results.csv
```
is expected the result:
```
header,file,abundance
>seq1 ka:f:30,0,29
>seq2 ka:f:30,0,29
>seq3 ka:f:1500,0,1450
>seq3 ka:f:1500,1,4
```


#### QUERY (results: colored graph FASTA)
With the command:
```
Reindeer2 --mode query --fasta test_files/file1Q.fa --index ../index_test --color true 
cat ../index_test/colored_graph.fa
```
is expected the result:
```
>seq1 ka:f:30 col:0:29 col:1:0
AAAAAAAAAAAAAAAAAAAAAACACAGATCA
>seq2 ka:f:30 col:0:29 col:1:0
AAAAAAAAAAAAAAAAAAAAACACAGATCAT
>seq3 ka:f:1500 col:0:1450 col:1:4
AAAAAAAAAAAAAAAAAAAAAACAAAAAGAA
```


