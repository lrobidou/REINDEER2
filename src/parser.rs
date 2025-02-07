use clap::{Arg, Command, ArgMatches};

pub fn parse_args() -> ArgMatches {
    Command::new("PACA")
        .version("0.1")
        .author("Camille Marchet")
        .about("k-mer quantification in large collections of samples")
        .arg(
            Arg::new("mode")
                .short('M')
                .long("mode")
                .value_name("MODE")
                .help("Mode: 'index' to build the index or 'query' to query a built index")
                .required(true),
        )
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .value_name("FILE")
                .help("For 'query' mode: Path to the FASTA file containing query sequences"),
        )
        .arg(
            Arg::new("index")
                .short('i')
                .long("index")
                .value_name("DIR")
                .help("For 'query' mode: Path to the directory containing the prebuilt index"),
        )
        .arg(
            Arg::new("fof")
                .short('F')
                .long("fof")
                .value_name("FILE")
                .help("For 'index' mode: Path to the input list of FASTA files (Logan format)"),
        )
        .arg(
            Arg::new("kmer")
                .short('k')
                .long("kmer")
                .value_name("SIZE")
                .help("For 'index' mode: Sets the k-mer size"),
        )
        .arg(
            Arg::new("bloomfilter")
                .short('b')
                .long("bloomfilter")
                .value_name("BF")
                .help("For 'index' mode: Sets the Bloom filter size (in log scale)"),
        )
        .arg(
            Arg::new("minimizer")
                .short('m')
                .long("minimizer")
                .value_name("MINSIZE")
                .help("For 'index' mode: Sets the minimizer size"),
        )
        .arg(
            Arg::new("partitions")
                .short('p')
                .long("partitions")
                .value_name("MINSIZE")
                .help("For 'index' mode: Sets the number of parititions"),
        )
        .arg(
            Arg::new("abundance")
                .short('a')
                .long("abundance")
                .value_name("ABUND")
                .help("For 'index' mode: Sets the abundance granularity"),
        )
        .arg(
            Arg::new("abundance_max")
                .short('A')
                .long("abundance-max")
                .value_name("ABUND_MAX")
                .help("For 'index' mode: Sets the maximal abundance to take into account"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("OUT")
                .help("For 'index' mode: Sets an output dir for the index"),
        )
        .arg(
            Arg::new("color")
                .short('c')
                .long("color")
                .value_name("COLOR")
                .help("For 'query' mode: bool for coloring a graph instead of regular query (default: false)"),
        )
        .get_matches()
}
