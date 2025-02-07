mod paca;
mod parser;

use std::io::{self};
use std::time::Instant;

use crate::paca::{main_index_construction, query_bloom_filter, read_fof_file, merge_multiple_indexes};
use crate::parser::parse_args;

fn main() -> io::Result<()> {
    let matches = parse_args();

    let mode = matches
        .get_one::<String>("mode")
        .expect("Required argument for mode (index or query)");

    match mode.as_str() {
        "index" => {
            let fof_file = matches
                .get_one::<String>("fof")
                .expect("Required argument for FOF file in index mode");
            
            let kmer = matches
                .get_one::<String>("kmer")
                .expect("Required argument for k-mer size")
                .parse::<usize>()
                .expect("Invalid k-mer size");

            let bloomfilter = matches
                .get_one::<String>("bloomfilter")
                .map(|s| s.parse::<usize>().expect("Invalid Bloom filter size"))
                .unwrap_or(26); // default size
            
            let minimizer = matches
                .get_one::<String>("minimizer")
                .map(|s| s.parse::<usize>().expect("Invalid minimizer size"))
                .unwrap_or(11); // default size

            let bf_size = 1u64 << bloomfilter; // Bloom filter size as a power of 2
            let partitions = matches
                .get_one::<String>("partitions")
                .map(|s| s.parse::<usize>().expect("Invalid number of partitions"))
                .unwrap_or(8192); // default number of partitions
            
            let abundance = matches
                .get_one::<String>("abundance")
                .map(|s| s.parse::<usize>().expect("Invalid abundance number"))
                .unwrap_or(256); // default abundance levels

            let abundance_max = matches
                .get_one::<String>("abundance_max")
                .map(|s| s.parse::<u16>().expect("Invalid maximal abundance"))
                .unwrap_or(65535);
            
            let output_dir = matches.get_one::<String>("output");

            
            let start_time = Instant::now();

            // read the file of files  and extract file paths and color count
            let (file_paths, color_nb) = read_fof_file(fof_file)?;

            // run the index construction process: build and fill BFs per partitions and in chunks, serialize, merge chunks
            main_index_construction(
                file_paths,
                kmer,
                minimizer,
                bf_size,
                partitions,
                color_nb,
                abundance,
                abundance_max,
                output_dir.as_deref().map(|x| x.as_str()), //pass optional output dir, todo in passed arguments
            )?;

            println!("Indexing complete in {:.2?}", start_time.elapsed());
        }

        "query" => {
            let fasta_file = matches
                .get_one::<String>("fasta")
                .expect("Required argument for FASTA file in query mode");
            
            let index_dir = matches
                .get_one::<String>("index")
                .expect("Required argument for index directory in query mode");

            let color_graph = matches
                .get_one::<String>("color")
                .map(|s| s.parse::<bool>().expect("Invalid color option"))
                .unwrap_or(false);


            println!("Index directory: {}", index_dir);


            let mut query_output = format!("{}/query_results.csv", index_dir);
            if color_graph {
                query_output = format!("{}/colored_graph.fa", index_dir);
            }

            let start_time = Instant::now();
            query_bloom_filter(fasta_file, index_dir, &query_output, color_graph).expect("Failed to query sequences");

            println!("Query complete in {:.2?}", start_time.elapsed());
        }
        "merge" => {
            // argument= path to a fof + output file
            let indexes_fof = matches
                .get_one::<String>("indexes")
                .expect("Required argument: indexes (file-of-index directories)");
            let output_dir = matches.get_one::<String>("output");
            merge_multiple_indexes(indexes_fof, output_dir.as_deref().map(|x| x.as_str()))?;
        }
        _ => {
            eprintln!("Invalid mode: {}. Use 'index', 'query', or 'merge'.", mode);
        }
    }

    Ok(())
}
