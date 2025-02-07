# Benchs on Metagraph
# 07/01/2025

import os
import sys


if __name__ == "__main__":

    OVERWRITE = False
    
    fofs = ["./fof1_fastq.lst"]  # list of file_of_files

    time_command = "/usr/bin/time -f %U,%S,%e,%M,%C -a -o times.csv"

    # PARAMS
    metagraph_conda = "metagraph"
    kmc = "kmc"
    min_count = 2
    k = 31
    threads = 32
    tmp_dir = "./tmp"
    max_ram = 25
    anno_type = "int_brwt"

    for fof in fofs:
        files = [f.strip() for f in open(fof, "r").readlines() if f != "\n"]
        assert(len(files)%2 == 0) # fastq paired

        fof_basename = fof.split("/")[-1].split('.')[0]
        accessions = []
        os.system(f"{time_command} echo 'Working now on {fof_basename}'")

        # build the general graph => {fof_basename}.graph.dbg
        os.system(f"{time_command} {metagraph_conda} build --min-count {min_count} --mode canonical -k {k} -p {threads} --disk-swap {tmp_dir} --mem-cap-gb {max_ram} -o {fof_basename}.graph {" ".join(files)}")

        # for each pair of fastq
        for i in range(0, len(files), 2):
            file1, file2 = files[i:i+2]
            accession = os.path.basename(file1).split("/")[-1].split("_")[0]
            accessions.append(accession)

            subgraphs_done = False
            kmc_done = False
            
            if OVERWRITE or not os.path.exists(f"metacounts/{accession}.fasta.gz"):
                
                if OVERWRITE or not (os.path.exists(f"{accession}.dbg") and os.path.exists(f"{accession}.dbg.weights")):
                    
                    if OVERWRITE or not os.path.exists(f"{accession}.kmc_suf"):
                        
                        # pre-count the kmers => kmc_files/{file_basename}.kmc_pre + kmc_files/{file_basename}.kmc_suf
                        open(f"fof_{accession}.txt", "w").write(f"{file1}\n{file2}\n")
                        os.system(f"{time_command} {kmc} -k{k} -m{max_ram} -ci{min_count} -fq -t{threads} @fof_{accession}.txt kmc_files/{accession} {tmp_dir}/")
                        os.system(f"rm fof_{accession}.txt")
                        kmc_done = True
                    
                    # build a small graph for just one sample (with paired fastq) => subgraphs/{file_basename}.dbg + subgraphs/{file_basename}.dbg.weights
                    os.system(f"{time_command} {metagraph_conda} build --mode canonical -k {k} --count-kmers --count-width 8 --disk-swap {tmp_dir} --mem-cap-gb {max_ram} --min-count {min_count} -o subgraphs/{accession} kmc_files/{accession}.kmc_suf")
                    subgraphs_done = True

                # clean the graph and output the fasta with abundances => metacounts/{accession}.fasta.gz
                os.system(f"{time_command} {metagraph_conda} clean --to-fasta --smoothing-window 1 -p {threads} --min-count {min_count} -o metacounts/{accession} subgraphs/{accession}.dbg")
            
            # delete "temporary" files
            if kmc_done:
                os.system(f"rm kmc_files/{accession}.kmc_pre kmc_files/{accession}.kmc_suf")
            if subgraphs_done:
                os.system(f"rm subgraphs/{accession}.dbg subgraphs/{accession}.dbg.weights")

        # annotate the general graph with each color => {fof_basename}.column.annodbg + {fof_basename}.column.annodbg.counts
        os.system(f"{time_command} {metagraph_conda} annotate -p {threads} -i {fof_basename}.graph.dbg --anno-filename --separately --count-kmers --count-width 8 -o {fof_basename} {' '.join('metacounts/'+acc+'.fasta.gz' for acc in accessions)}")

        # transform annotations
        os.system(f"{time_command} {metagraph_conda} transform_anno --anno-type int_brwt --count-kmers -p {threads} --greedy --disk-swap {tmp_dir} -o {fof_basename}_anno {fof_basename}/*.column.annodbg")

