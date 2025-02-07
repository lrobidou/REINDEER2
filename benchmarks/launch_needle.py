# Benchs on Needle
# 07/01/2025

import os
import sys 



if __name__ == "__main__":

    # fof_arg => un fof ou plusieurs séparés par une simple virgule
    fof_arg = sys.argv[1]
    if ',' in fof_arg :
        fof_list = fof_arg.split(",")
    else:
        fof_list = [fof_arg]

    # params
    needle = "needle"
    k = 31
    output_dir = "./"
    time_command = f"/usr/bin/time -f %U,%S,%e,%M,%C -a -o times.csv "
    threads = 32
    levels = 16
    false_positive_rate = 0.3


    for fof in fof_list:
        files = " ".join(path.strip() for path in open(fof, "r").readlines() if path != "\n")

        # testing with 2 window sizes : 41 and 51
        for w in [41, 51]:
            needle_minimiser_dir = f"{output_dir}needle_minimiser_{os.path.basename(fof).split('.')[0]}_k{k}_w{w}/"
            if not os.path.exists(needle_minimiser_dir):
                os.system("mkdir" + needle_minimiser_dir)

            needle_ibfmin_dir = f"{output_dir}needle_ibfmin_{os.path.basename(fof).split('.')[0]}_k{k}_w{w}/"
            if not os.path.exists(needle_ibfmin_dir):
                os.system("mkdir" + needle_ibfmin_dir)
            
            # Pre-processing minimisers from paired fastq files
            minimiser_command = f"{needle} minimiser {files} --paired -k {k} -w {w} -t {threads} -o {needle_minimiser_dir}"
            os.system(time_command + minimiser_command)

            # Indexing the minimisers
            ibfmin_command = f"{needle} ibfmin {needle_minimiser_dir}* -l {levels} -f {false_positive_rate} -t {threads} -o {needle_ibfmin_dir}"
            os.system(time_command + ibfmin_command)


