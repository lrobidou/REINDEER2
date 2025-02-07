# Benchs on Reindeer
# 07/01/2025

import os
import sys


if __name__ == "__main__":
    reindeer_release = "/Reindeer" # path to Reindeer
    time_command = "/usr/bin/time -f %U,%S,%e,%M,%C -a -o times.csv "
    filelist = ["fof1.txt", "fof2.txt", "fof4.txt", "fof8.txt", "fof16.txt", "fof32.txt", "fof64.txt", 
                "fof128.txt", "fof256.txt", "fof512.txt", "fof.txt"] # file_of_files

    for fof in filelist:
        command = f"{reindeer_release} --index -f {fof} -k 31 -t 32 --log-count -o ./{fof.split('.')[0]}_logcount/"
        os.system(time_command+command)
