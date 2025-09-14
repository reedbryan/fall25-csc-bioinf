from dbg import DBG
from utils import read_data
import sys
import os
import time
from n50 import parse_fasta_lengths, compute_n50  # <-- import functions

sys.setrecursionlimit(1000000)


if __name__ == "__main__":
    total_start_time = time.time()
    argv = sys.argv
    short1, short2, long1 = read_data(os.path.join('../data/', argv[1]))

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])
    # dbg.show_count_distribution()
    contig_path = os.path.join('../data/', argv[1], 'contig.fasta')
    with open(contig_path, 'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            print(i, len(c))
            f.write('>contig_%d\n' % i)
            f.write(c + '\n')

    # Compute and print N50
    lengths = parse_fasta_lengths(contig_path)
    n50 = compute_n50(lengths)
    print(f'N50: {n50}')
    
    # Print total execution time
    total_end_time = time.time()
    total_time = total_end_time - total_start_time
    print(f"Total execution time: {total_time:.4f} s")


