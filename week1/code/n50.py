def parse_fasta_lengths(fasta_path):
    lengths = []
    with open(fasta_path) as f:
        seq = ''
        for line in f:
            if line.startswith('>'):
                if seq:
                    lengths.append(len(seq))
                    seq = ''
            else:
                seq += line.strip()
        if seq:
            lengths.append(len(seq))
    return lengths

def compute_n50(lengths):
    lengths.sort(reverse=True)
    total = sum(lengths)
    half = total / 2
    running = 0
    for l in lengths:
        running += l
        if running >= half:
            return l