from Bio import AlignIO


def has_distinct_seqs(msa):
    formats = ["fasta", "phylip", "phylip_relaxed", "iphylip_relaxed", "phylip_interleaved"]
    for f in formats:
        try:
            align = AlignIO.read(msa, f)
        except:
            pass
    return len(set(map(lambda r: r.seq, align))) > 3
