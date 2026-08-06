from pathlib import Path

def parse_fasta(path):
    "Read a FASTA file and return the sequence"
    sequence = {}
    seq_id = None
    parts = []
    with open(path) as f:
        for line in f:              
            line = line.strip()
            if not line:
                continue           
            if line.startswith(">"):
                if seq_id is not None:
                    sequence[seq_id] = "".join(parts)
                seq_id = line[1:].split()[0]
                parts =[]
            else:
                parts.append(line.upper())
        if seq_id is not None:
            sequence[seq_id] = "".join(parts)
    return sequence

class Sequence:
    "A single DNA sequence with basic analysis methods"

    def __init__(self,seq_id,seq):
        self.seq_id = seq_id
        self.seq = seq.upper()

    def lenght(self):
        return len(self.seq)

    def gc_content(self, as_percent= True):
        "Return the GC fraction (or percent) of the sequnce"
        gc = self.seq.count("G") + self.seq.count("C")
        frac = gc / len(self.seq)
        return frac * 100 if as_percent else frac

    def __repr__(self):
        return f"Sequence({self.seq_id}, len={self.lenght()})"



if __name__ == "__main__":
    BASE_DIR = Path(__file__).resolve().parent
    fasta_file = BASE_DIR / "sequences.fasta"
    records = parse_fasta(fasta_file)
    for seq_id, seq in records.items():
        sequence = Sequence(seq_id,seq)
        print(sequence, "GC=", round(sequence.gc_content(),2))





            