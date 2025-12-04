import argparse
import numpy
import random
import RNA
import csv 
import os
import subprocess

def sequence_shuffle(sequence, n):
    sequence = sequence.upper()
    sequences = [sequence]
    for _ in range(n):
        fasta_input = f">seq\n{sequence}\n"
        seed = random.randint(1, 1_000_000)
        proc = subprocess.Popen(
            ["fasta-shuffle-letters", '-kmer', '2', '-rna', '-', '-seed', str(seed)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout,stderr = proc.communicate(input=fasta_input)

        if proc.returncode !=0:
            raise RuntimeError(f'fasta shuffles failes: {stderr}')

        shuffled_seq = ''.join([line.strip() for line in stdout.splitlines() if not line.startswith('>')])
        sequences.append(shuffled_seq)
    
    return sequences

def cal_z(sequence, n):
    sequences = sequence_shuffle(sequence, n)

    mfes = []
    for seq in sequences:
        fc = RNA.fold_compound(seq)
        structure, mfe_value = fc.mfe()
        mfes.append(mfe_value)
    
    original_mfe = mfes[0]
    shuffled_mfes = mfes[1:]       

    mean_shuf = numpy.mean(shuffled_mfes)
    std_shuf = numpy.std(shuffled_mfes)

    z_scores = [
        (mfe - mean_shuf) / std_shuf if std_shuf != 0 else 0
        for mfe in mfes
    ]

    return sequences, mfes, z_scores, original_mfe

def main():

    parser = argparse.ArgumentParser(
        description = "Shuffle RNA sequence, calculate MFEs and z-score."
    )

    parser.add_argument("sequence", type=str, help="The RNA sequence.")

    parser.add_argument(
        "-n", "--number",
        type=int,
        default=100,
        help = "Number of shuffles performed.(Default:100)"
    )

    parser.add_argument(
        "-o", "--output",
        type=str,
        help = "TSV output file containing the sequences, MFEs and the z-score."
    )

    parser.add_argument(
        "-f", "--folder",
        type=str,
        help = "Makes a folder containing all the output files within the directory."
    )
    
    args = parser.parse_args()

    sequences, mfes, z_scores, original_mfe = cal_z(
        args.sequence,
        args.number
    )

    if args.output:
        if args.folder:
            os.makedirs(args.folder, exist_ok=True)
            output_path = os.path.join(args.folder, args.output)
        else:
            output_path = args.output

        with open(output_path, 'w', newline='') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            writer.writerow(["seq", "mfe", "z_score"])
            writer.writerow([args.sequence.upper(), original_mfe, z_scores[0]])

            for seq, mfe, z_score in zip(sequences[1:], mfes[1:], z_scores[1:]):
                writer.writerow([seq, mfe, z_score])
    else:
        print("Sequence:", args.sequence)
        print("Original MFE:", original_mfe)
        print("z-score:", z_scores[0])

if __name__ == "__main__":
    main()
