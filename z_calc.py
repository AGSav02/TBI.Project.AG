import argparse
import random
import numpy
import RNA
import csv 
import os

def sequence_shuffle(sequence, n):
    sequence = sequence.upper()
    sequences = [sequence]
    for _ in range(n):
        shuffled = list(sequence)
        random.shuffle(shuffled)
        shuffled = "".join(shuffled)
        sequences.append(shuffled)
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
    z_score = (original_mfe - mean_shuf) / std_shuf if std_shuf != 0 else 0

    return original_mfe, sequences[1:], shuffled_mfes, z_score

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

    original_mfe, sequences, shuffled_mfes, z_score = cal_z(
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
            writer.writerow([args.sequence.upper(), original_mfe, z_score])
            for seq, mfe in zip(sequences, shuffled_mfes):
                writer.writerow([seq, mfe, ''])
        print(f"File saved to: {output_path}")
    else:
        print("Sequence:", args.sequence)
        print("Original MFE:", original_mfe)
        print("z-score:", z_score)

if __name__ == "__main__":
    main()
