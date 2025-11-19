import argparse
import random

def sequence_shuffle(sequence, n):
    sequence = sequence.upper()
    shuffled_sequences = []
    for j in range(n):
        shuffled = list(sequence)
        random.shuffle(shuffled)
        shuffled_seq = ''.join(shuffled)
    return shuffled_sequences

def cal_z(sequence, n):
    mfes = []
    shuffled_sequences = sequence_shuffle(sequence, n)
    for seq in shuffled_sequences:
        fc = RNA.fold_compound(seq)
        structure, mfe_value = fc.mfe()
        mfes.append(mfe_value)

    original_mfe = mfes[0]
    shuffled_mfes = mfes
    mean_shuf = numpy.mean(shuffled_mfes)
    std_shuf = numpy.std(shuffled_mfes)
    z_score = (original_mfe - mean_shuf) / std_shuf if std_shuf != 0 else 0

def main():

    parser = argparse.ArgumentParser(
        description = "Shuffle 200nt sequence. Calculates the mfe of the original sequence and the its z-score."
    )

    parser.add_argument("sequence", type=str, help="The 200nt RNA sequence.")

    parser.add_argument(
        "-n", "--number",
        type=int,
        default=100,
        help = "Number of shuffles performed.(Default:100)"
    )

    parser.add_argument(
        "-o", "--output",
        type=str,
        help = "Output file to save the original sequence, the shuffled sequences along with the mfes and the calculated z-score."
    )

    args = parser.parse_args()

    shuffle_results = sequence_shuffle(args.sequence, args.number)
    z_results = cal_z(args.sequence)

    if args.output:
        with open(args.output, 'w', newline ='') as tsvfile:
            writer = csv.writer(tsvfile, delimiter = '\t')
            writer.writerow(["seq", 'mfe', 'z_score'])
            writer.writerow([shuffled_sequences[0][1], original_mfe, z_score])
            for seq, mfe in zip(sequences[1:], mfes[1:]):
                writer.writerow([seq, mfe, ''])
        print(f"Protein saved to {args.output}!")    
    else:
        print("mfe/z-score:")
        print(z_results)

if __name__ == "__main__":
    main()
