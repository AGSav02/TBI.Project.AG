import argparse
import pandas as pd
import os


def gc_content(seq):
    seq = seq.upper()
    if len(seq) == 0:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def add_gc_to_file(input_path, output_path):
    df = pd.read_csv(input_path, sep="\t")
    df["gc"] = df["seq"].apply(gc_content)
    df.to_csv(output_path, sep="\t", index=False)


def main():

    parser = argparse.ArgumentParser(
        description="Add GC content as a 4th column to an existing TSV file."
    )

    parser.add_argument(
        "input",
        type=str,
        help="Input TSV file (expects columns: seq, mfe, zscore)."
    )

    parser.add_argument(
        "-o", "--output",
        type=str,
        help="Name of output TSV file."
    )

    parser.add_argument(
        "-f", "--folder",
        type=str,
        help="Output folder to store annotated files."
    )

    args = parser.parse_args()

    if args.folder:
        os.makedirs(args.folder, exist_ok=True)
        if args.output:
            output_path = os.path.join(args.folder, args.output)
        else:
            output_filename = os.path.basename(args.input)
            output_path = os.path.join(args.folder, output_filename)
    else:
        if args.output:
            output_path = args.output
        else:
            output_path = args.input.replace(".tsv", "_gc.tsv")

    add_gc_to_file(args.input, output_path)



if __name__ == "__main__":
    main()
