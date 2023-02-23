from typing import Callable

import pandas as pd
import re
import argparse


def get_sequences_from_file(filepath: str) -> pd.DataFrame:
    # is it .excel? or .csv?
    match filepath.split(".")[-1]:
        case "csv":
            reader = pd.read_csv
        case "xlsx":
            reader = pd.read_excel
        case x:
            raise ValueError(f"Unknown filetype! Can't parse .{x}")
    return reader(filepath)


def get_codons(seq: str, n=3) -> list[str]:
    """
    A function that takes an input DNA sequence representing an open
    reading frame (ORF)and splits that sequence into codons of length n
    (default=3)

    Parameters
    ----------
        str seq: string representing an ORF (DNA sequence)
        int n: length of codons (default=3)

    Returns
    -------
        list<str> codons: input sequence split into codons
    """
    # check that sequence is divisible by n
    if len(seq) % n != 0:
        raise ValueError(f"seq is not divisible by n ({n})")

    num_codons = int(len(seq) / n)
    codons = [
        seq[n * i:n * (i + 1)] for i in range(num_codons)
    ]

    return codons


def translate(seq: str) -> str:
    """
    A function that takes a DNA sequence as an input and returns its translation

    Parameters
    ----------
        str seq: string representing an ORF (DNA sequence)

    Returns
    -------
        str protein_sequence: amino acid sequence corresponding to translation of seq
    """
    standard_code = {
        'TTT': 'F', 'TTC': 'F',
        'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y',
        'TAA': '*', 'TAG': '*', 'TGA': '*',
        'TGT': 'C', 'TGC': 'C',
        'TGG': 'W',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H',
        'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
        'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N',
        'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S',
        'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D',
        'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    codons = get_codons(seq)
    aminos = [standard_code[c] for c in codons]
    protein_sequence = "".join(aminos)
    return protein_sequence


def parse_spec_element(elem: str) -> str:
    match list(elem):
        case ["[", *x, "]"]:
            options = "|".join(x)
            new_elem = f"[{options}]{{1}}"
        case aminos if any(aa.casefold() == "x".casefold() for aa in aminos):
            new_elem = f".{{{len(aminos)}}}"
        case x:
            new_elem = "".join(x)

    return new_elem


def parse_spec(spec: str, capture: bool = False) -> str:
    spec_element_list = spec.split("-")

    parsed_elem_list = [
        parse_spec_element(elem)
        for elem in spec_element_list
    ]
    if capture:
        parsed_elem_list = [
            f"({elem})"
            for elem in parsed_elem_list
        ]
    full_spec = "".join(parsed_elem_list)
    return full_spec

def give_prefixes(prefix: str = None, suffix: str =None) -> tuple[str, str]:
    """
    Start: aGGTCTCaA
    End:   GCTTtGAGACCt
    """
    if prefix is None:
        prefix = "aGGTCTCaA"
    if suffix is None:
        suffix = "GCTTtGAGACCt"

    return prefix, suffix


def add_restriction_sites(seq: str) -> str:
    prefix, suffix = give_prefixes()
    new_seq = prefix + seq + suffix
    return new_seq


def cut_codons(seq: list[str], indices: tuple[int, int], codon: str) -> list[str]:
    low, high = indices
    new_seq_list = seq[:low] + seq[high:]
    return new_seq_list


def replace_codons(seq: list[str], indices: tuple[int, int], codon: str) -> list[str]:
    low, high = indices
    replacement_codon = [codon for _ in range(high-low)]
    new_seq_list = seq[:low] + replacement_codon + seq[high:]
    return new_seq_list


def sequence_transform_factory(
        func: Callable[[list[str], tuple[int, int], str], list[str]],
        degron_spec: str, codon: str,
        add_prefix: bool = True, has_prefix: bool = True,
) -> Callable[[str], str]:
    def sequence_transform(seq: str) -> str:
        # optionally strip prefixes first
        if has_prefix:
            add_prefix = True
            prefix, suffix = give_prefixes()
            seq = seq[len(prefix):-len(suffix)]

        # translate to protein sequence
        prot_seq = translate(seq)

        # use regex to find instances of degron
        match = re.search(degron_spec, prot_seq)

        if match is None:
            new_seq = seq
        else:
            # split sequence into codons
            old_codon_list = get_codons(seq)
            indices = match.span()
            updated_codon_list = func(old_codon_list, indices, codon)
            new_seq = "".join(updated_codon_list)

        if add_prefix:
            new_seq = add_restriction_sites(new_seq)

        return new_seq
    return sequence_transform


def new_dataframe(df: pd.DataFrame, colname: str, func: Callable[[str], str]) -> pd.DataFrame:
    new_df = df.copy(deep=True)
    new_df[colname] = new_df[colname].map(func)
    return new_df


def new_filenames(old_filename: str, modifier: str, outfile_extension: str = "csv") -> str:
    parent, extension = old_filename.split(".")
    new_filename = f"{parent}_{modifier}.{outfile_extension}"
    return new_filename


if __name__ == "__main__":
    # handle user input
    parser = argparse.ArgumentParser(description='Remove degrons and replace with poly-Alanines')
    parser.add_argument('filepath', metavar='filepath', type=str,
                        help='filepath to DNA sequences (.xlsx or .csv)')
    parser.add_argument('--has_prefix', type=bool, default=False, help="Do the DNA sequences already have prefixes?")
    parser.add_argument('--colname', type=str, default="Sequence",
                        help='column name with DNA seqeunces (default: Bases', required=False)
    parser.add_argument('--degron', type=str, default="[VI]-GWPP-[VIHSK]-[GR]-xx-R",
                        help='degron sequence to remove (default: [VI]-GWPP-[VIHSK]-[GR]-xx-R)', required=False)
    parser.add_argument('--codon', type=str, default="GCA",
                        help='codon to replace degron sequence with (default: GCA)', required=False)
    args = parser.parse_args()


    # parse and process args
    filepath = args.filepath
    old_sequences = get_sequences_from_file(filepath)
    has_prefixes = args.has_prefix
    degron_spec = parse_spec(args.degron)
    codon = args.codon
    colname = args.colname

    # get transforms
    functions_to_use = {
        "deleted_degron": cut_codons,
        "poly_alanine": replace_codons
    }
    transforms = {
        name: sequence_transform_factory(func, degron_spec, codon, add_prefix=True)
        for name, func in functions_to_use.items()
    }
    new_dfs: dict[str, pd.DataFrame] = {
        name: new_dataframe(old_sequences, colname, func)
        for name, func in transforms.items()
    }
    # save new dataframes
    for name, df in new_dfs.items():
        outpath = new_filenames(filepath, name)
        df.to_csv(outpath)