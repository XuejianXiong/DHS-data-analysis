import pandas as pd


# Functions used to create sequence columns
def sequence_bounds(summit: int, start: int, end: int, length: int):
    """
    Calculate the sequence coordinates (bounds) for a given DHS.

    https://github.com/meuleman/SynthSeqs/blob/main/make_data/process.py
    """
    half = length // 2

    if (summit - start) < half:
        return start, start + length
    elif (end - summit) < half:
        return end - length, end

    return summit - half, summit + half


def add_sequence_column(df: pd.DataFrame, genome, length: int):
    """
    Query the reference genome for each DHS and add the raw sequences
    to the dataframe.
    Parameters
    ----------
    df : pd.DataFrame
        The dataframe of DHS annotations and NMF loadings.
    genome : ReferenceGenome(DataSource)
        A reference genome object to query for sequences.
    length : int
        Length of a DHS.

    https://github.com/meuleman/SynthSeqs/blob/main/make_data/process.py
    """
    seqs = []
    for rowi, row in df.iterrows():
        l, r = sequence_bounds(row['summit'], row['start'], row['end'], length)
        seq = genome.sequence(row['seqname'], l, r)

        seqs.append(seq)

    df['sequence'] = seqs
    return df
