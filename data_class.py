import numpy as np
import pandas as pd
from Bio import SeqIO


class DataSource:
    # Sourced from https://github.com/meuleman/SynthSeqs/blob/main/make_data/source.py

    def __init__(self, data, filepath):
        self.raw_data = data
        self.filepath = filepath

    @property
    def data(self):
        return self.raw_data


class ReferenceGenome(DataSource):
    """Object for quickly loading and querying the reference genome."""

    @classmethod
    def from_path(cls, path):
        genome_dict = {record.id: str(record.seq).upper() for record in SeqIO.parse(path, "fasta")}
        return cls(genome_dict, path)

    @classmethod
    def from_dict(cls, data_dict):
        return cls(data_dict, filepath=None)

    @property
    def genome(self):
        return self.data

    def sequence(self, chrom, start, end):
        chrom_sequence = self.genome[chrom]

        assert end < len(chrom_sequence), (
            f"Sequence position bound out of range for chromosome {chrom}. "
            f"{chrom} length {len(chrom_sequence)}, requested position {end}."
        )
        return chrom_sequence[start:end]
    

class FilteringData:
    def __init__(self, df: pd.DataFrame, cell_list: list):
        self.df = df
        self.cell_list = cell_list
        self._test_data_structure()

    def _test_data_structure(self):
        # Ensures all columns after the 11th are named cell names
        assert all("_ENCL" in x for x in self.df.columns[11:]), \
                        "_ENCL not in all columns after 11th"

    def filter_exclusive_replicates(self, sort: bool = False, 
                                    balance: bool = True):
        """
        Given a specific set of samples (one per cell type),
        capture the exclusive peaks of each samples (the ones matching 
        just one sample for the whole set)
        and then filter the dataset to keep only these peaks.

        Returns:
            pd.DataFrame: The original dataframe plus a column for each 
            cell type with the exclusive peaks
        """
        print("Filtering exclusive peaks between replicates")
        # Selecting the columns corresponding to the cell types
        subset_cols = self.df.columns[:11].tolist() + self.cell_list
        # Creating a new dataframe with only the columns corresponding to 
        # the cell types
        df_subset = self.df[subset_cols]
        # Creating a new column for each cell type with the exclusive peaks
        # or "NO_TAG" if not exclusive
        df_subset["TAG"] = df_subset[self.cell_list]\
                              .apply(lambda x: "NO_TAG" if x.sum() != 1 
                              else x.idxmax(), axis=1)

        # Creating a new dataframe with only the rows with exclusive peaks
        new_df_list = []
        for k, v in df_subset.groupby("TAG"):
            if k != "NO_TAG":
                cell, replicate = "_".join(k.split("_")[:-1]), k.split("_")[-1]

                v["additional_replicates_with_peak"] = (
                    self.df[self.df.filter(like=cell).columns]\
                        .apply(lambda x: x.sum(), axis=1) - 1)
                
                temp_df = self.df.filter(like=cell)
                print(f"Cell type: {cell}, Replicate: {replicate}, \
                      Number of exclusive peaks: {v.shape[0]}")
            else:
                v["additional_replicates_with_peak"] = 0

            new_df_list.append(v)

        new_df = pd.concat(new_df_list).sort_index()
        new_df["other_samples_with_peak_not_considering_reps"] = (
                  new_df["numsamples"] - \
                  new_df["additional_replicates_with_peak"] - 1)

        # Sorting the dataframe by the number of samples with the peak
        if sort:
            new_df = pd.concat(
                [x_v.sort_values(by=["additional_replicates_with_peak", 
                        "other_samples_with_peak_not_considering_reps"],
                        ascending=[False, True])
                        for x_k, x_v in new_df.groupby("TAG")],
                ignore_index=True)

        # Balancing the dataset
        if balance:
            lowest_peak_count = new_df.groupby("TAG").count()["sequence"].min()
            new_df = pd.concat([v_bal.head(lowest_peak_count) 
                        for k_bal, v_bal in new_df.groupby("TAG") 
                        if k_bal != "NO_TAG"])

        return new_df
