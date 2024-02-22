import os
import numpy as np
import pandas as pd
from typing import *
from Bio import SeqIO
from IPython.display import display

from data_class import *
from create_seq_columns import *


def download_data() -> None:

    db_path = "./database"
      
    #### Downloading human genome
    os.system(f"wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P {db_path}/")
    os.system(f"gzip -d {db_path}/hg38.fa.gz")

    #### Downloading 733 biosample metadata
    # Metadata files describing biosample characteristics and annotations
    os.system(f"wget https://www.meuleman.org/DHS_Index_and_Vocabulary_metadata.tsv -P {db_path}/")

    #### Downloading basis array
    # Non-negative Matrix Factorization (NMF) results, 
    # decomposing the presence/absence matrix (hg38) into 16 components.
    #
    # It contains a 733 row (biosample) x 16 (component) peak presence/abscence 
    # matrix (not a binary matrix)
    #
    # Used later to map component number within metadata dataframe and 
    # find proportion for given component
    os.system(f"wget https://zenodo.org/record/3838751/files/2018-06-08NC16_NNDSVD_Basis.npy.gz -P {db_path}/")
    os.system(f"gzip -d {db_path}/2018-06-08NC16_NNDSVD_Basis.npy.gz")

    #### Downloading matrix array
    # The mixture array contains 3.5M x 16 matrix of 
    # peak presence/absence decomposed into 16 components
    os.system(f"wget https://zenodo.org/record/3838751/files/2018-06-08NC16_NNDSVD_Mixture.npy.gz -P {db_path}/")
    os.system(f"gzip -d {db_path}/2018-06-08NC16_NNDSVD_Mixture.npy.gz")

    # Loading in DHS_Index_and_Vocabulary of ~3.6M DHSs that contains 
    # the following information:
    # seqname, start, end, identifier, mean_signal, numsaples, summit, 
    # core_start, core_end, component
    #
    # Tab-separated file with DNase I Hypersensitive Site (DHS) coordinates,
    # including DHS summits and core regions and assignments to regulatory components.
    os.system(f"wget https://www.meuleman.org/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz -P {db_path}/")
    os.system(f"gunzip -d {db_path}/DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz")
    
    # Downloading binary peak matrix file
    # Presence/absence matrix of DHSs (rows) versus biosamples (columns)
    os.system(f"wget https://zenodo.org/records/3838751/files/dat_bin_FDR01_hg38.txt.gz -P {db_path}/")
    os.system(f"gzip -d {db_path}/dat_bin_FDR01_hg38.txt.gz")
    


def master_dataset() -> None:

    db_path = "./database"
    d_path = "./data"

    #### Genome data
    GENOME_PATH = f"{db_path}/hg38.fa"    
    # Call data_class.py
    genome = ReferenceGenome.from_path(GENOME_PATH)
    #print(genome)

    
    #### DHS metadata
    # Remove the last row who is empty
    DHS_Index_and_Vocabulary_metadata = pd.read_table(
                            f"{db_path}/DHS_Index_and_Vocabulary_metadata.tsv").iloc[:-1]    


    #### Basis matrix
    # Contains a 733 row (biosample) x 16 (component) peak presence/abscence 
    # matrix (not a binary matrix)
    # Used later to map component number within metadata dataframe and 
    # find proportion for given component
    
    # Converting npy file to csv
    basis_array = np.load(f"{db_path}/2018-06-08NC16_NNDSVD_Basis.npy")
    np.savetxt(f"{db_path}/2018-06-08NC16_NNDSVD_Basis.csv", basis_array, delimiter=",")

    # Creating nmf_loadings matrix from csv
    nmf_loadings1 = pd.read_csv(f"{db_path}/2018-06-08NC16_NNDSVD_Basis.csv", header=None)
    nmf_loadings1.columns = ["C" + str(i) for i in range(1, 17)]
    #print(nmf_loadings1.shape)


    # Joining metadata with component presence matrix
    DHS_Index_and_Vocabulary_metadata = pd.concat(
                        [DHS_Index_and_Vocabulary_metadata, nmf_loadings1], 
                        axis=1)    
    #print(DHS_Index_and_Vocabulary_metadata.head(5))

    # Add a "component" column indicates the component number 
    # with the largest presence value
    COMPONENT_COLUMNS = ["C1",    "C2",    "C3",    "C4",    "C5",    "C6",
            "C7",    "C8",    "C9",    "C10",    "C11",    "C12",    "C13", 
            "C14",    "C15",    "C16",]
    DHS_Index_and_Vocabulary_metadata["component"] = (
        DHS_Index_and_Vocabulary_metadata[COMPONENT_COLUMNS].idxmax(axis=1)\
        .apply(lambda x: int(x[1:])))
    #print(DHS_Index_and_Vocabulary_metadata.head(5))


    #### Creating sequence metadata dataframe
    # The mixture array contains 3.5M x 16 matrix of 
    # peak presence/absence decomposed into 16 components

    # Turning npy file into csv
    mixture_array = np.load(f"{db_path}/2018-06-08NC16_NNDSVD_Mixture.npy").T
    np.savetxt(f"{db_path}/2018-06-08NC16_NNDSVD_Mixture.csv", mixture_array, 
                delimiter=",")

    # Creating nmf_loadings matrix from csv and renaming columns
    nmf_loadings2 = pd.read_csv(f"{db_path}/2018-06-08NC16_NNDSVD_Mixture.csv", 
                            header=None, names=COMPONENT_COLUMNS)
    #print(nmf_loadings2.shape)


    # DHS_Index_and_Vocabulary contains 
    # [seqname, start, end, identifier, mean_signal, numsaples, summit, 
    # core_start, core_end, component]

    # Loading sequence metadata
    sequence_metadata = pd.read_table(
                f"{db_path}/DHS_Index_and_Vocabulary_hg38_WM20190703.txt", sep="\t")
    #print(sequence_metadata.shape)

    # Dropping the component column that contains associated tissue 
    # rather than component number 
    # (We will use the component number from DHS_Index_and_Vocabulary_metadata)
    sequence_metadata = sequence_metadata.drop(columns=["component"], axis=1)
    #print(sequence_metadata.shape)

    # Join metadata with component presence matrix
    df = pd.concat([sequence_metadata, nmf_loadings2], axis=1, sort=False)
    #print(df.shape)

    # Recreating some of the columns from our original dataset
    df["component"] = df[COMPONENT_COLUMNS].idxmax(axis=1)\
                        .apply(lambda x: int(x[1:]))
    df["proportion"] = df[COMPONENT_COLUMNS].max(axis=1) / \
                        df[COMPONENT_COLUMNS].sum(axis=1)
    df["total_signal"] = df["mean_signal"] * df["numsamples"]
    df["proportion"] = df[COMPONENT_COLUMNS].max(axis=1) / \
                        df[COMPONENT_COLUMNS].sum(axis=1)
    df["dhs_id"] = df[["seqname", "start", "end", "summit"]]\
                        .apply(lambda x: "_".join(map(str, x)), axis=1)
    df["DHS_width"] = df["end"] - df["start"]

    # Creating sequence column
    df = add_sequence_column(df, genome, 200)
    # Changing seqname column to chr
    df = df.rename(columns={"seqname": "chr"})
    #print(df.shape)

    # Reordering and unselecting columns
    df = df[["dhs_id",        "chr",        "start",        "end",
            "DHS_width",        "summit",        "numsamples",
            "total_signal",        "component",        "proportion",
            "sequence"]]
    #print(df.shape)
    

    # Binary peak matrix file

    # Opening file
    binary_matrix = pd.read_table(f"{db_path}/dat_bin_FDR01_hg38.txt", header=None)
    #print(f"binary_matrix: {binary_matrix.shape}")

    # Collecting names of cells into a list with fromat celltype_encodeID
    celltype_encodeID = [row["Biosample name"] + "_" + row["DCC Library ID"] 
        for _, row in DHS_Index_and_Vocabulary_metadata.iterrows()]

    # Renaming columns using celltype_encodeID list
    binary_matrix.columns = celltype_encodeID
    #print(binary_matrix.head(5))
    #print(binary_matrix.shape)

    # Combining human DHS metadata with binary peak matrix
    # ["dhs_id", "chr", "start", "end", "DHS_width", "summit",
    #  "numsamples", "total_signal", "component", "proportion", 
    #  "sequence", 733 celltype_encodeIDs]
    master_dataset = pd.concat([df, binary_matrix], axis=1, sort=False)
    #print(master_dataset.head(5))
    #print(master_dataset.shape)

    # Save as feather file
    master_dataset.to_feather(f"{d_path}/master_dataset.ftr")




def filter_master(cell_list: List[str]) -> pd.DataFrame:
      
      d_path = "./data"
      
      # ["dhs_id", "chr", "start", "end", "DHS_width", "summit",
      #  "numsamples", "total_signal", "component", "proportion", 
      #  "sequence", 733 celltype_encodeIDs]
      df = pd.read_feather(f"{d_path}/master_dataset.ftr")

      sub_df = FilteringData(df, cell_list=cell_list)\
                        .filter_exclusive_replicates(sort=True, balance=True)

      return sub_df


def main() -> None:

    #### Download all necessary data files    
    #download_data()

    #### Generate a master dataset with information of DHS and biosamples
    master_dataset()

    #### Given a list of samples, filter the master dataset
    cell_list = ["fLeftVentricle_ENCLB231DPY", "fLeftAtrium_ENCLB226FNM", 
                        "fHeart_ENCLB491BID", "fRightVentricle_ENCLB608VQR"]
    sub_data = filter_master(cell_list)
    print(sub_data.head)





if __name__ == "__main__":
    main()      
