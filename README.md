# DHS-data-analysis

DNase I hypersensitive sites (DHSs) serve as valuable genetic markers of regulatory DNA regions and contain important genetic information essentially associated with cardiovascular diseases (CVDs). CVDs represent a significant global health burdan, requiring a deeper understanding of the underlying genetic and regulartory mechanisms.

<!-- blank line -->
<br>
<!-- blank line -->

Recently, the Meuleman team generated a comprehensive common coordinate system of DHSs for regulatory DNA (https://www.meuleman.org/research/dhsindex/). It integrated 733 human biosamples including 438 cell and tissue types and states, and numerically indexed ~3.6 million DHSs. 

<!-- blank line -->
<br>
<!-- blank line -->

This repository is created in order to process DHSs data, and apply advanced AI models for predicting cell- and tissue-specific regulatory elemetns implicated in CVD pathogenesis.


## Usage

The main script in this repository is **process_DHS_data.py**. It includes three functions:

1) download_data(): download all necessary DHSs data files.
2) master_dataset(): generate a master dataset with information of ~3.6M DHS and 733 biosamples.
3) filter_master(): given a list of samples, filtering the master dataset to collect useful information for further analysis.
<!-- blank line -->
<br>
<!-- blank line -->

In addition, there are two other python scripts:

1) data_class.py: define objects for setting data source, getting reference genome, and fitering master dataset.
2) create_seq_columns.py: define functions used to create columns in the sequence matrix.

<!-- blank line -->
<br>
<!-- blank line -->

Give a list of cell types related to CVDs, the number of exclusive peaks are as follows:

| Cell type | Number of exclusive pearks|
| ---- | ---- |
| fHeart, Replicate: ENCLB491BID | 19715 |
| fLeftAtrium, Replicate: ENCLB226FNM | 10466 |
| fLeftVentricle, Replicate: ENCLB231DPY | 23797 |
| fRightVentricle, Replicate: ENCLB608VQR | 15412 |
