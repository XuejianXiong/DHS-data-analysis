# DHS-data-analysis

DNase I hypersensitive sites (DHSs) serve as valuable genetic markers of regulatory DNA regions and contain important genetic information essentially associated with cardiovascular diseases (CVDs). CVDs represent a significant global health burdan, requiring a deeper understanding of the underlying genetic and regulartory mechanisms.

<!-- blank line -->
<br>
<!-- blank line -->

Recently, the Meuleman team generated a comprehensive common coordinate system of DHSs for regulatory DNA (https://www.meuleman.org/research/dhsindex/). It integrated 733 human biosamples including 438 cell and tissue types and states, and numerically indexed ~3.6 million DHSs. 

<!-- blank line -->
<br>
<!-- blank line -->

This repository is created in order to process DHSs data, apply advanced AI models for predicting cell- and tissue-specific regulatory elemetns implicated in CVD pathogenesis.


## Usage

The main script in this repository is **process_DHS_data.py**. It includes three functions:

1) download_data()
2) master_dataset()
3) filter_master()
<!-- blank line -->
<br>
<!-- blank line -->

In addition, there are two other python scripts:

1) data_class.py
2) create_seq_columns.py

<!-- blank line -->
<br>
<!-- blank line -->

Give a list of cell types related to CVDs, the number of exclusive peaks are as follows:

| Cell type | Number of exclusive pearks|
| ---- | ---- |
| Cell type: fHeart, Replicate: ENCLB491BID | Number of exclusive peaks: 19715 |
| Cell type: fLeftAtrium, Replicate: ENCLB226FNM | Number of exclusive peaks: 10466 |
| Cell type: fLeftVentricle, Replicate: ENCLB231DPY | Number of exclusive peaks: 23797 |
| Cell type: fRightVentricle, Replicate: ENCLB608VQR | Number of exclusive peaks: 15412 |
