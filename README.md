# DHS-data-analysis

DNase I hypersensitive sites (DHSs) are crucial genetic markers of regulatory DNA regions, harboring significant genetic information linked to cardiovascular diseases (CVDs). Given the profound impact of CVDs on global health, there is a pressing need for a deeper understanding of the underlying genetic and regulartory mechanisms.

<!-- blank line -->
<br>
<!-- blank line -->

The Meuleman team recently developed a comprehensive common coordinate system for DHSs, integrating data from 733 human biosamples including 438 cell and tissue types and states. This effort resulted in the numerical indexing of ~3.6 million DHSs, providing a valuable resource for regulatory DNA analysis (https://www.meuleman.org/research/dhsindex/).

<!-- blank line -->
<br>
<!-- blank line -->

This repository is created to facilitate the processing of DHSs data and the application of advanced AI models for predicting cell- and tissue-specific regulatory elements implicated in the pathogenesis of CVDs. By leveraging these cutting-edge techniques, we aim to uncover insights into the regulatory landsape underlying CVDs, ultimately contributing to the development of novel therapeutic interventions and diagnostic strategies.

## Usage

### DHS data processing

First, clone this repository into your local computer.

```
git clone git@github.com:XuejianXiong/DHS-data-analysis.git
```

Second, create "database" and "data" folders if not exist.

```
cd DHS-data-analysis
mkdir database
mkdir data
```

Next, run the main script **process_DHS_data.py**. 

```
python process_DHS_data.py
```

If any required python library, e.g. pandas, is missing, install it using pip3.

```
pip3 install pandas
``` 


**process_DHS_data.py** includes three functions as follows,

1) download_data(): download all necessary DHSs data files.
2) master_dataset(): generate a master dataset with information of ~3.6M DHS and 733 biosamples.
3) filter_master(): given a list of samples, filtering the master dataset to collect useful information for further analysis.
<!-- blank line -->
<br>
<!-- blank line -->
In addition, there are two other python scripts as follows,

1) data_class.py: define objects for setting data source, getting reference genome, and fitering master dataset.
2) create_seq_columns.py: define functions used to create columns in the sequence matrix.
<!-- blank line -->
<br>
<!-- blank line -->


### DNABERT model

DNABERT, a pre-trained tranasformer deep learning model (https://github.com/jerryji1993/DNABERT/tree/master?tab=readme-ov-file), is used to study the DHS dataset. 

First, create and activate a virtual environment named "dnabert".

```
conda create -n dnabert python=3.6
conda activate dnabert
```

Second, install the package and other python libraries.

```
conda install pytorch torchvision cudatoolkit=10.0 -c pytorch
git clone https://github.com/jerryji1993/DNABERT
cd DNABERT
python3 -m pip install --editable .
cd examples
python3 -m pip install -r requirements.txt
```

Next, download pre-train data "DNABERT6" (if set kmer=6).

```
cd PATH_TO_THE_PRETRAINED_MODEL

manually download from 
https://drive.google.com/file/d/1BJjqb5Dl2lNMg2warsFQ0-Xvn1xxfFXC/view?usp=sharing

unzip 6-new-12w-0.zip
```

Next, Fine-tune with pre-trained model.

```
cd examples
chmod +x fine_tune.sh
.fine_tune.sh
```

Next, after the model is fine-tuned, we an run predictions as follows,
```
cd examples
chmod +x prediction.sh
.prediction.sh
```

### Issues

DNABERT requires GPU resources and is optimized for high-performace computing (HPC), particularly when dealing with large-scale datasets. However, due to the current unavailability of adequate computational resources, we are unable to present any results from DNABERT modeling at this time.


## Results

Given a list of cell types related to CVDs, the number of exclusive peaks in each sample is as follows:

| Cell type | Number of exclusive pearks|
| ---- | ---- |
| Cell type: fHeart, Replicate: ENCLB491BID | Number of exclusive peaks: 19715 |
| Cell type: fLeftAtrium, Replicate: ENCLB226FNM | Number of exclusive peaks: 10466 |
| Cell type: fLeftVentricle, Replicate: ENCLB231DPY | Number of exclusive peaks: 23797 |
| Cell type: fRightVentricle, Replicate: ENCLB608VQR | Number of exclusive peaks: 15412 |
