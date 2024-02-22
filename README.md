# DHS-data-analysis

DNase I hypersensitive sites (DHSs) serve as valuable genetic markers of regulatory DNA regions and contain important genetic information essentially associated with cardiovascular diseases (CVDs). CVDs represent a significant global health burdan, requiring a deeper understanding of the underlying genetic and regulartory mechanisms.

<!-- blank line -->
<br>
<!-- blank line -->

Recently, the Meuleman team generated a comprehensive common coordinate system of DHSs for regulatory DNA (https://www.meuleman.org/research/dhsindex/). It integrated 733 human biosamples including 438 cell and tissue types and states, and numerically indexed ~3.6 million DHSs. 

<!-- blank line -->
<br>
<!-- blank line -->

This repository is created in order to process DHSs data, and apply advanced AI models for predicting cell- and tissue-specific regulatory elements implicated in CVD pathogenesis.


## Usage

### DHS data processing

First, clone this repository into your local computer:

```
git clone git@github.com:XuejianXiong/DHS-data-analysis.git
```

Second, create "database" and "data" folders if not exist:

```
cd DHS-data-analysis
mkdir database
mkdir data
```

Next, run the main script **process_DHS_data.py**. 

```
python process_DHS_data.py
```

If any required python library, e.g. pandas, is missing, install it using pip3

```
pip3 install pandas
``` 

**process_DHS_data.py** includes three functions:

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


### DNABERT model

DNABERT, a pre-trained tranasformer deep learning model (https://github.com/jerryji1993/DNABERT/tree/master?tab=readme-ov-file), is used to study the DHS dataset. 

First, create and activate a virtual environment named dnabert:

```
conda create -n dnabert python=3.6
conda activate dnabert
```

Second, install the package and other python libraries:
```

conda install pytorch torchvision cudatoolkit=10.0 -c pytorch
git clone https://github.com/jerryji1993/DNABERT
cd DNABERT
python3 -m pip install --editable .
cd examples
python3 -m pip install -r requirements.txt
```

Next, download pre-train data "DNABERT6" (if set kmer=6) from
```
cd PATH_TO_THE_PRETRAINED_MODEL
https://drive.google.com/file/d/1BJjqb5Dl2lNMg2warsFQ0-Xvn1xxfFXC/view?usp=sharing
unzip 6-new-12w-0.zip
```

Next, Fine-tune with pre-trained model
```
cd examples
export KMER=6
export MODEL_PATH=PATH_TO_THE_PRETRAINED_MODEL
export DATA_PATH=sample_data/ft/$KMER
export OUTPUT_PATH=./ft/$KMER

python run_finetune.py \
    --model_type dna \
    --tokenizer_name=dna$KMER \
    --model_name_or_path $MODEL_PATH \
    --task_name dnaprom \
    --do_train \
    --do_eval \
    --data_dir $DATA_PATH \
    --max_seq_length 100 \
    --per_gpu_eval_batch_size=32   \
    --per_gpu_train_batch_size=32   \
    --learning_rate 2e-4 \
    --num_train_epochs 5.0 \
    --output_dir $OUTPUT_PATH \
    --evaluate_during_training \
    --logging_steps 100 \
    --save_steps 4000 \
    --warmup_percent 0.1 \
    --hidden_dropout_prob 0.1 \
    --overwrite_output \
    --weight_decay 0.01 \
    --n_process 8
```

Next, after the model is fine-tuned, we an run predictions as follows:
```
export KMER=6
export MODEL_PATH=./ft/$KMER
export DATA_PATH=sample_data/ft/$KMER
export PREDICTION_PATH=./result/$KMER

python run_finetune.py \
    --model_type dna \
    --tokenizer_name=dna$KMER \
    --model_name_or_path $MODEL_PATH \
    --task_name dnaprom \
    --do_predict \
    --data_dir $DATA_PATH  \
    --max_seq_length 75 \
    --per_gpu_pred_batch_size=128   \
    --output_dir $MODEL_PATH \
    --predict_dir $PREDICTION_PATH \
    --n_process 48
```

DNABERT requires GPU and recommend HPC for large-scale datasets. Lacking the necessary computer resources currently, no results can be shown at and after fine-tune step.


## Results

Give a list of cell types related to CVDs, the number of exclusive peaks are as follows:

| Cell type | Number of exclusive pearks|
| ---- | ---- |
| Cell type: fHeart, Replicate: ENCLB491BID | Number of exclusive peaks: 19715 |
| Cell type: fLeftAtrium, Replicate: ENCLB226FNM | Number of exclusive peaks: 10466 |
| Cell type: fLeftVentricle, Replicate: ENCLB231DPY | Number of exclusive peaks: 23797 |
| Cell type: fRightVentricle, Replicate: ENCLB608VQR | Number of exclusive peaks: 15412 |
