## static function documentation

This directory contains a lot of implemented script under before pipelines. So if you are not familiar with this whole project, please do not use it. Or use it carefully following this doc.


##  Demultiplexing script

For easy handle demultiplexing the sequencing data with barcode,we implement a script for performing demultiplexing.

#### header of metadata
For demultiplex data, you must prepare a tab separator file named metadata. The header of metadata should unified because we need to parse 5 or more columns to extract the **barcode** and **primer** information of each sample. 

The header of the metadata should like below:
```python
id_col = 'SampleID'
fb_col = 'Forward_Barcode'
rb_col = 'Reverse_Barcode'
fp_col = 'Forward_Primer'
rp_col = 'Reverse_Primer'
```

Of course, you could also change the default setting at `config/default_params.py` to fit your metadata.

#### usage 

```python
python3 static/toolkit/demux.py --help
```
see help documentation of demux.py

```python
python3 static/toolkit/demux.py -m metadata.tsv -o output_directory -r1 "path_of_R1" -r2 "path_of_R2" -f -p 5
```

the `-r1` param could accept wildcard input, and it will auto coordinate the **r1** and **r2**. But if it failed to coordinate 

`-f` could force to overwrite the before output.
`-p` could control the thread which use to perform this demultiplexing. 0 for use all thread of this computer.

After demupltiexing data, it will generate all splitted data at assigned **output directory** and it will generate a stats of this process. Called `demux_stats.csv`.

For join `static\toolkit\demux.py` and `main.py`, of course it will generate a input metadata formatted as  required. It could directly used for run `main.py`


#### test set

```python
python3 static/toolkit/demux.py -m testset/seq_data2/metadata.tab -o ~/test_demux/ -r1 "testset/seq_data2/*_1.fastq.gz" -r2 "testset/seq_data2/*_2.fastq.gz" -f -p 5
```

using above command to perform test-set validation.