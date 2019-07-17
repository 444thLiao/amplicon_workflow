# 16s pipelines based on luigi

For better handle different pipelines with vsearch/usearch, qiime2-deblur or qiime2-dada2.
We present a unified pipelines which you could perforom different analysis with **single data_input.tab**


# installation

### qiime2
Because this pipelines are embed with qiime2 pipelines, it is better to install a qiime2 environment.

Just follow the official instruction [qiime2](https://qiime2.org/)
like that :
```bash
wget https://data.qiime2.org/distro/core/qiime2-2019.4-py36-linux-conda.yml
conda env create -n qiime2-2019.4 --file qiime2-2019.4-py36-linux-conda.yml
```

### vsearch == 2.6.2

Just follow the official instruction [vsearch](https://github.com/torognes/vsearch)
like that :
```bash
wget https://github.com/torognes/vsearch/archive/v2.6.2.tar.gz
tar xzf v2.6.2.tar.gz
cd vsearch-2.6.2
./autogen.sh
./configure
make
make install # (optional)  # as root or sudo make install
```

### script requires
Requirements could follow `requirements.txt`

Using `pip install -r requirements.txt`

or install its inside the environment of qiime2 
```bash
source active qiime2-2019.4
pip install -r requirements.txt
``` 

### conf
there are multiple config file need to be adjusted.

* conf of pipelines: `config/soft_db_path.py` 
* conf of fastq_screen db: `dir_of_fastq_screen/fastq_screen.conf`

# testing

Within environment of qiime2 (if you don't want to perform qiime2 relative analysis pipelines), you could ignore this.

Just type:

```bash
source activate qiime2-2019.4
python3 main.py test -o ~/temp/16s_testdata --local-scheduler
``` 

It will run all three pipelines with testing data contained at **testset** directory.

Because there have `--local-scheduler`, it will not monitor by `luigid`. More detailed about the **scheduler and luigid**, please follow [luigi Central Scheduler doc](https://luigi.readthedocs.io/en/stable/central_scheduler.html)


## QuickStart

After installation, you need to run testdata first to validate all software and required database is installed. 

When everything is ready, you may have your own **pair-end** sequencing data.

Following the header and separator of `config.data_input.template`, fulfill a new `data_input.tab`.

With this tab, you could run:

```bash
python3 main.py run -- workflow --tab data_input.tab --odir output_dir --analysis-type otu --workers 4 --log-path output_dir/cmd_log.txt
``` 

Besides the params `--tab`, `--odir`, `--analysis-type`, `--log-path`, other params are luigi implemented. 

Here describe a little bit about these params. For more detailed, you should check the documentation of luigi at [luigi doc](https://luigi.readthedocs.io/en/stable/)

* `--tab`: given a path(could be relative/absolute) of `input_data.tab`
* `--odir`: jus the path of output directory. a little be need to say is that, different pipelines like `otu, deblur, dada2`, it will separately located the final output of different pipelines. So **don't worry using same output dir will confuse the result**.
* `--analysis-type`: for now, three options including *otu, deblur, dada2* could be selected, if you want to perform all at once. You could pass `all` param to it. Because there are a lot of overlapped tasks among three different pipelines, it would save a lot of time than running these separately with different `odir`. 
* `--log-path`: it just record the cmd history.*(optional)*
* `--workers`: it could control how many tasks could be parallel.


## about the `input_data.tab`

If you look at the `config.data_input.template`, there are only three header. 

`Tab` is taken as separator of this `input_data.tab` for better handle some weird filename.

Inside this `iniput_data.tab`, you could append more columns besides the necessary `three columns(sample ID	R1	R2)`. This pipelines only check you have these three instead of only have these three.

## Problems

### 1. Could I use barcoded/multiplexed data?
If you have a multiplexed data and you want to use this pipelines. First, you need to use a `demux` tool to demutiplex your data, and fulfill a `data_input.tab` with generated/demutiplexed data.

If you are not sure which **demux tool** to use, you could use our demux pipelines instead of main pipelines.

please following `README.md` at `static` directory.

### 2. Could I tune the parameter of this pipelines?
All parameter have been embedd into a unify directory called `config`. If you want to change the path of software/database, you could see `config/soft_db_path.py`. If you want to change the parameter of otu/deblur/dada2, you could see `config/default_params.py`. 

`config/default_file_structures.py` is not yet used at pipelines, so changing it would not change the result.

### 3. Error raised by dada2?
Maybe because the version of `r-base` in global environment is higher than the installed `r-base` version inside the conda environment. So please united the version of R in global and env.

## Feedback

Please file questions, bugs or ideas 
to the [Issue Tracker](https://github.com/444thLiao/16s_workflow)

## Authors
Tianhua Liao
email: l0404th@gmail.com

 

