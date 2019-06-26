# 16s pipelines implemented with qiime2 API

Of course, first you need to install qiime2.
Just follow the official instruction [qiime2](https://qiime2.org/)
Beside install the qiime2, this project also have some external library need to install.

```bash
source activate qiime2-2019.1
pip install biopython
```


just run below code to test:

```bash
source activate qiime2-2019.1
python main.py deblur
# or
python main.py dada2
```

if you want to custom the parameter,you could copy the param.template to any path and modify it.
then pass it to the script like.

```bash
python main.py dada2 -p param.template
```