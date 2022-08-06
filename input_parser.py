import os

import pandas as pd

from config import input_template_path



class fileparser():
    def __init__(self, filename):
        filename = os.path.abspath(filename)

        self.df = pd.read_csv(filename, sep='\t', index_col=None, dtype=str)

        self.cols, self.df = validate_df(self.df, filename)
        self.df = self.df.set_index("sample ID")
        self.df = self.df.fillna('')
    def get_attr(self, col):
        if col == self.df.index.name:
            return list(self.df.index)
        if col not in self.cols:
            raise Exception("attr %s not in input df" % col)
        else:
            return self.df[col].to_dict()

    @property
    def sid(self):
        return self.get_attr("sample ID")

    @property
    def R1(self):
        return self.get_attr("R1")

    @property
    def R2(self):
        return self.get_attr("R2")


def validate_df(df, filename):
    template_file = input_template_path
    columns_values = open(template_file).read().strip('\n').split('\t')

    if set(columns_values).difference(set(df.columns)):
        missing_cols = set(columns_values).difference(set(df.columns))
        raise Exception("some required columns is missing.  "
                        "%s is missing " % ';'.join(missing_cols))

    if df["sample ID"].duplicated().any():
        raise Exception("sample_name has duplicated.")

    chdir = os.path.dirname(os.path.abspath(filename))
    # os.chdir(chdir)
    # don't use it chdir, it will confuse the abspath and realpath of os.path.
    # print('chdir',chdir)
    for idx, row in df.iterrows():
        # auto implement filepath
        # so easy~~~
        row["R1"] = row["R1"] if (pd.isna(row["R1"]) or os.path.exists(row["R1"])) else os.path.join(chdir, row["R1"])
        row["R2"] = row["R2"] if (pd.isna(row["R2"]) or os.path.exists(row["R2"])) else os.path.join(chdir, row["R2"])
    return columns_values, df
