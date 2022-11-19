from os.path import dirname,join
from . import soft_db_path
from . import default_file_structures
import luigi
from tasks.basic_tasks import base_luigi_task
from toolkit import run_cmd, valid_path,get_validate_path
from Bio import SeqIO
import gzip,os


filepath = __file__
input_template_path = join(dirname(filepath),
                           "data_input.template")

from input_parser import fileparser


