import os,sys
from os.path import abspath
from subprocess import check_call
from glob import glob

def run_cmd(cmd, dry_run=False, log_file=None, **kwargs):
    outstream = None
    if type(log_file) == str:
        if os.path.isfile(log_file):
            if os.path.getsize(log_file) != 0:
                outstream = open(log_file, 'a')
        if outstream is None:
            valid_path(log_file,check_ofile=True)
            outstream = open(log_file, 'w')
    elif log_file is None:
        outstream = sys.stdout
    else:
        outstream = log_file

    executable = "/usr/bin/zsh"
    if not os.path.exists(executable):
        executable = "/bin/bash"
    print(cmd, file=outstream)
    outstream.flush()
    if not dry_run:
        check_call(cmd,
                   shell=True,
                   executable=executable,
                   stdout=outstream,
                   stderr=outstream,
                   **kwargs)
        outstream.flush()

def get_validate_path(pth):
    if not pth.startswith('/'):
        pth = './' + pth
    pth = abspath(pth)
    return pth

def valid_path(in_pth,
               check_size=False,
               check_dir=False,
               check_glob=False,
               check_odir=False,
               check_ofile=False):
    if type(in_pth) == str:
        in_pths = [in_pth]
    else:
        in_pths = in_pth[::]
    for in_pth in in_pths:
        in_pth = os.path.abspath(os.path.realpath(in_pth))
        if in_pth is None:
            continue
        if check_glob:
            query_list = glob(in_pth)
            if not query_list:
                raise Exception('Error because of input file pattern %s' % in_pth)
        if check_dir:
            if not os.path.isdir(in_pth):
                raise Exception("Error because %s doesn't exist" % in_pth)
        if check_size:
            if os.path.getsize(in_pth) <= 0:
                raise Exception("Error because %s does not contain content." % in_pth)
        if check_odir:
            if not os.path.isdir(in_pth):
                os.makedirs(in_pth, exist_ok=True)
        if check_ofile:
            odir_file = os.path.dirname(in_pth)
            if not os.path.isdir(odir_file):
                os.makedirs(odir_file, exist_ok=True)
    return True

def get_dir_path(path,num=1):
    path = os.path.abspath(os.path.realpath(path))
    for _ in range(num):
        path = os.path.dirname(path)
    if path == "/":
        raise Exception("reach root path....")
    return path
