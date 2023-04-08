import os as _os
import sys as _sys
import pickle as _pickle
import platform as _platform
import subprocess as _subprocess
from Bio import SeqIO as _SeqIO

def get_filelist(path):
    Filelist = []
    for home, dirs, files in _os.walk(path):
        for dirname in dirs:
            Filelist.append(_os.path.join(home, dirname))
    return Filelist

def dependency_check(software_name, subpath=False):
    soft_path = None
    paths = [_os.environ.get('PATH', _os.defpath).split(_os.pathsep)] # current environment PATH
    if subpath==True:
        path1 = get_filelist(_os.path.dirname(_os.path.abspath(_sys.executable))) #path of as module file in GUI
        path2 = get_filelist(_os.path.dirname(_os.path.abspath(__file__))), #path of as module file in CMD
        path3 = get_filelist(_os.path.abspath(_sys.path[0])), #path of as script in CMD
        #print('_sys.executable', _os.path.dirname(_os.path.abspath(_sys.executable)))
        #print('__file__', _os.path.dirname(_os.path.abspath(__file__)))
        #print('_sys.path[0]', _os.path.abspath(_sys.path[0]))
        paths.extend(path1)
        paths.extend(path2)
        paths.extend(path3)
    for path in [j for i in paths for j in i]:
        if _os.path.exists(_os.path.join(path, software_name)):
            soft_path =  _os.path.join(path, software_name)
    if soft_path == None:
        print('Please install ' + software_name + ' software.')
    else:
        return soft_path

def _pickle_fasta(input_file, output_file):
    seq_db = {}
    for record in _SeqIO.parse(input_file, 'fasta'):
        seq_db.setdefault(record.id, record.seq)
    _pickle.dump(seq_db, open(output_file, 'wb'))
    return None

def _makeblastdb(input_file, database_name):
    if _sys.platform.startswith('linux'):
        cmd = [dependency_check('makeblastdb', subpath=True), '-in', input_file,
               '-dbtype', 'nucl', '-title', database_name, '-parse_seqids', '-out', database_name]
    elif _sys.platform.startswith('win32'):
        cmd = [dependency_check('makeblastdb.exe', subpath=True), '-in', input_file, 
               '-dbtype', 'nucl', '-title', database_name, '-parse_seqids', '-out', database_name]
    _subprocess.run(' '.join(cmd), shell=True)
    return None

def _pickle_anno(input_file, output_file):
    db = [None, {}]
    with open(input_file, 'r') as f:
        header = next(f).strip().split('\t')
        db[0] = header[1:]
        for l in f:
            l = l.strip().split('\t')
            db[1].setdefault(l[0],l[1:])
    _pickle.dump(db, open(output_file, 'wb'))
    return None
    
def create_blastdb(cds, annotation, database_name):
    cds = _os.path.abspath(cds)
    annotation = _os.path.abspath(annotation)
    _os.mkdir(database_name)
    _os.chdir(database_name)
    _makeblastdb(cds, database_name)
    _pickle_fasta(cds, database_name+'.cds')
    _pickle_anno(annotation, database_name+'.anno')
    return None

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Create a blast database.', add_help=False, 
                                     epilog='date:2023/03/23 author:guisen chen email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-c','--cds', metavar='str', help='A file of fasta format include of cds seqence.', required=True)
    required.add_argument('-a','--annotation', metavar='str', help='A file of tsv format inclue of annotation information.', required=True)
    required.add_argument('-d','--database_name', metavar='str', help='The name of the generated database.', required=True)
    optional.add_argument('-h','--help',action='help',help='show this help message and exit')
    optional.add_argument('-v','--version',action='version',version='v1.00')
    args = parser.parse_args()
    create_blastdb(args.cds, args.annotation, args.database_name)