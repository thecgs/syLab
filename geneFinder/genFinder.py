import os as _os
import sys as _sys
import pickle as _pickle
import platform as _platform
import subprocess as _subprocess
from Bio import Seq as _Seq
from Bio import SeqIO as _SeqIO

def get_filelist(path):
    Filelist = []
    for home, dirs, files in _os.walk(path):
        for dirname in dirs:
            Filelist.append(_os.path.join(home, dirname))
    return Filelist

def database_check(index):
    _stat = False
    index_path = None
    if _os.path.exists(index):
        index = _os.path.abspath(index)
        prefix = _os.path.split(index)[1]
        index_path = _os.path.join(index, prefix)
        for suffix in ['.anno','.cds','.nhr','.nin','.nog','.nsd','.nsi','.nsq']:
            if not _os.path.exists(_os.path.join(index, prefix + suffix)):
                print('Missing file '+ prefix + suffix +' please check database.')
                _stat = True
    else:
        Built_in_database_path = [_os.path.join(_os.path.dirname(_os.path.abspath(_sys.executable)), 'database', index),
                                  _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'database', index),
                                  _os.path.join(_os.path.abspath(_sys.path[0]), 'database', index)
                                 ]
        for path in Built_in_database_path:
            if _os.path.exists(path):
                _stat = False
                prefix = _os.path.split(index)[1]
                index_path = _os.path.join(path, prefix)
                for suffix in ['.anno','.cds','.nhr','.nin','.nog','.nsd','.nsi','.nsq']:
                    if not _os.path.exists(_os.path.join(path, prefix + suffix)):
                        print('Missing file '+ prefix + suffix +' please check database.')
                        _stat = True
    if index_path ==None:
        _stat = True
        print('Database is not exists, please check.')
    if _stat == True:
        _sys.exit()
    return index_path

def dependency_check(software_name, subpath=False):
    soft_path = None
    paths = [_os.environ.get('PATH', _os.defpath).split(_os.pathsep)] # current environment PATH
    if subpath==True:
        path1 = get_filelist(_os.path.dirname(_os.path.abspath(_sys.executable))) #path of as module file in GUI
        #path2 = get_filelist(_os.path.dirname(_os.path.abspath(__file__))), #path of as module file in CMD
        path3 = get_filelist(_os.path.abspath(_sys.path[0])), #path of as script in CMD
        #print('_sys.executable', _os.path.dirname(_os.path.abspath(_sys.executable)))
        #print('__file__', _os.path.dirname(_os.path.abspath(__file__)))
        #print('_sys.path[0]', _os.path.abspath(_sys.path[0]))
        paths.append(path1)
        #paths.extend(path2)
        paths.extend(path3)
    for path in [j for i in paths for j in i]:
        if _os.path.exists(_os.path.join(path, software_name)):
            soft_path =  _os.path.join(path, software_name)
    if soft_path == None:
        print('Please install ' + software_name + ' software.')
    else:
        return soft_path
    
def run_blast(query, prefix, db, evalue, stype, thread, max_target_seqs):
    if _sys.platform.startswith('linux'):
        if stype == 'prot':
            cmd = [dependency_check('tblastn', subpath=True), '-query', query, '-db', db,
                   '-num_threads', str(thread), '-evalue', evalue,
                   '-outfmt', '6', '-max_target_seqs', str(max_target_seqs),
                   '-out', prefix+'.blast_out_genFinder.xls']
        if stype == 'nucl':
            cmd = [dependency_check('blastn', subpath=True), '-query', query, '-db', db,
                   '-num_threads', str(thread), '-evalue', evalue,
                   '-outfmt', '6', '-max_target_seqs', str(max_target_seqs),
                   '-out', prefix+'.blast_out_genFinder.xls']
    elif _sys.platform.startswith('win32'):
        if stype == 'prot':
            cmd = [dependency_check('tblastn.exe', subpath=True), '-query', query, '-db', db,
                   '-num_threads', str(thread), '-evalue', evalue,
                   '-outfmt', '6', '-max_target_seqs', str(max_target_seqs),
                   '-out', prefix+'.blast_out_genFinder.xls']
        if stype == 'nucl':
            cmd = [dependency_check('blastn.exe', subpath=True), '-query', query, '-db', db,
                   '-num_threads', str(thread), '-evalue', evalue,
                   '-outfmt', '6', '-max_target_seqs', str(max_target_seqs),
                   '-out', prefix+'.blast_out_genFinder.xls']
    _subprocess.run(' '.join(cmd), shell=True)
    print('Software of blast running completed!')
    return None

def translate(input_file, output_file):
    out = open(output_file, 'w')
    for record in _SeqIO.parse(input_file, 'fasta'):
        out.write('>'+record.id+'\n'+str(_Seq.Seq(record.seq).translate())+'\n')
    out.close()
    return None

def get_unique(prefix):
    unique = {}
    with open(prefix+'.blast_out_genFinder.xls', 'r') as f:
        for l in f:
            l = l.strip().split('\t')
            if l[1] not in unique:
                unique.setdefault(l[1], {l[0]})
            else:
                unique[l[1]].add(l[0])
    return unique
    
def get_seq(prefix, db):
    db = _pickle.load(open(db+'.cds', 'rb'))
    out_cds = open(prefix+'.cds_genFinder.fasta', 'w')
    unique = get_unique(prefix)
    for key in unique:
        out_cds.write('>'+key+'\n'+str(db[key])+'\n')
    out_cds.close()        
    translate(prefix+'.cds_genFinder.fasta', prefix+'.pep_genFinder.fasta')
    print('Seqence extraction completed!')
    return None

def get_anno(prefix, db):
    db = _pickle.load(open(db+'.anno', 'rb'))
    unique = get_unique(prefix)
    out_anno = open(prefix+'.annoation_genFinder.xls', 'w')
    out_anno.write('\t'.join(['Target ID', 'Query ID'] + db[0])+'\n')
    for key in unique:
        out_anno.write(key + '\t' + '|'.join(unique[key]) +'\t'+ '\t'.join(db[1][key])+'\n')
    out_anno.close()
    print('Anntation information extraction completed!')
    return None

def add_header(prefix):
    blast_header = ['Query ID', 'Target ID', 'Percent of Identity', 'Length', 'Misss Match', 'Gap Open', 
                    'Query Start', 'Query End', 'Target Start', 'Target End', 'E-value', 'Bit Score']
    with open(prefix+'.blast_out_genFinder.xls', 'r') as f:
        out = f.read()
    with open(prefix+'.blast_out_genFinder.xls', 'w') as f:
        f.write('\t'.join(blast_header)+'\n')
        f.write(out)
    return None

def main(query, prefix, database_name, evalue, stype, thread, max_target_seqs):
    db = database_check(database_name)
    run_blast(query, prefix, db, evalue, stype, thread, max_target_seqs)
    get_seq(prefix, db)
    get_anno(prefix, db)
    add_header(prefix)
    return None

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Blast genes base on database.', add_help=False, 
                                     epilog='date:2023/03/23 author:guisen chen email:thecgs001@foxmail.com')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i','--input_file', metavar='str', help='A file of fasta format.', required=True)
    required.add_argument('-o','--output_prefix', metavar='str', help='Prefix output file.', required=True)
    required.add_argument('-d','--database_name', metavar='str', help='The name of the generated database.', required=True)
    required.add_argument('-s','--seqence_type', metavar='str', choices=['nucl','prot'], help='Type of input seqence.', required=True)
    optional.add_argument('-e','--evalue', metavar='str', default='1e-5', help='Alignmemt evalue. default: 1e-5')
    optional.add_argument('-t','--thread', metavar='int', default=1, type=int, help='Alignmemt thread number. default: 1')
    optional.add_argument('-m','--max_target_seqs', metavar='int', default=1, help='Alignmemt max target seqence number. default: 1')
    optional.add_argument('-h','--help',action='help',help='show this help message and exit')
    optional.add_argument('-v','--version',action='version',version='v1.00')
    args = parser.parse_args()
    main(args.input_file, args.output_prefix, args.database_name, args.evalue, args.seqence_type, args.thread, args.max_target_seqs)