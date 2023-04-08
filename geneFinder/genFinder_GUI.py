#!/usr/bin/env python
# coding: utf-8

import os
import sys
import textwrap
import genFinder
from Bio import SeqIO
import PySimpleGUI as sg

def get_database_dict():
    database_dict = {}
    path1 = os.path.join(os.path.dirname(os.path.abspath(sys.executable)), 'database')
    #path2 =_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'database')
    path3 =os.path.join(os.path.abspath(sys.path[0]), 'database')
    for path in [path1, path3]:
        if os.path.exists(path):
            for root, folders, files in os.walk(os.path.realpath(path)):
                for folder in folders:
                    database_dict.setdefault(folder, os.path.join(root,folder))
    return database_dict

def get_desktop_path(prefix=None):
    if prefix == None:
        return os.path.join(os.path.expanduser("~"), 'Desktop')
    else:
        return os.path.join(os.path.expanduser("~"), 'Desktop', prefix)
    
def comlist(list1, list2):
    tmp = []
    for index, value in enumerate(list1):
        if list2[index] > value:
            tmp.append(list2[index])
        else:
            tmp.append(value)
    return tmp

def table(file, width=4):
    with open(file, 'r') as f:
        max_list = list(map(len, next(f).strip().split('\t')))
        for l in f:
            l = list(map(len,l.strip().split('\t')))
            max_list = comlist(max_list, l)
    max_list = [i + width for i in max_list]
    text = ''        
    with open(file, 'r') as f:
        for l in f:
            l = l.strip().split('\t')
            line = ''
            for index, value in enumerate(l):
                string = '{:<'+ str(max_list[index]) + '}'
                line += string.format(value)
            text += line + '\n'
    return text

def fasta(file, width=70):
    text = ''
    for record in SeqIO.parse(file, 'fasta'):
        r  = '>' + record.id + '\n' + textwrap.fill(str(record.seq), width=width) + '\n'
        text += r
    return text


sg.theme('Default1')
database_dict = get_database_dict()
layout_right = [[sg.Text('文件类型', size=(10,0)), sg.Radio('prot', 'R1', key='prot', default=True), sg.Radio('nucl', 'R1', key='nucl')], 
                [sg.Text('数据库', size=(10,0)), sg.Listbox(list(database_dict.keys()), key='database', size=(28,5), default_values=list(database_dict.keys())[0]), 
                 sg.Button('run', key='run', size=10)]]
layout_left =  [[sg.Text('E值', size=(10,0)), sg.InputText(key='evalue', size=30, default_text='1e-5')], 
                [sg.Text('线程', size=(10,0)), sg.InputText(key='thread', size=30, default_text='1')],
                [sg.Text('最大匹配数', size=(10,0)), sg.InputText(key='max_target_seqs', size=30, default_text='1')],
                [sg.Text('输出文件名', size=(10,0)), sg.InputText(key='output_prefix', size=30)],
                [sg.Text('输出文件夹', size=(10,0)), sg.InputText(key='output_folder', default_text=get_desktop_path(), size=(30,0)), 
                 sg.FolderBrowse(button_text='open', size=10)],
                [sg.Text('输入文件', size=(10,0)), sg.InputText(key='input_file',enable_events=True, size=(30,0)), 
                 sg.FileBrowse(button_text='open', size=10)]
               ]

layout_args = sg.Frame(layout=[[sg.Col(layout_left), sg.Text("",size=10), sg.Col(layout_right)]],
                       title='参数', size=(1000, 200), title_color='#0A59F7')

layout_out = sg.Frame(layout=[[sg.Text('annotion information', size=(72,0), font=("宋体", 10)), sg.Text('blast result', size=(70,0), font=("宋体", 10))],
          [sg.Multiline(key='annotion', size=(70, 20), font=("宋体", 10), horizontal_scroll=True),sg.Multiline(key='blast_out', size=(70, 20), font=("宋体", 10), horizontal_scroll=True)],
          [sg.Text('cds seqence', size=(72,0), font=("宋体", 10)), sg.Text('protein seqence', size=(70,0), font=("宋体", 10))],
          [sg.Multiline(key='cds', size=(70, 20), font=("宋体", 10), horizontal_scroll=True), sg.Multiline(key='pep', size=(70, 20), font=("宋体", 10), horizontal_scroll=True)], ],
         title='输出', size=(1000, 650), title_color='#0A59F7')

layout = [[layout_args], [layout_out],
          [sg.Text('Author: Guisen Chen; Version: 1.0; Email: thecgs001@foxmail.com; Datetime: 2023/04/01', size=(140,0), font=("Arial", 10))]]


window = sg.Window('genFinder GUI', layout, size=(1032,895))

while True:
    event, values = window.read()
    if event == 'input_file':
        prefix = os.path.splitext(os.path.basename(values['input_file']))[0]
        window['output_prefix'].update(prefix)
    if event == 'run':
        if values['nucl']:
            stype = 'nucl'
        elif values['prot']:
            stype = 'prot'
        try:
            genFinder.main(values['input_file'], 
                             os.path.join(values['output_folder'], prefix),
                             database_dict[values['database'][0]],
                             values['evalue'],
                             stype,
                             values['thread'],
                             values['max_target_seqs'])
            annotion = os.path.join(values['output_folder'], values['output_prefix']) + '.annoation_genFinder.xls'
            blast_out = os.path.join(values['output_folder'], values['output_prefix']) + '.blast_out_genFinder.xls'
            cds = os.path.join(values['output_folder'], values['output_prefix']) + '.cds_genFinder.fasta'
            pep = os.path.join(values['output_folder'], values['output_prefix']) + '.pep_genFinder.fasta'
            window['annotion'].update(table(annotion))
            window['blast_out'].update(table(blast_out))
            window['cds'].update(fasta(cds))
            window['pep'].update(fasta(pep))
        except:
            print('check input file')
    if event == None:
        break
window.close()

