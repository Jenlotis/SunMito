import dash
import subprocess
import dash_daq as daq # was needed for switches
from .utilities import render_layout
import os
from dash import html, dcc, callback, Input, Output, State, ctx, dash_table
import pandas as pd
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from io import BytesIO
import base64
import plotly.graph_objects as go
import dash_bootstrap_components as dbc

dash.register_page(__name__, path='/analytics')

def fig_to_url(in_fig, close_all=True, **save_args):
    out_img = BytesIO()
    in_fig.savefig(out_img, format='png', **save_args)
    if close_all:
        in_fig.clf()
        plt.close('all')
    out_img.seek(0)
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)


contents = html.Div(children=[
    dcc.Store(id='raw-df'),
    dcc.Store(id='target-df'),
    dcc.Store(id='trace-dict'),
    dcc.Store(id='annotations'),
    # dcc.Dropdown(
    #     options=[{'label': f'{specimen}'.replace('_', ' '), 'value': f"{data_path}/{specimen}"} for specimen in specimen_dirs[0]],
    #     id='specimen-dropdown',
    #     clearable=False,
    #     placeholder='Select data to analyze'
    # ),
    dcc.Tabs([
        dcc.Tab(children=[
            html.Br(),
            daq.ToggleSwitch(
                id = 'path_switch',
                value = False,
                label = "Short Reads             Long Reads",
                labelPosition='bottom',
                style={'white-space':'pre'}
            ),
            dbc.RadioItems(
                ['Genome', 'Chromosome'], # removed protein view bc it is unnecessary now
                inline=True,
                style={'display': 'none'},
                id='viz-level'
            ),
            html.Div([
                html.B('Select a chromosome:'),
                html.Br(),
                dcc.Dropdown(
                    id='chromosome-select'
                ),
                html.Br(),
            ], style={'display': 'none'}, id='chromosome-select-container'),
            dcc.Graph(id='single-chromosome-graph',
                      style={'display': 'none'}),
            html.Div([
                html.B('Select a protein:'), # merged from protein window
                html.Br(),
                dcc.Dropdown(multi=True, # allows you to choose multiple proteins at the same time
                             clearable=True,
                             searchable=True,
                             id='gene-select'
                ),
                dcc.Dropdown(multi=True, # storage for align proteins
                             id='Aligment',
                             style={'display':'none'}
                ),
                dbc.Row([ # buttons to activate new fucionalities
                    dbc.Col(daq.BooleanSwitch(label="Alignment", labelPosition="top", id='aligment_switch', style={'display': 'none'}, on=False, disabled=False, color='#68CDA3')),
                    dbc.Col(daq.BooleanSwitch(label="Sort", labelPosition="top", id='sort_switch', style={'display': 'none'}, on=False, disabled=True, color='#68CDA3')),
                    dbc.Col(daq.BooleanSwitch(label="Scale", labelPosition="top", id='scale_switch', style={'display': 'none'}, on=False, disabled=True, color='#68CDA3'))
                ]),

                html.Br(),
                html.Div([ # parameters for scoring of alignment
                    dbc.Row([ 
                        dbc.Col([
                            html.B('Score treshold:'),
                            dbc.Input(
                                id='score-threshold',
                                value=0,
                                required=True,
                                style={'width':125}
                            )]),
                        dbc.Col([
                            html.B('Open gap score:'),
                            dbc.Input(
                                id='open-gap',
                                value=-14,
                                required=True,
                                style={'width':125}
                            )]),
                        dbc.Col([
                            html.B('Extend gap score:'),
                            dbc.Input(
                                id='extend-gap',
                                value=-4,
                                required=True,
                                style={'width':125}
                            )])
                        ])
                ], style={'display': 'none'}, id='score-threshold-container'),
               html.Br(),

            ], style={'display': 'none'}, id='gene-select-container'),
            html.Div(id='results'),
            html.Div(id='align-res'), # div for aligned proteins
            html.Br(),
        ], label='Visualization', value='viz'),
        dcc.Tab(children=[
            html.Br(),
            html.B('Select the scope of statistics:'),
            html.Br(),
            dbc.RadioItems(
                ['Full genome', 'Distinct chromosomes'],
                id='stats-source-select',
                inline=True
            ),
            html.Br(),
            html.B('Select chromosome(s) for statistics:'),
            html.Br(),
            dcc.Dropdown(
                id='stats-chromosome-select',
                disabled=True,
                multi=True,
                clearable=True
            ),
            html.Br(),
            dbc.Button(
                'Calculate statistics',
                id='gen-stats-button',
                n_clicks=0,
                class_name='col-2'
            ),
            html.Br(),
            html.Br(),
            html.Div(id='stats-results')
        ], label='Statistics', value='stat')
    ], id='result-tabs', style={'display': 'none'}, value='viz'),
])

def folder_making(path_path_to_directory_with_fasta):
    subprocess.call(["mkdir", path_path_to_directory_with_fasta + "/cleaned"])
    subprocess.call(["mkdir", path_path_to_directory_with_fasta + "/"])

    
    
    
def qc_check(path_path_to_directory_with_fasta):
    subprocess.call(["fasqc", path_path_to_directory_with_fasta])
    subprocess.call(["multiqc", path_path_to_directory_with_fasta])# mozna dodac -o aby miec kontrole na tym gdzie pliki wyladuja
    
    # zrob widoczne okienko ktore bedzie zawierac wyniki z multiqc
    # return

def clean_ilu(path_path_to_directory_with_fasta):
    ilu_file_name = ("_".split(subprocess.call(["ls", path_path_to_directory_with_fasta + "/*fastq.gz"], capture_output=True, text=True)))[0]
    subprocess.call(["TrimmomaticPE",
                     path_path_to_directory_with_fasta + "*1.fastq.gz", path_path_to_directory_with_fasta + "*2.fastq.gz", path_path_to_directory_with_fasta + "/cleaned/" + ilu_file_name + "_1_out.fastq.gz", path_path_to_directory_with_fasta + "/cleaned/" + ilu_file_name + "_1_un.fastq.gz", path_path_to_directory_with_fasta + "/cleaned/" + ilu_file_name + "_2_out.fastq.gz", path_path_to_directory_with_fasta + "/cleaned/" + ilu_file_name + "_2_un.fastq.gz", "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:True", "SLIDINGWINDOW:4:20"])
    return(ilu_file_name)
    
# def downsam_s2s(path_path_to_directory_with_fasta, ilu_file_name):
#     procenty = subprocess.call(["ls", path_path_to_directory_with_fasta, "|", "grep", ilu_file_name + ".down_" ], capture_output=True, text=True)
#     print(procenty)
    

def mitobim(path_path_to_directory_with_fasta):

	# starts a docker(if docker doesn't exist ona a computer crates it) then run MITObim, after mitobim end
    ilu_file_name = ("_".split(subprocess.call(["ls", path_path_to_directory_with_fasta + "/*fastq.gz"], capture_output=True, text=True)))[0]
    p = subprocess.call(["pwd"], capture_output=True, text=True)
    subprocess.call(["sudo", "docker", "run", "-d", "-it", "-v", p + "/" + ilu_file_name + "/cleaned/:/home/data/input/", "-v", p + "/" + ilu_file_name + "/output/:/home/data/output/", "-v", p+"/reference/:/home/data/reference/",  "chrishah/mitobim", "/bin/bash"])

    kontener = subprocess.call(["sudo", "docker", "ps", "|", "awk", "'$0", "~", "\"chrishah\"", "{print" "$1}'"], capture_output=True, text=True)

    subprocess.call(["sudo docker", "exec", kontener, "/home/src/scripts/MITObim.pl", "-sample", ilu_file_name, "-ref", ilu_file_name, "-readpool", "/home/data/input/$i.Out_inter.fastq.gz", "--quick", "/home/data/reference/$REFERENCE_B", "-end", "10", "--clean", "--redirect_tmp", "/home/data/output/"])

    subprocess.call(["sudo", "docker", "exec", kontener, "cp", "-r", "./iteration*", "./data/output/"])

    subprocess.call(["sudo", "docker", "stop", kontener])
    subprocess.call(["sudo", "docker", "rm", kontener])
    
    
    
def layout():
    return render_layout('Analytics', contents)