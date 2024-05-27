import dash
import subprocess
import dash_daq as daq
from .utilities import render_layout  # Ensure this is correctly imported
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
    dcc.Tabs([
        dcc.Tab(children=[
            html.Br(),
            daq.ToggleSwitch(
                id='path_switch',
                value=False,
                label="Short Reads             Long Reads",
                labelPosition='bottom',
                style={'white-space':'pre'}
            ),
            dbc.RadioItems(
                ['Genome', 'Chromosome'], 
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
                html.B('Select a protein:'),
                html.Br(),
                dcc.Dropdown(multi=True,
                             clearable=True,
                             searchable=True,
                             id='gene-select'
                ),
                dcc.Dropdown(multi=True,
                             id='Alignment',
                             style={'display':'none'}
                ),
                dbc.Row([
                    dbc.Col(daq.BooleanSwitch(label="Alignment", labelPosition="top", id='alignment_switch', style={'display': 'none'}, on=False, disabled=False, color='#68CDA3')),
                    dbc.Col(daq.BooleanSwitch(label="Sort", labelPosition="top", id='sort_switch', style={'display': 'none'}, on=False, disabled=True, color='#68CDA3')),
                    dbc.Col(daq.BooleanSwitch(label="Scale", labelPosition="top", id='scale_switch', style={'display': 'none'}, on=False, disabled=True, color='#68CDA3'))
                ]),
                html.Br(),
                html.Div([
                    dbc.Row([
                        dbc.Col([
                            html.B('Score threshold:'),
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
            html.Div(id='align-res'),
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

def run_subprocess(command):
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error running command {command}: {e}")
        return None

def folder_making(path_to_directory_with_fasta):
    os.makedirs(os.path.join(path_to_directory_with_fasta, "cleaned"), exist_ok=True)
    os.makedirs(path_to_directory_with_fasta, exist_ok=True)

def qc_check(path_to_directory_with_fasta):
    run_subprocess(["fastqc", path_to_directory_with_fasta])
    run_subprocess(["multiqc", path_to_directory_with_fasta])

def clean_ilu(path_to_directory_with_fasta, sliding_window, sw_treshold):
    ilu_file_name = run_subprocess(["ls", f"{path_to_directory_with_fasta}/*fastq.gz"]).split("_")[0]
    run_subprocess([
        "TrimmomaticPE",
        f"{path_to_directory_with_fasta}/*1.fastq.gz",
        f"{path_to_directory_with_fasta}/*2.fastq.gz",
        f"{path_to_directory_with_fasta}/cleaned/{ilu_file_name}_1_out.fastq.gz",
        f"{path_to_directory_with_fasta}/cleaned/{ilu_file_name}_1_un.fastq.gz",
        f"{path_to_directory_with_fasta}/cleaned/{ilu_file_name}_2_out.fastq.gz",
        f"{path_to_directory_with_fasta}/cleaned/{ilu_file_name}_2_un.fastq.gz",
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:True",
        f"SLIDINGWINDOW:{sliding_window}:{sw_treshold}"
    ])
    return ilu_file_name

def downsam_s2s(path_path_to_directory_with_fasta, ilu_file_name):
    procenty = run_subprocess(["ls", path_path_to_directory_with_fasta, "|", "grep", ilu_file_name + ".down_" ], capture_output=True, text=True)
    szacunek = run_subprocess(["seqkit", "stats", f"./cleaned/{ilu_file_name}.Out.fasta", "4", "|", "awk", "-v", f"dolari={ilu_file_name}", "'$1~\"./cleaned/\"dolari\".Out1.fastq", "\{print\}'", "|", "sed", "'s/,//g", "|", "awk", "'\{print", "7000000/$1*100}'"])
    run_subprocess(["python2", "./github/MITObim/misc_scripts/downsample.py", "-s", procenty, "-r", f"./cleaned/{ilu_name_fasta}.Out.fasta", "|", "gzip", ">", f"./downsampling/{ilu_file_name}.downsam_{procenty}.fastq.gz"])
    print(procenty)
    return procenty

def mitfi_pair(path_to_directory_with_fasta, procenty, referencja_m, organism):
    ilu_file_name = run_subprocess(["ls", f"{path_to_directory_with_fasta}/*fastq.gz"]).split("_")[0]
    run_subprocess([
      "python2", "./github/MitoFinder/mitofinder", "-j", f"{ilu_file_name}.{procenty}", "-1", f"./downsampling/{ilu_file_name}.down_pair{procenty}.1.fastq.gz", "-2", f"./downsampling/{ilu_file_name}.down_pair{procenty}.2.fastq.gz", "-r", referencja_m, "-o", organism, "--override"])

def mitobim(path_to_directory_with_fasta, reference):
    ilu_file_name = run_subprocess(["ls", f"{path_to_directory_with_fasta}/*fastq.gz"]).split("_")[0]
    
    run_subprocess([
        "sudo", "docker", "run", "-d", "-it",
        "-v", f"{path_to_directory_with_fasta}/{ilu_file_name}/cleaned/:/home/data/input/",
        "-v", f"{path_to_directory_with_fasta}/{ilu_file_name}/output/:/home/data/output/",
        "-v", f"{path_to_directory_with_fasta}/reference/:/home/data/reference/",
        "chrishah/mitobim", "/bin/bash"
    ])
    
    container_id = run_subprocess([
        "sudo", "docker", "ps", "|", "awk", "'$0 ~ \"chrishah\" {print $1}'"
    ])
    
    run_subprocess([
        "sudo", "docker", "exec", container_id,
        "/home/src/scripts/MITObim.pl", "-sample", ilu_file_name, "-ref", ilu_file_name,
        "-readpool", f"/home/data/input/{ilu_file_name}.Out_inter.fastq.gz",
        "--quick", f"/home/data/reference/{reference}", "-end", "10", "--clean",
        "--redirect_tmp", "/home/data/output/"
    ])
    
    run_subprocess([
        "sudo", "docker", "exec", container_id,
        "cp", "-r", "./iteration*", "./data/output/"
    ])
    
    run_subprocess(["sudo", "docker", "stop", container_id])
    run_subprocess(["sudo", "docker", "rm", container_id])

def layout():
    return render_layout('Analytics', contents)
