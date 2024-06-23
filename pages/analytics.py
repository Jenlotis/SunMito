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
import time
import datetime
import io
import pandas as pd

dash.register_page(__name__, path='/')

def fig_to_url(in_fig, close_all=True, **save_args):
    out_img = BytesIO()
    in_fig.savefig(out_img, format='png', **save_args)
    if close_all:
        in_fig.clf()
        plt.close('all')
    out_img.seek(0)
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)

path_to_directory = "./data"
samples_dir = [dir_names for (dir_path, dir_names, file_names) in os.walk(path_to_directory) if dir_names]

path_to_ref = "./reference"
ref_files = [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_ref) if file_names]

path_to_ss = "./data/long/ss"
ss_files = [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_ss) if file_names]

path_to_cleaned = "./cleaned"

contents = html.Div(children=[
    html.H5("Select the type of data you want to analyse"),
    dcc.Dropdown(
        options = [
            {'label': 'Short Reads', 'value': 'short'},
            {'label': 'Long Reads', 'value': 'long'}
            ],
        id = 'tod_dropdown',
        placeholder = "Select a type of Data you have",
        style = {'white-space':'pre'}
    ),
    html.Br(),
    html.Div(
        id="path_box",
        style = {"display":"none"},
        children=[ 
        html.H5("Choose which Path you want to take(can be multiple)"),
        dcc.Dropdown(
            options = [
                {'label': 'MITObim', 'value': 'mitobim'},
                {'label': 'NOVOPlasty', 'value': 'novoplasty'},
                {"label": "MitoFinder", "value": "mitofinder"}
            ],
            id = 'path_dropdown',
            placeholder = "Select a Path",
            multi = True
        )
    ]),
    html.Div(id="qc_ilu_box",
             style={"display":"none"},
             children=[
        html.Br(),
        html.H5("Quality controll"),
        html.Br(),
        html.P("Choose data to analyze"),
        dcc.Dropdown(
            options=[],
            id='qc_ilu_dropdown',
            placeholder='Select data to analyze',
            multi = True
        ),
        html.Br(),
        html.Button(
            "Run QC",
            id = "qc_ilu_button",
            n_clicks = 0
        )
    ]),
    html.Div(id="merge_button_box",
             style = {"display":"none"},
             children=[
        html.H5("Data merging"),
        html.P(["If yor data is separated beetwen many fastq files ", 
                html.B("click"),
                " here to merge them together to one file."]),
        html.Button(
            "Run merger",
            id = "merg_button",
            n_clicks = 0,
        ),
        html.Br(),
        html.Br()
    ]),
    html.Div(id="qc_nano_box",
             style = {"display":"none"},
             children=[
        html.Br(),
        html.H5("Quality controll"),
        html.P(["If you have seqencing summary file in ",
                html.B("data/nanopor/ss"),
                " folder, then push this button to do QC based on that data"]),
        dcc.Dropdown(
            options=[{'label': f'{ss}', 'value': f"{path_to_ss}/{ss}"} for ss in ss_files[0]],
            id='ss_dropdown',
            placeholder='Select sequencing summary to analyze'),
        html.Br(),
        html.Button(
            "Run good QC",
            id = "qc_good_nano_button",
            n_clicks = 0
        ),
        html.Br(),
        html.Br(),
        html.P(["If you ",
               html.B("don't"),
               " have seqencing summary files push this button to do QC"]),
        dcc.Dropdown(
            options=[{'label': f'{dane}', 'value': f"data/long/{dane}"} for dane in [file_names for (dir_path, dir_names, file_names) in os.walk("data/long") if file_names][0]],
            id='qc_dropdown',
            placeholder='Select data to analyze',
            multi=True),
        html.Br(),
        html.Button(
            "Run FastQC",
            id = "qc_nano_button",
            n_clicks = 0
        )
    ]),
    html.Div(id="trimming_nano_box",
             style={"display":"none"},
             children=[
        html.Br(),
        html.Br(),
        html.H5("Trimming"),
        html.P("Choose data(pozniej bedzie dropdown, teraz trzeba wpisac nazwe: nazwa.fastq.gz"),
        dcc.Dropdown(
            options=[{'label': f'{clean}', 'value': f"data/long/{clean}"} for clean in [file_names for (dir_path, dir_names, file_names) in os.walk("data/long") if file_names][0]],
            id='nano_data_name',
            placeholder='Select data to analyze'),
        html.Br(),
        html.P("Quality treshold"),
        dcc.Input(id="nano_quality",
                  type="number",
                  value=15,
                  placeholder=15),
        html.Br(),
        html.P("Minimal length of the read"),
        dcc.Input(id="nano_minlen",
                  type="number",
                  value=300,
                  placeholder=300),
        html.Br(),
        html.P("Maximal length of the read"),
        dcc.Input(id="nano_maxlen",
                  type="number",
                  value=100000,
                  placeholder=100000),
        html.Br(),
        html.Br(),
        html.Button(
            "Run trimming",
            id = "trimming_nano_button",
            n_clicks = 0
        )
    ]),
    html.Div(id="trim_ilu_box",
             style={"display":"none"},
             children=[
        html.Br(),
        html.H5("Trimming"),
        html.P("Choose sliding window score treshold"),
        dcc.Input(id="ilu_sw_tresh",
                  type="number",
                  value=20,
                  placeholder=20),
        html.Br(),
        html.Br(),
        html.P("Choose minimum length of the sequence"),
        dcc.Input(id="ilu_minlen",
                  type="number",
                  value=60,
                  placeholder=60),
        html.Br(),
        html.Br(),
        html.P("Select two paired files to trim"),
        dcc.Dropdown(
            options=[],
            id='trim_ilu_dropdown',
            placeholder='Select files to trim',
            multi=True),
        html.Br(),
        html.Button(
            "Run trimming",
            id = "trimming_ilu_button",
            n_clicks = 0
        ),
        html.Br(),
        html.Br(),
        html.P("After trimming chosen files will be interlaved to prepare them for MITOBim"),
    ]),
    html.Div(id="mitobim_box",
             style={"display":"none"},
             children=[
        html.Br(),
        html.Br(),
        html.H5("MITOBim"),
        html.Br(),
        html.P("Select bait/reference for MITObim"),
        dcc.Dropdown(
            options=[{'label': f'{ref}', 'value': f"{path_to_ref}/{ref}"} for ref in ref_files[0]],
            id='ref_dropdown',
            placeholder='Select reference to analyze'),
        html.Br(),
        html.P("Select a cleaned sequence for analysis"),
        dcc.Dropdown(
            options=[],
            id='cleaned_dropdown',
            placeholder='Select file to analyze'),
        html.Br(),
        html.P("Choose length of bait sequence"),
        dcc.Input(id='kbait_nano',
                  type="number",
                  placeholder=31,
                  value=31
                  ),
        html.Br(),
        html.Br(),
        html.P("Choose number of iteratons that MITOBim will try to do"),
        dcc.Input(id='iter_nano',
                  type="number",
                  placeholder=10,
                  value=10
                  ),
        html.Br(),
        html.Br(),
        html.Button(
            "Run MITObim",
            id = "mitobim_button",
            n_clicks = 0,
        ) 
    ]),
    html.Div(id="downsampling_box",
             style={"display":"none"},
             children=[
        html.Br(),         
        html.H5("Downsampling"),
        html.P("Calculate to what percentage the cleaned data must be downsampled to fit within the functional limits for the MitoFinder"),
        html.Button(
            "Calculate estimate",
            id = "downsampling_check_button",
            n_clicks = 0
        ),
        html.Br(),
        html.P(id="down_text",
               children=[]
        ),
        html.Br(),
        dcc.Input(id="percent",
                  type="number",
                  placeholder="",
                  value="5",
                  style={"display":"none"}
        ),
        html.Button(
            "Run downsampling",
            id = "downsampling_button",
            n_clicks = 0,
            style={"display":"none"}
        )
    ]),
    html.Div(id="mitfi_box",
             style={"display":"none"},
             children=[
        html.Br(),
        html.H5("MitoFinder"),
        html.P("Select reference to analyze"),
        dcc.Dropdown(
            options=[{'label': f'{ref}', 'value': f"{path_to_ref}/{ref}"} for ref in ref_files[0]],
            id='ref_mitfi_dropdown',
            placeholder='Select reference to analyze'),
        html.Br(),
        html.P("Choose the type of genetic code/organism you have"),
        dcc.Dropdown(
            id="org_mitfi_dropdown",
            placeholder="Choose the type of code you have",
            options=[
                {"label": "The Standard Code", "value": 1},
                {"label": "The Vertebrate Mitochondrial Code", "value": 2},
                {"label": "The Yeast Mitochondrial Code", "value": 3},
                {"label": "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code", "value": 4},
                {"label": "The Invertebrate Mitochondrial Code", "value": 5},
                {"label": "The Ciliate, Dasycladacean and Hexamita Nuclear Code", "value": 6},
                {"label": "The Echinoderm and Flatworm Mitochondrial Code", "value": 9},
                {"label": "The Euplotid Nuclear Code", "value": 10},
                {"label": "The Bacterial, Archaeal and Plant Plastid Code", "value": 11},
                {"label": "The Alternative Yeast Nuclear Code", "value": 12},
                {"label": "The Ascidian Mitochondrial Code", "value": 13},
                {"label": "The Alternative Flatworm Mitochondrial Code", "value": 14},
                {"label": "Chlorophycean Mitochondrial Code", "value": 16},
                {"label": "Trematode Mitochondrial Code", "value": 21},
                {"label": "Scenedesmus obliquus Mitochondrial Code", "value": 22},
                {"label": "Thraustochytrium Mitochondrial Code", "value": 23},
                {"label": "Pterobranchia Mitochondrial Code", "value": 24},
                {"label": "Candidate Division SR1 and Gracilibacteria Code", "value": 25}]
        ),
        html.Br(),
        html.P("Select cleaned paired files to analyze(depending on the size downsompled or not)"),
        dcc.Dropdown(
            options=[],
            id='mitfi_dropdown',
            placeholder='Select file to analyze',
            multi=True),
        html.Br(),
        html.Button(
            "Run MitoFinder",
            id = "mitfi_button",
            n_clicks = 0
        )
    ]),
    html.Div(id="mitobim_ilu_box",
             style={"display":"none"},
             children=[
        html.Br(),
        html.H5("MITOBim"),
        html.P("Select bait/reference for MITObim"),
        dcc.Dropdown(
            options=[{'label': f'{ref}', 'value': f"{path_to_ref}/{ref}"} for ref in ref_files[0]],
            id='ref_ilu_dropdown',
            placeholder='Select reference to analyze'),
        html.Br(),
        html.P("Select one file to analyze, must be interlaved"),
        dcc.Dropdown(
            options=[],
            id='cleaned_ilu_dropdown',
            placeholder='Select file to analyze'),
        html.Br(),
        html.P("Choose length of bait sequence"),
        dcc.Input(id='kbait_ilu',
                  type="number",
                  placeholder=31,
                  value=31
                  ),
        html.Br(),
        html.P("Choose number of iteratons that MITOBim will try to do"),
        dcc.Input(id='iter_ilu',
                  type="number",
                  placeholder=10,
                  value=10
                  ),
        html.Br(),
        html.Br(),
        html.Button(
            "Run MITObim",
            id = "mitobim_ilu_button",
            n_clicks = 0,
        )
    ]), 
    html.Div(id="novopla_box",
             style={"display":"none"},
             children=[
        html.Br(),
        html.H5("NOVOPlasty"),
        html.P(["Choose data to analyze, must be ",
                html.B("paired data,"),
                " two files that have not been cleaned"]),
        dcc.Dropdown(
            options=[{'label': f'{dane}', 'value': f"data/short/{dane}"} for dane in [file_names for (dir_path, dir_names, file_names) in os.walk("data/short/") if file_names][0]],
            id="novo_data_dropdown",
            placeholder="Choose data to analyze",
            multi=True
        ),
        html.Br(),
        html.P("Choose your reference"),
        dcc.Dropdown(
            options=[{'label': f'{ref}', 'value': f"{path_to_ref}/{ref}"} for ref in ref_files[0]],
            id='ref_novo_dropdown',
            placeholder='Select reference to analyze'
        ),
        html.Br(),
        html.Button(
            "Run NOVOPlasty",
            id = "novo_button",
            n_clicks = 0
        )
    ])
])

def run_subprocess(command):
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error running command {command}: {e}")
        return None

@callback(
    Output("path_box", "style"),
    Input("tod_dropdown", "value"),
    prevent_initial_call=True
)
def ilu_starter(tod):
    if tod == "short":
        return {"display":"block"}
    else:
        return {"display":"none"}

@callback(
    Output("merge_button_box", "style"),
    Output("qc_nano_box", "style"),
    Input("tod_dropdown", "value"),
    prevent_initial_call=True 
)
def nano_starter(tod):
    if tod == "long":
        return {"display":"block"}, {"display":"block"}
    else:
        return {"display":"none"}, {"display":"none"}

@callback(
    Output("qc_ilu_box", "style"),
    Output("novopla_box", "style"),
    Output("qc_ilu_dropdown", "options"),
    Input("path_dropdown", "value"),
    prevent_initial_call=True 
)
def path_starter(path):
    if ("mitofinder" in path) or ("mitobim" in path):
        options = [{'label': f'{dane}', 'value': f"data/short/{dane}"} for dane in [file_names for (dir_path, dir_names, file_names) in os.walk("data/short/") if file_names][0]]
        if "novoplasty" in path:
            return {"display":"block"}, {"display":"block"}, options
        else:
            return {"display":"block"}, {"display":"none"}, options
    elif "novoplasty" in path:
        return {"display":"none"}, {"display":"block"}, dash.no_update
    else:
        return {"display":"none"}, {"display":"none"}, dash.no_update
        
@callback(
    Output("trimming_nano_box", "style"),
    Input("qc_nano_button", "n_clicks"),
    State("qc_dropdown","value"),
    prevent_initial_call=True
)
def qc_nano_check(n_clicks, data):
    a=[]
    for i in data:
        subprocess.run(["fastqc", i, "-o", "qc"])
        a.append((i.split("/")[-1]).split(".")[0])
    b = "\n".join(["qc/" + s + "_fastqc.zip" for s in a])
    with open("programs/multiqc_ilu.sh","w") as file:
        file.write(b)
    czas = time.strftime("%d-%m-%Y_%H:%M:%S")
    subprocess.run(["multiqc","-o", f"qc/multiqc_long_{czas}", "-l", "programs/multiqc_ilu.sh"])
    return {"display":"block"}

@callback(
    Output("trim_ilu_box", "style"),
    Output("trim_ilu_dropdown", "options"),
    Input("qc_ilu_button", "n_clicks"),
    State("qc_ilu_dropdown", "value"),
    prevent_initial_call=True
)
def qc_ilu_check(n_clicks, chosen):
    a=[]
    for i in chosen:
        subprocess.run(["fastqc", i, "-o", "qc"])
        a.append((i.split("/")[-1]).split(".")[0])
    b = "\n".join(["qc/" + s + "_fastqc.zip" for s in a])
    with open("programs/multiqc_ilu.sh","w") as file:
        file.write(b)
    czas = time.strftime("%d-%m-%Y_%H:%M:%S")
    subprocess.run(["multiqc","-o", f"qc/multiqc_short_{czas}", "-l", "programs/multiqc_ilu.sh"])
    options = [{'label': f'{dane}', 'value': f"data/short/{dane}"} for dane in [file_names for (dir_path, dir_names, file_names) in os.walk("data/short/") if file_names][0]]
    return {"display":"block"}, options

@callback(
    Output("marge_button", "style"),
    Input("merge_button", "n_clicks"),
    prevent_initial_call=True
)
def nano_one_file(run_id):
    with open("programs/merge_nano.sh", "w") as file:
        file.write(f"zcat data/long/fastq*gz | gzip > {run_id}.fasq.gz")
    subprocess.run(["chmod","+x","programs/merge_nano.sh"])
    subprocess.run(["bash", "./programs/merge_nano.sh"])
    return dash.no_update

@callback(
    Output("trimming_nano_box", "style", allow_duplicate=True),
    Input("qc_good_nano_button", "n_clicks"),
    State("ss_dropdown","value"),
    prevent_initial_call=True
)
def qc_nano_good(n_clicks, seq_sum):
    if "pycoQC" not in run_subprocess(["conda", "env", "list"]).split():
        subprocess.run(["conda", "create", "--name", "pycoQC"])
        subprocess.run(["conda", "activate", "pycoQC"])
        subprocess.run(["pip3", "install", "pycoQC"])
    else:
        subprocess.run(["conda", "init"])
        subprocess.run(["conda", "activate", "pycoQC"])
    czas = time.strftime("%d-%m-%Y_%H:%M:%S")
    subprocess.run(["pycoQC", "-f", seq_sum, "-o",  f"qc/pyco_{czas}.html"])
    subprocess.run(["Rscript", "programs/MinIONQC.R", "-i", seq_sum, "-o", "qc"])
    return {"display":"block"}

@callback(
    Output("mitobim_box","style"),
    Output("cleaned_dropdown", "options"),
    Input("trimming_nano_button","n_clicks"),
    State("nano_data_name", "value"),
    State("nano_quality", "value"),
    State("nano_minlen", "value"),
    State("nano_maxlen", "value"),
    prevent_initial_call=True
)
def clean_nano(n_clicks, name, quality=15, min_len=300, max_len=50000000):
    name_clea=name.split("/")[-1].split(".")[0]
    with open("programs/trim_nano.sh", "w") as file:
        file.write(f"gunzip -c {name} | chopper -q {quality} -l {min_len} --maxlength {max_len} | gzip > cleaned/{name_clea}.cleaned.fastq.gz")
    subprocess.run(["chmod", "+x", "programs/trim_nano.sh"])
    subprocess.run(["bash", "./programs/trim_nano.sh"])
    cleaned_files = [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_cleaned) if file_names]
    options = [{'label': f'{clean}', 'value': f"{path_to_cleaned}/{clean}"} for clean in cleaned_files[0]]
    return {"display":"block"}, options

@callback(
    Output("cleaned_ilu_dropdown", "options"),
    Output("mitobim_ilu_box", "style"),
    Output("downsampling_box", "style"),
    Input("trimming_ilu_button","n_clicks"),
    State("trim_ilu_dropdown", "value"),
    State("path_dropdown", "value"),
    State("ilu_sw_tresh", "value"),
    State("ilu_minlen", "value"),
    prevent_initial_call=True
)
def clean_ilu(n_clicks, dane, path, sw_treshold=20, minlen=60):
    dane_1=(dane[0].split("/")[-1]).split(".")[0]
    dane_2=(dane[1].split("/")[-1]).split(".")[0]
    subprocess.run([
        "TrimmomaticPE",
        dane[0],
        dane[1],
        f"cleaned/{dane_1}_out.fastq.gz",
        f"cleaned/{dane_1}_un.fastq.gz",
        f"cleaned/{dane_2}_out.fastq.gz",
        f"cleaned/{dane_2}_un.fastq.gz",
        f"SLIDINGWINDOW:4:{sw_treshold}",
        f"MINLEN:{minlen}"
    ])
    options = [{'label': f'{clean}', 'value': f"{path_to_cleaned}/{clean}"} for clean in [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_cleaned) if file_names][0]]
    if ("mitofinder" in path) and ("mitobim" in path):
        dane_0=dane_1.split("_")[0]
        subprocess.run(["reformat.sh", f"in1=cleaned/{dane_1}_out.fastq.gz", f"in2=cleaned/{dane_2}_out.fastq.gz", f"out=cleaned/{dane_0}.Out_inter.fastq.gz", "overwrite=true"])
        return options, {"display":"block"}, {"display":"block"}
    elif "mitobim" in path:
        dane_0=dane_1.split("_")[0]
        subprocess.run(["reformat.sh", f"in1=cleaned/{dane_1}_out.fastq.gz", f"in2=cleaned/{dane_2}_out.fastq.gz", f"out=cleaned/{dane_0}.Out_inter.fastq.gz", "overwrite=true"])
        return options, {"display":"block"}, {"display":"none"}
    else:
        return dash.no_update, dash.no_update, {"display":"block"}
        
@callback(
    Output("down_text", "children"),
    Output("percent","style"),
    Output("percent", "placeholder"),
    Output("percent", "value"),
    Output("downsampling_button", "style"),
    Output("mitfi_box", "style"),
    Output("mitfi_dropdown","options"),
    Input("downsampling_check_button", "n_clicks"),
    State("trim_ilu_dropdown", "value"),
    prevent_initial_call=True
)
def downsam_check(path_path_to_directory_with_fasta, dane):
    with open("programs/downsam_check.sh", "w") as file:
        file.write(f"seqkit stats {dane[0]} | awk -v dolari=\"{dane[0]}\" '$1~\"\"dolari\"\" {{print $4}}' | sed 's/,//g' | awk '{{print 7000000/$1*100}}'")
    subprocess.run(["chmod", "+x", "programs/downsam_check.sh"])
    procenty = run_subprocess(["bash", "./programs/downsam_check.sh"])
    # print(procenty)
    if float(procenty) >= 80:
        return f"Calculated percentage {procenty}% informs us that data is small enough for MitoFinder(≤7 000 000)", {"display":"none"}, dash.no_update, dash.no_update, dash.no_update, {"display":"block"}, [{'label': f'{clean}', 'value': f"{path_to_cleaned}/{clean}"} for clean in [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_cleaned) if file_names][0]]
    else:
        return f"Calculated percentage is {procenty}%", {"display":"block"}, round(procenty), round(procenty), {"display":"block"}, dash.no_update, dash.no_update
    
@callback(
    Output("mitfi_box", "style", allow_duplicate=True),
    Output("mitfi_dropdown","options", allow_duplicate=True),
    Input("downsampling_button", "n_clicks"),
    State("percent", "value"),
    State("trim_ilu_dropdown", "value"),
    prevent_initial_call=True
)
def downsam_do(n_clicks, percent, data):
    name=((data[1].split("/")[-1]).split(".")[0]).split("_")[0]
    with open("programs/down.sh", "w") as file:
        file.write(f"""python2 programs/downsample.py -s {percent} --interleave -r {data[0]} -r {data[1]} | gzip > cleaned/{name}_{percent}.fastq.gz
        reformat.sh int=t in=cleaned/{((data[1].split("/")[-1]).split(".")[0]).split("_")[0]}_{percent}.fastq.gz out1=cleaned/{name}.down_pair{percent}.1.fastq.gz out2=cleaned/{name}.down_pair{percent}.2.fastq.gz overwrite=true""")
    subprocess.run(["chmod", "+x", "programs/down.sh"])
    subprocess.run(["bash", "./programs/down.sh"])
    return {"display":"block"}, [{'label': f'{clean}', 'value': f"{path_to_cleaned}/{clean}"} for clean in [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_cleaned) if file_names][0]]
          
@callback(
    Output("mitfi_button", "style"),
    Input("mitfi_button", "n_clicks"),
    State("mitfi_dropdown", "value"),
    State("ref_mitfi_dropdown", "value"),
    State("org_mitfi_dropdown", "value"),
    prevent_initial_call=True
)
def mitfi_pair(n_clicks, data, ref, org):
    data_1=data[0]
    data_2=data[1]
    name=((data[1].split("/")[-1]).split(".")[0]).split("_")[0]
    subprocess.run(["python2", "programs/MitoFinder-master/mitofinder", "-j", name, "-1", data_1, "-2", data_2, "-r", ref, "-o", org, "--override"])
    return dash.no_update

@callback(
    Output("trimming_nano_button", "style"),
    Input("mitobim_button", "n_clicks"),
    State('ref_dropdown', "value"),
    State('cleaned_dropdown', "value"),
    State('kbait_nano', "value"),
    State('iter_nano', "value"),
    prevent_initial_call=True
)
def mitobim_nano(n_clicks, reference, file_name, kbait=31, iterations=10):
    file_name=file_name.split("/")[-1]
    path=os.getcwd()
    with open("programs/make_mitobim.sh", "w") as file:
        file.write(f"sudo docker run -d -it -v {path}/cleaned/:/home/data/input/ -v {path}/output/:/home/data/output/ -v {path}/reference/:/home/data/reference/ chrishah/mitobim /bin/bash")
    subprocess.run(["chmod", "+x", "programs/make_mitobim.sh"])
    subprocess.run(["bash", "./programs/make_mitobim.sh"]) 

    with open("programs/name_mitobim.sh", "w") as file:
        file.write("sudo docker ps | awk '$0 ~ \"chrishah\" {print $1}'")
    subprocess.run(["chmod", "+x", "programs/name_mitobim.sh"])
    container_id=run_subprocess(["bash", "./programs/name_mitobim.sh"])[0]
    reference=reference.split("/")[-1]

    with open("programs/run_mitobim.sh", "w") as file:
        file.write(f"sudo docker exec {container_id} /home/src/scripts/MITObim.pl -sample {file_name} -ref {file_name} -readpool /home/data/input/{file_name} --quick /home/data/reference/{reference} -end {iterations} --kbait {kbait} --clean --redirect_tmp /home/data/output/")
    subprocess.run(["chmod", "+x", "programs/run_mitobim.sh"])
    subprocess.run(["bash", "./programs/run_mitobim.sh"])

    with open("programs/move_mitobim.sh","w") as file:
        file.write(f"sudo docker exec {container_id} cp -r ./iteration* ./data/output/")
    subprocess.run(["chmod", "+x", "programs/move_mitobim.sh"])
    subprocess.run(["bash", "./programs/move_mitobim.sh"])
    
    subprocess.run(["sudo", "docker", "stop", container_id])
    subprocess.run(["sudo", "docker", "rm", container_id])
    return dash.no_update
    
@callback(
    Output("trimming_nano_button", "style", allow_duplicate=True),
    Input("mitobim_ilu_button", "n_clicks"),
    State('ref_ilu_dropdown', "value"),
    State('cleaned_ilu_dropdown', "value"),
    State('kbait_ilu', "value"),
    State('iter_ilu', "value"),
    prevent_initial_call=True
)
def mitobim_ilu(n_clicks, reference, file_name, kbait=31, iterations=10):
    file_name=(file_name[0]).split("/")[-1]
    path=os.getcwd()
    with open("programs/make_mitobim.sh", "w") as file:
        file.write(f"sudo docker run -d -it -v {path}/cleaned/:/home/data/input/ -v {path}/output/:/home/data/output/ -v {path}/reference/:/home/data/reference/ chrishah/mitobim /bin/bash")
    subprocess.run(["chmod", "+x", "programs/make_mitobim.sh"])
    subprocess.run(["bash", "./programs/make_mitobim.sh"]) 

    with open("programs/name_mitobim.sh", "w") as file:
        file.write("sudo docker ps | awk '$0 ~ \"chrishah\" {print $1}'")
    subprocess.run(["chmod", "+x", "programs/name_mitobim.sh"])
    container_id=run_subprocess(["bash", "./programs/name_mitobim.sh"])[0]
    reference=reference.split("/")[-1]

    with open("programs/run_mitobim.sh", "w") as file:
        file.write(f"sudo docker exec {container_id} /home/src/scripts/MITObim.pl -sample {file_name} -ref {file_name} -readpool /home/data/input/{file_name} --quick /home/data/reference/{reference} -end {iterations} --kbait {kbait} --clean --redirect_tmp /home/data/output/")
    subprocess.run(["chmod", "+x", "programs/run_mitobim.sh"])
    subprocess.run(["bash", "./programs/run_mitobim.sh"])

    with open("programs/move_mitobim.sh","w") as file:
        file.write(f"sudo docker exec {container_id} cp -r ./iteration* ./data/output/")
    subprocess.run(["chmod", "+x", "programs/move_mitobim.sh"])
    subprocess.run(["bash", "./programs/move_mitobim.sh"])
    
    subprocess.run(["sudo", "docker", "stop", container_id])
    subprocess.run(["sudo", "docker", "rm", container_id])
    return dash.no_update

@callback(
    Output("novo_button", "style"),
    Input("novo_button", "n_clicks"),
    State("novo_data_dropdown", "value"),
    State("ref_novo_dropdown", "value"),    
    prevent_initial_call=True
)    
def novpla(n_clics, data, reference):
    czas=time.strftime("%d-%m-%Y_%H:%M:%S")
    name=((data[0].split("/")[-1]).split(".")[0]).split("_")[0]
    b = ["Project:",
     "-----------------------",
     f"Project name          = {name}",
     "Type                  = mito",
     "Genome Range          = 12000-22000",
     "K-mer                 = 33",
     "Max memory            = 14",
     "Extended log          = 0",
     "Save assembled reads  = no",
     f"Seed Input            = {reference}"
     "Extend seed directly  = no",
     "Reference sequence    =",
     "Variance detection    =",
     "Chloroplast sequence  =",
     "",
     "Dataset 1:",
     "-----------------------",
     "Read Length           = 151",
     "Insert size           = 300",
     "Platform              = illumina",
     "Single/Paired         = PE",
     "Combined reads        =",
     f"Forward reads         = {data[0]}",
     f"Reverse reads         = {data[1]}",
     "Store Hash            =",
     "",
     "Heteroplasmy:",
     "-----------------------",
     "MAF                   =",
     "HP exclude list       =",
     "PCR-free              =",
     "",
     "Optional:",
     "-----------------------",
     "Insert size auto      = yes",
     "Use Quality Scores    = no",
     f"Output path           = output/{czas}_novo"]
    with open(f"programs/{czas}_Nconfig.txt", "w") as file:
        file.write("\n".join(b))
    subprocess.run(["perl","programs/NOVOplasty-master/NOVOPlasty4.3.1.pl", "-c",f"programs/{czas}_Nconfig.txt"])
    return dash.no_update

def layout():
    return render_layout('SunMito', contents)
