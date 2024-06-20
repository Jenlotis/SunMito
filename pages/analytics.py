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

path_to_directory = "./data"
samples_dir = [dir_names for (dir_path, dir_names, file_names) in os.walk(path_to_directory) if dir_names]

path_to_ref = "./reference"
ref_files = [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_ref) if file_names]

path_to_ss = "./data/long/ss"
ss_files = [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_ss) if file_names]

path_to_cleaned = "./cleaned"

contents = html.Div(children=[
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
        dcc.Dropdown(
            options=[{'label': f'{specimen}', 'value': f"{path_to_directory}/{specimen}"} for specimen in samples_dir[0]],
            id='sample_dropdown',
            placeholder='Select data folder to analyze',
            style={"display":"none"}
        ),
    ),
    html.Br(),
    html.Div(
        dcc.Dropdown(
            options = [
                {'label': 'MITObim', 'value': 'mitobim'},
                {'label': 'NOVOPlasty', 'value': 'novoplasty'},
                {"label": "MitoFinder", "value": "mitofinder"}
            ],
            id = 'path_dropdown',
            placeholder = "Select a Path",
            style = {"display":"none"},
            multi = True
        )
    ),
    html.Br(),
    html.Div(id="qc_ilu_box",
             style={"display":"none"},
             children=[
        dcc.Dropdown(
            options=[],
            id='qc_ilu_dropdown',
            placeholder='Select data to analyze',
            multi = True
        ),
        html.Button(
            "Run QC",
            id = "qc_ilu_button",
            n_clicks = 0
        )
    ]),
    html.Div(id="merge_button_box",
             style = {"display":"none"},
             children=[
        html.P("If yor data is separated beetwen many fastq files click here to merge them together to one file."),
        html.Button(
            "Run merger",
            id = "merg_button",
            n_clicks = 0,
        )
    ]),
    html.Br(),
    html.Div(id="qc_nano_box",
             style = {"display":"none"},
             children=[
        html.P("If you have seqencing summary file in data/nanopor/ss folder then push this button to do QC based on that data"),
        dcc.Dropdown(
            options=[{'label': f'{ss}', 'value': f"{path_to_ss}/{ss}"} for ss in ss_files[0]],
            id='ss_dropdown',
            placeholder='Select sequencing summary to analyze'),
        html.Button(
            "Run good QC",
            id = "qc_good_nano_button",
            n_clicks = 0
        ),
        html.Br(),
        html.P("If you have don't seqencing summary files push this button to do QC"),
        html.Button(
            "Run FastQC",
            id = "qc_nano_button",
            n_clicks = 0
        )
    ]),
    html.Br(),
    html.Div(id="trimming_nano_box",
             style={"display":"none"},
             children=[
        html.P("Choose data(pozniej bedzie dropdown, teraz trzeba wpisac nazwe: nazwa.fastq.gz"),
        dcc.Input(id="nano_data_name",
                  type="text",
                  placeholder="Data name"),
        html.P("Quality treshold"),
        dcc.Input(id="nano_quality",
                  type="number",
                  value=15,
                  placeholder=15),
        html.P("Minimal length of the read"),
        dcc.Input(id="nano_minlen",
                  type="number",
                  value=300,
                  placeholder=300),
        html.P("Maximal length of the read"),
        dcc.Input(id="nano_maxlen",
                  type="number",
                  value=100000,
                  placeholder=100000),
        html.Button(
            "Run trimming",
            id = "trimming_nano_button",
            n_clicks = 0
        )
    ]),
    html.Br(),
    html.Div(id="trim_ilu_box",
             style={"display":"none"},
             children=[
        dcc.Input(id="ilu_sw_tresh",
                  type="number",
                  value=20,
                  placeholder=20),
        dcc.Input(id="ilu_minlen",
                  type="number",
                  value=60,
                  placeholder=60),
        dcc.Dropdown(
            options=[],
            id='trim_ilu_dropdown',
            placeholder='Select files to trim',
            multi=True),
        html.Button(
            "Run trimming",
            id = "trimming_ilu_button",
            n_clicks = 0
        )
    ]),
    html.Br(),
    html.Div(id="mitobim_box",
             style={"display":"none"},
             children=[
        dcc.Dropdown(
            options=[{'label': f'{ref}', 'value': f"{path_to_ref}/{ref}"} for ref in ref_files[0]],
            id='ref_dropdown',
            placeholder='Select reference to analyze'),
        dcc.Dropdown(
            options=[],
            id='cleaned_dropdown',
            placeholder='Select file to analyze'),
        html.P("Choose length of bait sequence"),
        dcc.Input(id='kbait_nano',
                  type="number",
                  placeholder=31,
                  value=31
                  ),
        html.P("Choose number of iteratons that MITOBim will try to do"),
        dcc.Input(id='iter_nano',
                  type="number",
                  placeholder=10,
                  value=10
                  ),
        html.Button(
            "Run MITObim",
            id = "mitobim_button",
            n_clicks = 0,
        ) 
    ]),
    html.Br(),
    html.Div(id="downsampling_box",
             style={"display":"none"},
             children=[
        html.Button(
            "Calculate estimate",
            id = "downsampling_check_button",
            n_clicks = 0
        ),
        html.P(id="down_text",
               children=[]
        ),
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
    html.Br(),
    html.Div(id="mitfi_box",
             style={"display":"none"},
             children=[
        dcc.Dropdown(
            options=[{'label': f'{ref}', 'value': f"{path_to_ref}/{ref}"} for ref in ref_files[0]],
            id='ref_mitfi_dropdown',
            placeholder='Select reference to analyze'),
        dcc.Dropdown(
            id="org_mitfi_dropdown",
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
        dcc.Dropdown(
            options=[],
            id='mitfi_dropdown',
            placeholder='Select file to analyze',
            multi=True),
        html.Button(
            "Run MitoFinder",
            id = "mitfi_button",
            n_clicks = 0
        )
    ]),
    html.Br(),
    html.Div(id="mitobim_ilu_box",
             style={"display":"none"},
             children=[
        dcc.Dropdown(
            options=[{'label': f'{ref}', 'value': f"{path_to_ref}/{ref}"} for ref in ref_files[0]],
            id='ref_ilu_dropdown',
            placeholder='Select reference to analyze'),
        dcc.Dropdown(
            options=[],
            id='cleaned_ilu_dropdown',
            placeholder='Select file to analyze',
            multi=True),
        html.P("Choose length of bait sequence"),
        dcc.Input(id='kbait_ilu',
                  type="number",
                  placeholder=31,
                  value=31
                  ),
        html.P("Choose number of iteratons that MITOBim will try to do"),
        dcc.Input(id='iter_ilu',
                  type="number",
                  placeholder=10,
                  value=10
                  ),
        html.Br(),
        html.Button(
            "Run MITObim",
            id = "mitobim_ilu_button",
            n_clicks = 0,
        )
    ]),
    html.Br(),
    html.Div(id="novopla_box",
             style={"display":"none"},
             children=[
        dcc.Dropdown(
            options=[{'label': f'{ref}', 'value': f"{path_to_ref}/{ref}"} for ref in ref_files[0]],
            id='ref_novo_dropdown',
            placeholder='Select reference to analyze'),
        html.Button(
            "Run NOVOPlasty",
            id = "novopla_button",
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
    Output("sample_dropdown","style"),
    Input("tod_dropdown", "value"),
    prevent_initial_call=True 
)
def view_starter(tod):
    if tod == "short":
        return {"display":"block"}
    elif tod == "long":
        return {"display":"block"}
    else:
        return {"display":"none"}

@callback(
    Output("path_dropdown", "style"),
    Input("sample_dropdown","value"),
    State("tod_dropdown", "value"),
    prevent_initial_call=True
)
def ilu_starter(dane, tod):
    if tod == "short":
        return {"display":"block"}
    else:
        return {"display":"none"}

@callback(
    Output("merge_button_box", "style"),
    Output("qc_nano_box", "style"),
    Input("sample_dropdown","value"), 
    State("tod_dropdown", "value"),
    prevent_initial_call=True 
)
def nano_starter(folder, tod):
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
    elif path == "novoplasty":
        return {"display":"none"}, {"display":"block"}, dash.no_update
    else:
        return {"display":"none"}, {"display":"none"}, dash.no_update
        
@callback(
    Output("trimming_nano_box", "style"),
    State("sample_dropdown", "value"),
    Input("qc_nano_button", "n_clicks"),
    prevent_initial_call=True
)
def qc_nano_check(path_to_directory_with_fasta,n_clicks):
    for i in run_subprocess(["find", path_to_directory_with_fasta, "-name", "*gz"]).split():
        subprocess.run(["fastqc", i, "-o", "qc"])
    subprocess.run(["multiqc","-o", "qc", "qc"])
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
    subprocess.run(["multiqc","-o", f"qc/multiqc_{czas}", "-l", "programs/multiqc_ilu.sh"])
    options = [{'label': f'{dane}', 'value': f"data/short/{dane}"} for dane in [file_names for (dir_path, dir_names, file_names) in os.walk("data/short/") if file_names][0]]
    return {"display":"block"}, options

@callback(
    Output("marge_button", "style"),
    Input("merge_button", "n_clicks"),
    State("sample_dropdown", "value"),
    prevent_initial_call=True
)
def nano_one_file(path_to_directory_with_fasta, run_id):
    with open("programs/merge_nano.sh", "w") as file:
        file.write(f"zcat {path_to_directory_with_fasta}/fastq*gz | gzip > {run_id}.fasq.gz")
    subprocess.run(["chmod","+x","programs/merge_nano.sh"])
    subprocess.run(["bash", "./programs/merge_nano.sh"])
    return dash.no_update

@callback(
    Output("trimming_nano_box", "style", allow_duplicate=True),
    Input("qc_good_nano_button", "n_clicks"),
    State("sample_dropdown", "value"),
    State("ss_dropdown","value"),
    prevent_initial_call=True
)
def qc_nano_good(n_clicks, path_to_directory_with_fasta, seq_sum):
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
    State("sample_dropdown","value"),
    State("nano_data_name", "value"),
    State("nano_quality", "value"),
    State("nano_minlen", "value"),
    State("nano_maxlen", "value"),
    prevent_initial_call=True
)
def clean_nano(n_clicks, path_to_directory_with_fasta, name, quality=15, min_len=300, max_len=50000000):
    with open("programs/trim_nano.sh", "w") as file:
        file.write(f"gunzip -c {path_to_directory_with_fasta}/{name}.fastq.gz | chopper -q {quality} -l {min_len} --maxlength {max_len} | gzip > cleaned/{name}.cleaned.fastq.gz")
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
        return options, {"display":"block"}, {"display":"block"}, options
    elif "mitobim" in path:
        dane_0=dane_1.split("_")[0]
        subprocess.run(["reformat.sh", f"in1=cleaned/{dane_1}_out.fastq.gz", f"in2=cleaned/{dane_2}_out.fastq.gz", f"out=cleaned/{dane_0}.Out_inter.fastq.gz", "overwrite=true"])
        return options, {"display":"block"}, {"display":"none"}
    else:
        return dash.no_update, dash.no_update, {"display":"block"}, options
        
@callback(
    Output("down_text", "children"),
    Output("percent","style"),
    Output("percent", "placeholder"),
    Output("percent", "value"),
    Output("downsampling_button", "style"),
    Output("mitfi_box", "style"),
    Input("downsampling_check_button", "n_clicks"),
    State("trim_ilu_dropdown", "value"),
    prevent_initial_call=True
)
def downsam_check(path_path_to_directory_with_fasta, dane):
    with open("programs/downsam_check.sh", "w") as file:
        file.write(f"seqkit stats {dane[0]} | awk -v dolari=\"{dane[0]}\" '$1~\"\"dolari\"\" {{print $4}}' | sed 's/,//g' | awk '{{print 7000000/$1*100}}'")
    subprocess.run(["chmod", "+x", "programs/downsam_check.sh"])
    procenty = run_subprocess(["bash", "./programs/downsam_check.sh"])
    print(procenty)
    if float(procenty) >= 80:
        return f"Calculated percentage {procenty}% informs us that data is small enough for MitoFinder(≤7 000 000)", {"display":"none"}, dash.no_update, dash.no_update, dash.no_update, {"display":"block"}
    else:
        return f"Calculated percent is {procenty}%", {"display":"block"}, round(procenty), round(procenty), {"display":"block"}, dash.no_update
    
@callback(
    Output("mitfi_box", "style"),
    Output("mitfi_dropdown","options"),
    Input("downsampling_button", "n_clicks"),
    State("percent", "value"),
    State("trim_ilu_dropdown", "value"),
    prevent_initial_call=True
)
def downsam_do(n_clicks, percent, data):
    name=((data[1].split("/")[-1]).split(".")[0]).split("_")[0]
    with open("programs/down.sh", "w") as file:
        file.write(f"python2 programs/downsample.py -s {percent} --interleave -r {data[0]} -r {data[1]} | gzip > cleaned/{name}_{percent}.fastq.gz
        reformat.sh int=t in=cleaned/{((data[1].split("/")[-1]).split(".")[0]).split("_")[0]}_{percent}.fastq.gz out1=cleaned/{name}.down_pair{percent}.1.fastq.gz out2=cleaned/{name}.down_pair{percent}.2.fastq.gz overwrite=true")
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
    run_subprocess(["python2", "programs/MitoFinder-master/mitofinder", "-j", name, "-1", data_1, "-2", data_2, "-r", ref, "-o", org, "--override"])
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
    
    
def novpla(path_to_directory_with_fasta, reference):
    ilu_file_name = run_subprocess(["ls", f"{path_to_directory_with_fasta}/*fastq.gz"]).split("_")[0]
    b = ["Project:",
     "-----------------------",
     f"Project name          = {ilu_file_name}",
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
     f"Forward reads         = {path_to_directory_with_fasta}/{ilu_file_name}.1.fastq.gz",
     f"Reverse reads         = {path_to_directory_with_fasta}/{ilu_file_name}.2.fastq.gz",
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
     f"Output path           = {path_to_directory_with_fasta}/{ilu_file_name}/",
     "",
     "",
     "Project:",
     "-----------------------",
     "Project name         = Choose a name for your project, it will be used for the output files.",
     "Type                 = (chloro/mito/mito_plant) \"chloro\" for chloroplast assembly, \"mito\" for mitochondrial assembly and",
                           "\"mito_plant\" for mitochondrial assembly in plants.",
     "Genome Range         = (minimum genome size-maximum genome size) The expected genome size range of the genome.",
                            "Default value for mito: 12000-20000 / Default value for chloro: 120000-200000",
                            "If the expected size is know, you can lower the range, this can be useful when there is a repetitive",
                            "region, what could lead to a premature circularization of the genome.",
     "K-mer                = (integer) This is the length of the overlap between matching reads (Default: 33).",
                            "If reads are shorter then 90 bp or you have low coverage data, this value should be decreased down to 23.",
                            "For reads longer then 101 bp, this value can be increased, but this is not necessary.",
     "Max memory           = You can choose a max memory usage, suitable to automatically subsample the data or when you have limited",
                            "memory capacity. If you have sufficient memory, leave it blank, else write your available memory in GB",
                            "(if you have for example a 8 GB RAM laptop, put down 7 or 7.5 (don't add the unit in the config file))",
     "Extended log         = Prints out a very extensive log, could be useful to send me when there is a problem  (0/1).",
     "Save assembled reads = All the reads used for the assembly will be stored in seperate files (yes/no)",
     "Seed Input           = The path to the file that contains the seed sequence.",
     "Extend seed directly = This gives the option to extend the seed directly, in stead of finding matching reads. Only use this when your seed",
                            "originates from the same sample and there are no possible mismatches (yes/no)",
    "Reference (optional) = If a reference is available, you can give here the path to the fasta file.",
                           "The assembly will still be de novo, but references of the same genus can be used as a guide to resolve",
                           "duplicated regions in the plant mitochondria or the inverted repeat in the chloroplast.",
                           "References from different genus haven't beeen tested yet.",
     "Variance detection   = If you select yes, you should also have a reference sequence (previous line). It will create a vcf file",
                            "with all the variances compared to the give reference (yes/no)",
     "Chloroplast sequence = The path to the file that contains the chloroplast sequence (Only for mito_plant mode).",
                            "You have to assemble the chloroplast before you assemble the mitochondria of plants!",
     "",
     "Dataset 1:",
     "-----------------------",
     "Read Length          = The read length of your reads.",
     "Insert size          = Total insert size of your paired end reads, it doesn't have to be accurate but should be close enough.",
     "Platform             = illumina/ion - The performance on Ion Torrent data is significantly lower",
     "Single/Paired        = For the moment only paired end reads are supported.",
     "Combined reads       = The path to the file that contains the combined reads (forward and reverse in 1 file)",
     "Forward reads        = The path to the file that contains the forward reads (not necessary when there is a merged file)",
     "Reverse reads        = The path to the file that contains the reverse reads (not necessary when there is a merged file)",
     "Store Hash           = If you want several runs on one dataset, you can store the hash locally to speed up the process (put \"yes\" to store the hashes locally)",
                             "To run local saved files, goto te wiki section of the github page",
     "",
     "Heteroplasmy:",
     "-----------------------",
     "MAF                  = (0.007-0.49) Minor Allele Frequency: If you want to detect heteroplasmy, first assemble the genome without this option. Then give the resulting",
                            "sequence as a reference and as a seed input. And give the minimum minor allele frequency for this option",
                            "(0.01 will detect heteroplasmy of >1%)",
     "HP exclude list      = Option not yet available",
     "PCR-free             = (yes/no) If you have a PCR-free library write yes",
     "",
     "Optional:",
     "-----------------------",
     "Insert size auto     = (yes/no) This will finetune your insert size automatically (Default: yes)",
     "Use Quality Scores   = It will take in account the quality scores, only use this when reads have low quality, like with the",
                            "300 bp reads of Illumina (yes/no)",
     "Output path          = You can change the directory where all the output files wil be stored."]
    with open(f"{path_to_directory_with_fasta}/{ilu_file_name}_Nconfig.txt", "w") as file:
        file.write("\n".join(b))


def novpla(path_to_directory_with_fasta,ilu_file_name):
    run_subprocess(["perl",
                    f"{path_to_directory_with_fasta}/github/NOVOplasty/NOVOPlasty4.3.1.pl",
                    "-c",
                    f"{path_to_directory_with_fasta}/{ilu_file_name}_Nconfig.txt"])
 
def layout():
    return render_layout('Analytics', contents)
