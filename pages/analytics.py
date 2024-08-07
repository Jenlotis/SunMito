import dash
import subprocess
import dash_daq as daq
from .utilities import render_layout  # Ensure this is correctly imported
import os
from dash import html, dcc, callback, Input, Output, State, ctx, dash_table
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
    html.Div(id="path_box",
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
        ),
        html.Br(),
        #html.Div(html.Img(src="assets/loading.gif", height=100, width=100),
        #         id="qc_ilu"#,style={"display":"none"}
        #),
        html.Br(),
        html.Div(id="html_view_ilu",
                 children=[])
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
        html.Div(html.Img(src="assets/loading.gif", height=100, width=100),
                 id="merge",
                 style={"display":"none"}
        ),
        html.Br(),
        html.Br()
    ]),
    html.Div(id="qc_nano_box",
             style = {"display":"none"},
             children=[
        html.Br(),
        html.H5("Quality controll"),
        html.P("Choose file/-s for QC"),
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
        ),
        html.Br(),
        html.Div(html.Img(src="assets/loading.gif", height=100, width=100),
                 id="qc_nano",
                 style={"display":"none"}
        ),
        html.Div(id="html_view_nano",
                 children=[])
    ]),
    html.Div(id="trimming_nano_box",
             style={"display":"none"},
             children=[
        html.Br(),
        html.Br(),
        html.H5("Trimming"),
        html.P("Choose data to analyze"),
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
        html.Br(),
        html.P("Minimal length of the read"),
        dcc.Input(id="nano_minlen",
                  type="number",
                  value=300,
                  placeholder=300),
        html.Br(),
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
        ),
        html.Div(html.Img(src="assets/loading.gif", height=100, width=100),
                 id="trim_nano",
                 style={"display":"none"}
        ),
    ]),
    html.Div(id="trim_ilu_box",
             style={"display":"none"},
             children=[
        html.Br(),
        html.H5("Trimming"),
        html.P("There will also be step that removes ilumina adapters, it will use deafult values for the programm"),
        html.P("Remove leading bases that are chosen quality(when value is 0, bases won't be removed)"),
        dcc.Input(id="ilu_leading",
                  type="number",
                  value=3,
                  placeholder=3),
        html.Br(),
        html.P("Remove trailing bases that are below chosen quality(when value is 0, bases won't be removed)"),
        dcc.Input(id="ilu_trailing",
                  type="number",
                  value=3,
                  placeholder=3),
        html.Br(),
        html.P("Choose sliding window score treshold"),
        dcc.Input(id="ilu_sw_tresh",
                  type="number",
                  value=15,
                  placeholder=15),
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
        html.Div(html.Img(src="assets/loading.gif", height=100, width=100),
                 id="trim_ilu",
                 style={"display":"none"}
        ),
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
        ),
        html.Div(html.Img(src="assets/loading.gif", height=100, width=100),
                 id="mitobim_nano",
                 style={"display":"none"}
        ),
    ]),
    html.Div(id="downsampling_box",
             style={"display":"none"},
             children=[
        html.Br(),         
        html.H5("Downsampling"),
        #psutil.virtual_memory()[0]/1000000000
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
        html.Br(),
        html.Button(
            "Run downsampling",
            id = "downsampling_button",
            n_clicks = 0,
            style={"display":"none"}
        ),
        html.Div(html.Img(src="assets/loading.gif", height=100, width=100),
                 id="down_do",
                 style={"display":"none"}
        )
    ]),
    html.Div(id="mitfi_box",
             style={"display":"none"},
             children=[
        html.Br(),
        html.P(id="ram_text",
            children=[]
        ),
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
        html.P("Select cleaned paired files to analyze(depending on the size of the data downsampled or not)"),
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
        ),
        html.Div(html.Img(src="assets/loading.gif", height=100, width=100),
                 id="mitfi_ilu",
                 style={"display":"none"}
        ),
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
        html.Br(),
        html.P("Choose number of iteratons that MITObim will try to do"),
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
        ),
        html.Div(html.Img(src="assets/loading.gif", height=100, width=100),
                 id="mitobim_ilu",
                 style={"display":"none"}
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
        ),
        html.Div(html.Img(src="assets/loading.gif", height=100, width=100),
                 id="novopla",
                 style={"display":"none"}
        )
    ]),
    html.Br(),
    html.Br(),
    html.Div(id="empty"),
    html.Br(),
    html.Br(),
])

def run_subprocess(command):
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error running command {command}: {e}")
        return None

@callback(
    Output("qc_ilu","style"),
    Output("merge","style"),
    Output("qc_nano","style"),
    Output("trim_nano","style"),
    Output("trim_ilu","style"),
    Output("mitobim_nano","style"),
    Output("down_che","style"),
    Output("down_do","style"),
    Output("mitfi_ilu","style"),
    Output("mitobim_ilu","style"),
    Output("novopla","style"),
    Input("empty", "style"),
    prevent_initial_call=True
)
def loading_gif(on, box):
    out = []
    states = ["qc_ilu", "merge", "qc_nano", "trim_nano", "trim_ilu", "mitobim_nano", "down_che","down_do", "mitfi_ilu", "mitobim_ilu", "novopla"] #, "mitfi_nano"
    if on == "on":
        for i in states:
            if i == box:
                out.append({"display":"block"})
            else:
                out.append({"display":"none"})
        return out[0], out[1], out[2], out[3], out[4], out[5], out[6], out[7], out[8], out[9],out[10]#, out[11]
    else:
        return {"display":"none"}, {"display":"none"}, {"display":"none"}, {"display":"none"}, {"display":"none"}, {"display":"none"}, {"display":"none"}, {"display":"none"}, {"display":"none"}, {"display":"none"}, {"display":"none"}#, {"display":"none"}

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
    Output("html_view_nano", "children"),
    Input("qc_nano_button", "n_clicks"),
    State("qc_dropdown","value"),
    prevent_initial_call=True
)
def qc_nano_check(n_clicks, data):
    loading_gif("on", "qc_nano")
    a=[]
    for i in data:
        subprocess.run(["fastqc", i, "-o", "qc"])
        a.append((i.split("/")[-1]).split(".")[0])
    b = "\n".join(["qc/" + s + "_fastqc.zip" for s in a])
    with open("programs/multiqc_ilu.sh","w") as file:
        file.write(b)
    czas = time.strftime("%d-%m-%Y_%H.%M.%S")
    subprocess.run(["multiqc","-o", f"qc/multiqc_long_{czas}", "-l", "programs/multiqc_ilu.sh"])
    subprocess.run(["cp", "-r", f"qc/multiqc_long_{czas}", f"assets/multiqc_long_{czas}"])
    loading_gif("off", "qc_nano")
    return {"display":"block"}, html.Iframe(width="100%", height="500" ,src=f"assets/multiqc_long_{czas}/multiqc_report.html")

@callback(
    Output("trim_ilu_box", "style"),
    Output("trim_ilu_dropdown", "options"),
    Output("html_view_ilu", "children"),
    Input("qc_ilu_button", "n_clicks"),
    State("qc_ilu_dropdown", "value"),
    prevent_initial_call=True
)
def qc_ilu_check(n_clicks, chosen):
    loading_gif("on", "qc_ilu")
    a=[]
    for i in chosen:
        subprocess.run(["fastqc", i, "-o", "qc"])
        a.append((i.split("/")[-1]).split(".")[0])
    b = "\n".join(["qc/" + s + "_fastqc.zip" for s in a])
    with open("programs/multiqc_ilu.sh","w") as file:
        file.write(b)
    czas = time.strftime("%d-%m-%Y_%H.%M.%S")
    subprocess.run(["multiqc","-o", f"qc/multiqc_short_{czas}", "-l", "programs/multiqc_ilu.sh"])
    options = [{'label': f'{dane}', 'value': f"data/short/{dane}"} for dane in [file_names for (dir_path, dir_names, file_names) in os.walk("data/short/") if file_names][0]]
    subprocess.run(["cp", "-r", f"qc/multiqc_short_{czas}", f"assets/multiqc_long_{czas}"])
    loading_gif("off", "qc_ilu")
    return {"display":"block"}, options, html.Iframe(width="100%", height="500" ,src=f"assets/multiqc_long_{czas}/multiqc_report.html")

@callback(
    Output("marge_button", "style"),
    Input("merge_button", "n_clicks"),
    prevent_initial_call=True
)
def nano_one_file(run_id):
    loading_gif("on", "merge")
    with open("programs/merge_nano.sh", "w") as file:
        file.write(f"zcat data/long/fastq*gz | gzip > {run_id}.fasq.gz")
    subprocess.run(["chmod","+x","programs/merge_nano.sh"])
    subprocess.run(["bash", "./programs/merge_nano.sh"])
    loading_gif("off", "merge")
    return dash.no_update

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
    loading_gif("on", "trim_nano")
    name_clea=name.split("/")[-1].split(".")[0]
    with open("programs/trim_nano.sh", "w") as file:
        file.write(f"gunzip -c {name} | chopper -q {quality} -l {min_len} --maxlength {max_len} | gzip > cleaned/{name_clea}.cleaned.fastq.gz")
    subprocess.run(["chmod", "+x", "programs/trim_nano.sh"])
    subprocess.run(["bash", "./programs/trim_nano.sh"])
    cleaned_files = [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_cleaned) if file_names]
    options = [{'label': f'{clean}', 'value': f"{path_to_cleaned}/{clean}"} for clean in cleaned_files[0]]
    loading_gif("off", "trim_nano")
    return {"display":"block"}, options

@callback(
    Output("cleaned_ilu_dropdown", "options"),
    Output("mitobim_ilu_box", "style"),
    Output("downsampling_box", "style"),
    Output("ram_text", "children"),
    Output("mitfi_box", "style"),
    Input("trimming_ilu_button","n_clicks"),
    State("trim_ilu_dropdown", "value"),
    State("path_dropdown", "value"),
    State("ilu_sw_tresh", "value"),
    State("ilu_minlen", "value"),
    State("ilu_leading","value"),
    State("ilu_trailing","value"),
    prevent_initial_call=True
)
def clean_ilu(n_clicks, dane, path, sw_treshold=15, minlen=60, leading=3, trailing=3):
    loading_gif("on", "trim_ilu")
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
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True",
        f"LEADING:{leading}",
        f"TRAILING:{trailing}",
        f"SLIDINGWINDOW:4:{sw_treshold}",
        f"MINLEN:{minlen}"
    ])
    
    if ("mitofinder" in path) and ("mitobim" in path):
        dane_0=dane_1.split("_")[0]
        subprocess.run(["reformat.sh", f"in1=cleaned/{dane_1}_out.fastq.gz", f"in2=cleaned/{dane_2}_out.fastq.gz", f"out=cleaned/{dane_0}.out_inter.fastq.gz", "overwrite=true"])
        loading_gif("off", "trim_ilu")
        if float(subprocess.run("free", capture_output=True, text=True).stdout.strip().split()[7])/1024**2 > 16:
            return [{'label': f'{clean}', 'value': f"{path_to_cleaned}/{clean}"} for clean in [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_cleaned) if file_names][0]], {"display":"block"}, dash.no_update, "This computer have more then 16GB of RAM so downsampling is not necessary", {"display":"block"}, 
        else:
            return [{'label': f'{clean}', 'value': f"{path_to_cleaned}/{clean}"} for clean in [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_cleaned) if file_names][0]], {"display":"block"}, {"display":"block"}, dash.no_update, dash.no_update
    elif "mitobim" in path:
        dane_0=dane_1.split("_")[0]
        subprocess.run(["reformat.sh", f"in1=cleaned/{dane_1}_out.fastq.gz", f"in2=cleaned/{dane_2}_out.fastq.gz", f"out=cleaned/{dane_0}.out_inter.fastq.gz", "overwrite=true"])
        loading_gif("off", "trim_ilu")
        return [{'label': f'{clean}', 'value': f"{path_to_cleaned}/{clean}"} for clean in [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_cleaned) if file_names][0]], {"display":"block"}, {"display":"none"}, dash.no_update, dash.no_update,
    else:
        if float(subprocess.run("free", capture_output=True, text=True).stdout.strip().split()[7])/1024**2 > 16:
            return [{'label': f'{clean}', 'value': f"{path_to_cleaned}/{clean}"} for clean in [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_cleaned) if file_names][0]], {"display":"block"}, dash.no_update, "This computer have more then 16GB of RAM so downsampling is not necessary", {"display":"block"}, 
        return dash.no_update, dash.no_update, {"display":"block"}, dash.no_update, dash.no_update
        
@callback(
    Output("down_text", "children"),
    Output("percent","style"),
    Output("percent", "placeholder"),
    Output("percent", "value"),
    Output("downsampling_button", "style"),
    Output("mitfi_box", "style", allow_duplicate=True),
    Output("mitfi_dropdown","options"),
    Input("downsampling_check_button", "n_clicks"),
    State("trim_ilu_dropdown", "value"),
    prevent_initial_call=True
)
def downsam_check(n_clicks, dane):
    dane_1=(dane[0].split("/")[-1]).split(".")[0]
    with open("programs/downsam_check.sh", "w") as file:
        file.write(f"seqkit stats cleaned/{dane_1}_out.fastq.gz | awk -v dolari=\"{dane_1}\" '$1~dolari {{print $4}}' | sed 's/,//g' | awk '{{print 7000000/$1*100}}'")
    subprocess.run(["chmod", "+x", "programs/downsam_check.sh"])
    procenty = float(run_subprocess(["bash", "./programs/downsam_check.sh"]))
    # print(procenty)
    if procenty >= 80:
        loading_gif("off", "down_che")
        return f"Calculated percentage {procenty}% informs us that data is small enough for MitoFinder(≤7 000 000)", {"display":"none"}, dash.no_update, dash.no_update, dash.no_update, {"display":"block"}, [{'label': f'{clean}', 'value': f"{path_to_cleaned}/{clean}"} for clean in [file_names for (dir_path, dir_names, file_names) in os.walk(path_to_cleaned) if file_names][0]]
    else:
        loading_gif("off", "down_che")
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
    loading_gif("on", "down_do")
    name=((data[1].split("/")[-1]).split(".")[0]).split("_")[0]
    with open("programs/down.sh", "w") as file:
        file.write(f"""python2 programs/downsample.py -s {percent} --interleave -r {data[0]} -r {data[1]} | gzip > cleaned/{name}_{percent}.fastq.gz
        reformat.sh int=t in=cleaned/{((data[1].split("/")[-1]).split(".")[0]).split("_")[0]}_{percent}.fastq.gz out1=cleaned/{name}.down_pair_{percent}.1.fastq.gz out2=cleaned/{name}.down_pair_{percent}.2.fastq.gz overwrite=true""")
    subprocess.run(["chmod", "+x", "programs/down.sh"])
    subprocess.run(["bash", "./programs/down.sh"])
    loading_gif("off", "down_do")
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
    loading_gif("on", "mitfi_ilu")
    data_1=data[0]
    data_2=data[1]
    name=((data[1].split("/")[-1]).split(".")[0]).split("_")[0]
    subprocess.run(["cd", "output"])
    subprocess.run(["conda", "activate", "MitFi"])
    subprocess.run(["python", "../programs/MitoFinder-master/mitofinder", "-j", name, "-1", f"../{data_1}", "-2", f"../{data_2}", "-r", f"../{ref}", "-o", org, "--override"])
    subprocess.run(["conda", "deactivate"])
    subprocess.run(["cd", ".."])
    loading_gif("off", "mitfi_ilu")
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
    loading_gif("on", "mitobim_nano")
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
    loading_gif("off", "mitobim_nano")
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
    loading_gif("on", "mitobim_ilu")
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
    loading_gif("off", "mitobim_ilu")
    return dash.no_update

@callback(
    Output("novo_button", "style"),
    Input("novo_button", "n_clicks"),
    State("novo_data_dropdown", "value"),
    State("ref_novo_dropdown", "value"),    
    prevent_initial_call=True
)    
def novpla(n_clics, data, reference):
    loading_gif("on", "novopla")
    czas=time.strftime("%d-%m-%Y_%H-%M-%S")
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
    loading_gif("off", "novopla")
    return dash.no_update

def layout():
    return render_layout('SunMito', contents)
