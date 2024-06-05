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
samples_dirs = [dir_names for (dir_path, dir_names, file_names) in os.walk(path_to_directory) if dir_names]
contents = html.Div(children=[
    html.Br(),
    dcc.Dropdown(
        options = [
            {'label': 'Short Reads', 'value': 'short'},
            {'label': 'Long Reads', 'value': 'long'}
            ],
        id = 'tod_dropdown',
        placeholder = "Select a type of Data you have",
        style = {'white-space':'pre'}
    ),
    
    dcc.Dropdown(
        options=[{'label': f'{sample}'.replace('_', ' '), 'value': f"{path_to_directory}/{sample}"} for sample in samples_dirs[0]],
        id='sample-dropdown',
        placeholder='Select a samplee to analyze'
    ),
    
    html.Div(
        dcc.Dropdown(
            options = [
                {'label': 'MITObim', 'value': 'mitobim'},
                {'label': 'NOVOPlasty', 'value': 'novoplasty'},
                {"label": "MitoFinder", "value": "mitofinder"}
            ],
            id = 'path_dropdown',
            placeholder = "Select a Path",
            style={"display":"none"}
        )
    ),
    
    html.Br(),
    html.Button(
        "Run NOVOPlasty",
        id = "novopla_button",
        n_clicks = 0,
        style={"display":"none"}
    ),
    
    html.Br(),
    html.Button(
        "Run QC",
        id = "qc_ilu_button",
        n_clicks = 0,
        style={"display":"none"}
    ),
    
    html.Br(),
    html.Button(
        "Run trimming",
        id = "trimming_ilu_button",
        n_clicks = 0,
        style={"display":"none"}
    ),
    
    html.Br(),
    html.Button(
        "Run MITObim",
        id = "mitobim_button",
        n_clicks = 0,
        style={"display":"none"}
    ),
    
    html.Br(),
    html.Button(
        "Run downsampling",
        id = "downsampling_button",
        n_clicks = 0,
        style={"display":"none"}
    ),
    
    html.Br(),
    html.Button(
        "Run MitoFinder",
        id = "mitofi_button",
        n_clicks = 0,
        style={"display":"none"}
    ),
    
    html.Br(),
    html.Button(
        "Run merger",
        id = "merg_button",
        n_clicks = 0,
        style={"display":"none"}
    ),
    
    html.Br(),
    html.Button(
        "Run QC",
        id = "qc_nano_button",
        n_clicks = 0,
        style={"display":"none"}
    ),
 

    
    html.Br(),
    html.Button(
        "Run trimming",
        id = "trimming_nano_button",
        n_clicks = 0,
        style={"display":"none"}
    )
    
    # html.Br(),
    # html.Button(
    #     "Run",
    #     id = "_button",
    #     n_clicks = 0,
    #     style={"display":"none"}
    # )

])

def run_subprocess(command):
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error running command {command}: {e}")
        return None

def qc_check(path_to_directory_with_fasta):
    run_subprocess(["fasqc", path_to_directory_with_fasta])
    run_subprocess(["multiqc", path_to_directory_with_fasta])

def nano_one_file(path_to_directory_with_fasta, run_id):
    run_subprocess(["zcat", f"{path_to_directory_with_fasta}/fasq*gz", "|", "gzip", ">", f"{run_id}.fasq.gz"])

def qc_nano_pyco(path_to_ss):
    if "pycoQC" not in run_subprocess(["conda", "env", "list"]).split():
        run_subprocess(["conda", "create", "--name", "pycoQC"])
        run_subprocess(["conda", "activate", "pycoQC"])
        run_subprocess(["pip3", "install", "pycoQC"])
    else:
        run_subprocess(["conda", "activate", "pycoQC"])
    czas = time.strftime("%d-%m-%Y_%H:%M:%S")
    run_subprocess(["pycoQC", "-f", path_to_ss, "-o",  f"qc/pyco_{czas}.html"])

def qc_nano_min(path_to_ss):
    script_directory = os.path.dirname(os.path.abspath(__file__))
    run_subprocess(["Rscript", f"{script_directory}/../programs/MinIONQC.R", "-i", path_to_ss, "-o", f"{script_directory}/../qc"])

def clean_nano(path_to_directory_with_fasta, name, quality=15, min_len=300, max_len=50000000):
    run_subprocess(["gunzip", "-c", f"{name}.fastq.gz", "|", "chopper", "-q", quality, "-l", min_len, "--maxlength", max_len, "|", "gzip", ">", f"{name}.cleaned.fastq.gz"])

def clean_ilu(path_to_directory_with_fasta, sliding_window=4, sw_treshold=20):
    ilu_file_name = run_subprocess(["ls", f"{path_to_directory_with_fasta}/*fastq.gz"]).split("_")[0]
    run_subprocess([
        "TrimmomaticPE",
        f"{path_to_directory_with_fasta}/*1.fastq.gz",
        f"{path_to_directory_with_fasta}/*2.fastq.gz",
        f"{path_to_directory_with_fasta}/cleaned/{ilu_file_name}_1_out.fastq.gz",
        f"{path_to_directory_with_fasta}/cleaned/{ilu_file_name}_1_un.fastq.gz",
        f"{path_to_directory_with_fasta}/cleaned/{ilu_file_name}_2_out.fastq.gz",
        f"{path_to_directory_with_fasta}/cleaned/{ilu_file_name}_2_un.fastq.gz",
        f"SLIDINGWINDOW:{sliding_window}:{sw_treshold}"
    ])
    return ilu_file_name

def downsam_s2s(path_path_to_directory_with_fasta, ilu_file_name):
    procenty = run_subprocess(["ls", path_path_to_directory_with_fasta, "|", "grep", ilu_file_name + ".down_" ], capture_output=True, text=True)
    print(procenty)
    return procenty

def mitfi_pair(path_to_directory_with_fasta):
    ilu_file_name = run_subprocess(["ls", f"{path_to_directory_with_fasta}/*fastq.gz"]).split("_")[0]
    run_subprocess([
      "python2", "./github/MitoFinder/mitofinder", "-j", f"{ilu_file_name}.$XXX", "-1", f"./downsampling/{ilu_file_name}.down_pair$XXX.1.fastq.gz", "-2", f"./downsampling/{ilu_file_name}.down_pair$XXX.2.fastq.gz", "-r", "$REFERENCE_M", "-o", "$ORGANISM", "--override"])

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
