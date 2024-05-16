import subprocess
import importlib.util
import dash
import dash_daq as daq
from .utilities import render_layout
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO
import base64
import plotly.graph_objects as go
import dash_bootstrap_components as dbc
import dash
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

dash.register_page(__name__, path='/install')

def install_program(program_name):
    subprocess.call(['sudo', 'apt-get', 'install', program_name])

def is_program_installed(program_name):
    process = subprocess.Popen(['which', program_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return not stderr

def is_library_installed(library_name):
    spec = importlib.util.find_spec(library_name)
    return spec is not None

# List of libraries and programs to check and are needed to run the whole program
libraries_and_programs = [
    ("fastqc",  "Bash Program"),
    ("mawk",    "Bash Program"),
    ("perl",    "Bash Program"),
    ("bbmap",   "Bash Program"),
    ("dash",    "Python Library"),
    ("multiqc", "Python Library"),
]

not_detected = {'packages': [], 'programs': []}

# Iterate over the list and check the installation status
table_data = []
for name, type_ in libraries_and_programs:
    if type_ == "Python Library":
        is_installed = is_library_installed(name)
    else:
        is_installed = is_program_installed(name)

    if not is_installed:
        if type_ == "Python Library":
            not_detected['packages'].append(name)
        else:
            not_detected['programs'].append(name)
        installed = "❌"
    else:
        installed = "✔️"
    
    table_data.append((name, installed, type_))

# Define table layout for detected and not detected packages and programs
def generate_table(data, header):
    rows = [html.Tr([html.Th(col) for col in header])]
    rows += [html.Tr([html.Td(cell) for cell in row]) for row in data]
    return html.Table(rows)

table_header = ["Name", "Installed", "Type"]
detected_table = generate_table(table_data, table_header)

contents = html.Div([
    html.H1("Library and Program Checker"),
    detected_table
])

def layout():
    return render_layout('Install', contents)
