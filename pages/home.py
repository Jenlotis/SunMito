import dash
from dash import html, dcc
from .utilities import render_layout

dash.register_page(__name__, path='/')

contents = html.Div(children=[
    html.Div(children='''
        This will be small manual and link to GitHub with whole manual and master thesis
    '''),
])

def layout():
    return render_layout('SunMito', contents)