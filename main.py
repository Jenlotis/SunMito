from dash import Dash, html, dcc
import dash
import dash_bootstrap_components as dbc
import base64


def b64_image(image_filename):
    with open(image_filename, 'rb') as f:
        image = f.read()
    return 'data:image/png;base64,' + base64.b64encode(image).decode('utf-8')

app = Dash(__name__, use_pages=True, suppress_callback_exceptions=True, external_stylesheets=[dbc.themes.BOOTSTRAP])

sidebar = dbc.Nav(
    [
        dbc.NavLink(
            [
                html.Div(name, className='ms-0'),
            ],
            href=page,
            active='exact'
        )
        for (name, page) in [('Home', '/'), ('Analytics', '/analytics'), ('Installation', '/install')]
    ],
    vertical=True,
    pills=True,
    style={'height': '100vh'},
    className='bs-0'
)

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.Div([
                html.Center(
                    html.Img(src=b64_image('./Sun.jpg'),
                             height=120,
                             width=120)
                ),
                html.Br(),
                sidebar
            ])
        ], xs=4, sm=4, md=2, lg=2, xl=2, xxl=2, className='g-0', style={'background-color': '#261300'}),
        dbc.Col([
            dash.page_container
        ], xs=8, sm=8, md=10, lg=10, xl=10, xxl=10, style={'border-left': '2px groove black'})
    ])
], fluid=True, className='bg-light bs-0')

if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8050, debug=False)