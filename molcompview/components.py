import dash_bootstrap_components as dbc
from dash import dcc
from dash import html


from molcompview.actions import generate_figure_from_data

def header_alt():
    return dbc.Navbar(
        dbc.Col(
            html.H3("MolCompassView: Visual analysis of chemical space and AD for QSAR/QSPR modelling",
                    className="display-5, text-center, text-white",
                    style={"white-space": "nowrap", "text-align":"center"}),
            width={"size": 6, "offset": 3},
        ),
        color="#0063a6ff",
    )
# def header_alt():
#     return dbc.Navbar(
#         [
#             dbc.Col(html.Img(src = "assets/univie_white.svg", height= "40px"), width={"size": "auto"}),
#             dbc.Col(
#                 dbc.Row(
#                     dbc.Col(
#                         html.H3("MolCompassView: Visual analysis of chemical space and AD for QSAR/QSPR modelling", className="display-5, text-center, text-white",style={"white-space": "nowrap"}),
#                     ),
#                     justify="center",
#                 ),
#             width={"size": 6, "offset": 3},
#             ),
#         ],  color="#0063a6ff",
#         dark=False)


def manual_layout():
    return dbc.Row([
        dbc.Col([
            html.H3("Manual"),
            html.P("This is the manual"),
        ], width=10, align="center")
    ], id="manual-layout-content")

def select_property_dropdown(property_options):
    return dbc.Col(
        dbc.Row([
            dbc.Label("Select property:",width='auto'),
            dbc.Col(
                dcc.Dropdown(
                    id="molcompass-select-property-dropdown",
                    options={f:f for f in property_options},
                )
            )
        ])
    , width=6, align="center", className='g-2 p-2')



def molcompass_figure():
    wrapper = dbc.Col([
        dcc.Graph(id='molcompass-graph', figure=generate_figure_from_data(), style={'height': '90vh', 'width': '98vw'}),
        dcc.Tooltip(id="molcompass-graph-tooltip")
    ], width={"size": 10}, align="center")
    return wrapper

def range_selector(min=0, max=1, step=0.01, value=[0, 1]):
    #Get vertical height of graph
    return html.Div(dcc.RangeSlider(
            id='molcompass-range-slider',
            min=min,
            max=max,
            step=step,
            value=value,
            # marks={i: str(i) for i in range(min, max, step)},
            marks=None,
            vertical=True,
            # verticalHeight=800,
        #Caluclate the height of the slider from the height of the figure
        ),
       id='molcompass-range-slider-container',
       style={'visibility':'hidden'})

def analysis_layout():
    return dbc.Offcanvas(id="analysis-layout", is_open=False, placement="end")

def molcompass_layout(selectable_columns):
    #Make a layout for the molcompass figure with molcompass figure and range slider
    # return dbc.Row([
    #     dbc.Row(select_property_dropdown(selectable_columns), id="molcompass-show-property-dropdown", justify="center"),
    #             dbc.Row([
    #                 [dbc.Col(molcompass_figure()),dbc.Col(range_selector(), width=1),analysis_layout()],
    #         ],justify='center')  # Tabs content
    # ],id="molcompass-layout-content", className='g-0')
    #
    return dbc.Row([
        dbc.Row(select_property_dropdown(selectable_columns), id="molcompass-show-property-dropdown", justify="center"),
        dbc.Row([
            dbc.Col(molcompass_figure(), width=10),
            dbc.Col(range_selector(), width=1),
            analysis_layout()
        ],justify='center')  # Tabs content
    ],id="molcompass-layout-content", className='g-0')



