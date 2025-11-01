import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objs as go
from pyhyperiso.core.Core.HyperisoMaster import PyHyperisoMaster
from pyhyperiso.core.Core.Config import PyConfig, ExternalFlag, Model
from pyhyperiso.core.PhysicalModel.WilsonInterface import PyWilsonInterface, PyWilsonBuildConfig, WGroup, WCoeff, QCDOrder, ContributionType, PyWilsonRequest
from pyhyperiso.core.Core.ParameterSetter import PyParameterSetter, PyParamId, ParameterType
from pyhyperiso.phyperiso.pyhyperiso.core import APIAdapter as _CppAPIAdapter
from pyhyperiso.phyperiso.pyhyperiso.core import APIPath as _CppAPIPath
from pyhyperiso.core.Common.General import PyBlockName, PyLhaID
from pyhyperiso.core.Math.scalar import Scalar
from enum import Enum
from pyhyperiso.core.Core.APIAdapter import PyAPIAdapter
from pathlib import Path

# class PyAPIAdapter:
#     def __init__(self):
#         self._cpp_obj = _CppAPIAdapter()

#     def get_block_infos(self, block_name: str, param_type: ParameterType = ParameterType.SM):
#         raw_block_infos = self._cpp_obj.get_block_infos(PyBlockName(block_name)._cpp_obj, param_type.value)
#         return {PyLhaID(x): Scalar().from_cpp(y).real() for x, y in raw_block_infos.items()}

last_model = Model.THDM
last_lha_path = "lha/testinput_thdm.lha"

# Initialize Hyperiso
hyp = PyHyperisoMaster()
config = PyConfig(
    flags={
        ExternalFlag.IS_LHA_SPECTRUM: False,
        ExternalFlag.HAS_WILSON_INPUT: False,
        ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
        # ExternalFlag.USE_MARTY: False
    },
    model=Model.THDM,
)
hyp.init(lha_file="lha/testinput_thdm.lha", config=config)

interface = PyWilsonInterface()
wilson_config = PyWilsonBuildConfig(groups={WGroup.B, WGroup.BPrime, WGroup.BScalar}, matching_scale=81.0, hadronic_scale=4.18, order=QCDOrder.NNLO)
interface.build(wilson_config)

coefficients = [
    (WGroup.B, [WCoeff[f'C{i}'] for i in range(1, 11)]),
    (WGroup.BPrime, [WCoeff[f'CP{i}'] for i in range(1, 11)] + [WCoeff.CPQ1, WCoeff.CPQ2]),
    (WGroup.BScalar, [WCoeff.CQ1, WCoeff.CQ2]),
    # (WGroup.BPrime, [WCoeff.CPQ1, WCoeff.CPQ2])
]

api_adapter = PyAPIAdapter()

app = dash.Dash(__name__, title="Wilson Dashboard")

app.layout = html.Div(style={'backgroundColor': '#121212', 'color': 'white', 'padding': '20px', 'font-family': 'Arial'}, children=[
    html.H1("📊 Wilson Coefficients & LHA Analysis Dashboard", style={'textAlign': 'center'}),

    html.Div([
        html.Label("Model:"),
        dcc.Dropdown(id='model-selector', options=[
            {'label': 'THDM', 'value': 'THDM'},
            {'label': 'SUSY', 'value': 'SUSY'},
            {'label': 'CUSTOM', 'value': 'CUSTOM'}
        ], value='THDM'),
        dcc.Input(id='custom-model-name', type='text', placeholder='Custom model name', style={'display': 'none'}),
        dcc.Input(id='custom-model-path', type='text', placeholder='Custom model path', style={'display': 'none'}),
        dcc.Input(id='new-lha-path', type='text', value='lha/testinput_thdm.lha', placeholder='Path to LHA file'),

        html.Label("Matching Scale:"),
        dcc.Slider(id='matching-scale', min=50, max=200, step=1, value=81, marks={i: str(i) for i in range(50, 201, 25)}),

        html.Label("Hadronic Scale:"),
        dcc.Slider(id='hadronic-scale', min=1, max=10, step=0.1, value=4.18, marks={i: str(i) for i in range(1, 11)}),

        html.Label("BSM Parameter ID (LHA):"),
        dcc.Input(id='bsm-param-id', type='number', value=37),

        html.Label("BSM Parameter Range:"),
        dcc.RangeSlider(id='bsm-param-range', min=100, max=1000, step=10, value=[200, 500], marks={i: str(i) for i in range(100, 1001, 100)}),

        html.Br(),
        html.Button('Apply', id='apply-button', n_clicks=0)
    ], style={'padding': '15px', 'background-color': '#222', 'border-radius': '10px', 'margin-bottom': '20px'}),

    html.Div([
        html.Div([dcc.Graph(id=f'{group[0].name}-coeff-graph')], style={'width': '48%', 'display': 'inline-block'})
        for group in coefficients
    ] + [
        html.Div([dcc.Graph(id='lha-blocks-graph')], style={'width': '48%', 'display': 'inline-block'}),
        html.Div([dcc.Graph(id='bsm-param-variation-graph')], style={'width': '48%', 'display': 'inline-block'}),
        html.Div([dcc.Graph(id='order-comparison-heatmap')], style={'width': '100%', 'display': 'inline-block'}),
        html.Div([dcc.Graph(id='order-comparison-sm')], style={'width': '48%', 'display': 'inline-block'}),
        html.Div([dcc.Graph(id='order-comparison-bsm')], style={'width': '48%', 'display': 'inline-block'})
    ])
])

@app.callback(
    [Output('custom-model-name', 'style'), Output('custom-model-path', 'style')],
    Input('model-selector', 'value')
)
def toggle_custom_fields(model_value):
    if model_value == 'CUSTOM':
        return {'display': 'block'}, {'display': 'block'}
    return {'display': 'none'}, {'display': 'none'}

@app.callback(
    [Output(f'{group[0].name}-coeff-graph', 'figure') for group in coefficients] +
    [Output('lha-blocks-graph', 'figure'), Output('bsm-param-variation-graph', 'figure'), Output('order-comparison-heatmap', 'figure'),
     Output('order-comparison-sm', 'figure'), Output('order-comparison-bsm', 'figure')],
    Input('apply-button', 'n_clicks'),
    State('model-selector', 'value'),
    State('custom-model-name', 'value'),
    State('custom-model-path', 'value'),
    State('new-lha-path', 'value'),
    State('matching-scale', 'value'),
    State('hadronic-scale', 'value'),
    State('bsm-param-id', 'value'),
    State('bsm-param-range', 'value')
)
def update_graphs(n_clicks, model_value, custom_name, custom_path, lha_path, match_scale, had_scale, bsm_id, bsm_range):
    from pyhyperiso.core.Core.Config import PyConfig
    model_enum = getattr(Model, model_value)

    cfg_args = {'flags': {
        ExternalFlag.IS_LHA_SPECTRUM: False,
        ExternalFlag.HAS_WILSON_INPUT: False,
        ExternalFlag.HAS_TH_OBSERVABLE_INPUT: False,
        # ExternalFlag.USE_MARTY: False
    }, 'model': model_enum}
    if model_enum == Model.CUSTOM:
        cfg_args['mty_model_name'] = custom_name
        cfg_args['mty_model_path'] = Path(custom_path)

    global last_model, last_lha_path

    if model_value != last_model.name or lha_path != last_lha_path:
        print(model_value, lha_path)
        print(last_model, last_lha_path)
        import sys
        sys.exit()
        config = PyConfig(**cfg_args)
        hyp.switch_lha(lha_file=lha_path, config=config)
        last_model = model_value
        last_lha_path = lha_path
    # interface.build(PyWilsonBuildConfig(groups=set([x[0].value for x in coefficients]), matching_scale=match_scale, hadronic_scale=had_scale, order=QCDOrder.NNLO))

    py_set = PyParameterSetter()
    py_set.mutate(PyParamId(ParameterType.WILSON, "B_SCALE", 1), had_scale)
    interface.set_matching_scale(match_scale)
    interface.set_hadronic_scale(had_scale)

    fig_outputs = []
    for group, coeff_list in coefficients:
        bars = []
        for contrib in ContributionType:
            values = [interface.get_FR(PyWilsonRequest(group, coeff, QCDOrder.NNLO, contrib)).real()
                      for coeff in coeff_list]
            bars.append(go.Bar(name=contrib.name, x=[c.name for c in coeff_list], y=values))
        fig = go.Figure(data=bars)
        fig.update_layout(title=f"{group.name} - Coefficients", barmode='group', plot_bgcolor='#222', paper_bgcolor='#121212', font_color='white')
        fig_outputs.append(fig)

    blocks = api_adapter.get_all_blocks()
    block_infos = {str(blk): list(api_adapter.get_block_infos(str(blk), api_adapter.get_type_of_block(str(blk))[0]).keys())[:3] for blk in blocks if blk != ''}
    lha_fig = go.Figure([go.Bar(x=list(block_infos.keys()), y=[len(v) for v in block_infos.values()], marker_color='orange')])
    lha_fig.update_layout(title="LHA Block Summary", plot_bgcolor='#222', paper_bgcolor='#121212', font_color='white')

    bsm_x = list(range(bsm_range[0], bsm_range[1]+1, 10))
    bsm_y = []
    for val in bsm_x:
        py_set.mutate(PyParamId(ParameterType.BSM, "MASS", bsm_id), val)
        bsm_y.append(interface.get_FR(PyWilsonRequest(WGroup.B, WCoeff.C9, QCDOrder.NNLO, ContributionType.BSM)).real())

    bsm_fig = go.Figure(go.Scatter(x=bsm_x, y=bsm_y, mode='lines+markers', line=dict(color='lime')))
    bsm_fig.update_layout(title=f"BSM Variation of C9 (Param {bsm_id})", plot_bgcolor='#222', paper_bgcolor='#121212', font_color='white')

    orders = [QCDOrder.LO, QCDOrder.NLO, QCDOrder.NNLO]
    contribs = [ContributionType.SM, ContributionType.TOTAL, ContributionType.BSM]
    heat_data = []
    for c in contribs:
        row = []
        for o in orders:
            val = interface.get_FR(PyWilsonRequest(WGroup.B, WCoeff.C9, o, c)).real()
            row.append(val)
        heat_data.append(row)

    heatmap_fig = go.Figure(data=go.Heatmap(z=heat_data, x=[o.name for o in orders], y=[c.name for c in contribs], colorscale='Viridis'))
    heatmap_fig.update_layout(title="Order vs Contribution Heatmap for C9", plot_bgcolor='#222', paper_bgcolor='#121212', font_color='white')

    def total_contrib(group, coeff, contrib, max_order_idx):
        orders = [QCDOrder.LO, QCDOrder.NLO, QCDOrder.NNLO]
        return sum(interface.get_sep_order_running(group, coeff, contrib).get(o, Scalar(0)).real() for o in orders[:max_order_idx+1])

    labels = ['LO', 'LO+NLO', 'LO+NLO+NNLO']
    sm_vals = [total_contrib(WGroup.B, WCoeff.C9, ContributionType.SM, i) for i in range(3)]
    bsm_vals = [total_contrib(WGroup.B, WCoeff.C9, ContributionType.BSM, i) for i in range(3)]

    sm_bar = go.Figure([go.Bar(x=labels, y=sm_vals, marker_color='blue')])
    sm_bar.update_layout(title='Order Comparison (SM)', plot_bgcolor='#222', paper_bgcolor='#121212', font_color='white')

    bsm_bar = go.Figure([go.Bar(x=labels, y=bsm_vals, marker_color='red')])
    bsm_bar.update_layout(title='Order Comparison (BSM)', plot_bgcolor='#222', paper_bgcolor='#121212', font_color='white')

    return fig_outputs + [lha_fig, bsm_fig, heatmap_fig, sm_bar, bsm_bar]

if __name__ == '__main__':
    app.run(debug=True, port=8055)