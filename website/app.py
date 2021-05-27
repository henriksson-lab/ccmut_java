#!/usr/bin/env python3
import sqlite3
import json
import flask
import re
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import re
import time
from plotly.validators.scatter.marker import SymbolValidator
import math

import dash_table

from dash.exceptions import PreventUpdate

import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots

import pandas as pd
import numpy as np


##################################################################################################################
# loading files
##################################################################################################################

print("\n========= loading files =========\n")


#################
# List of strains
conn = sqlite3.connect("file:data/alignments.sqlite?mode=ro", uri=True)
list_strains = pd.read_sql_query("SELECT DISTINCT strain from alignment ORDER BY strain", conn)["strain"].tolist()
conn.close()
dropdown_strains_sorted = [ {'label': s, 'value':s} for s in list_strains ]

#################
# List of genes that we have analyzed
conn = sqlite3.connect("file:data/alignments.sqlite?mode=ro", uri=True)
list_wbid = pd.read_sql_query("SELECT DISTINCT wbid from alignment ORDER BY strain", conn)["wbid"].tolist()
conn.close()

#################
# Mapping wbid - gene symbol
data_genelist = pd.read_csv('data/map_wbid_genesym.csv',sep="\t")[["wbid","symbol"]]
data_genelist = data_genelist[data_genelist["wbid"].isin(list_wbid)]

print("genes in database:")
print(data_genelist.shape)

## List of unique genes
all_genes = data_genelist['symbol'].unique()

#################
# create the gene dropdown for genes with known names (ie. name != ensembl_id)
# data_genelist is: {"ensg":"gene symbol"}
genes_dropdown_options = { v[0]: v[1] for v in data_genelist[ ['wbid', 'symbol'] ].values }
genes_dropdown_options_inv = {v: k for k, v in genes_dropdown_options.items()}

dropdown_genesymbols_sorted = [ {'label': genes_dropdown_options[key], 'value':key} for key in genes_dropdown_options ]
dropdown_genesymbols_sorted = sorted(dropdown_genesymbols_sorted, key=lambda t : t["label"])



# add genes dictionary for genes display (this is a duplicate)   .... TODO, we should not have this twice ideally.
# Genes dict is: {"ensg":"gene symbol - long name"}  TODO was. kicked out the loop that made it
genes_dict = genes_dropdown_options

## Load coordinate system names
coordinate_list_data = pd.read_csv("data/list_coordinates.csv")
coordinate_list = {v[0]: v[1] for v in coordinate_list_data[["coord_id","coord_name"]].values}



##################################################################################################################
# Helper functions
##################################################################################################################
# def parse_genes(genes_textbox):
#
#     '''
#     sanitizes the list of genes provided in the textbox
#     converting it into a python array of valid ensembl genes
#     '''
#
#     # split at commas and remove trailing spaces
#     genes = [ g.strip() for g in genes_textbox.split(',') ]
#
#     # replace gene names with their ensembl ids
#     genes = [ genes_dropdown_options[g] if g in genes_dropdown_options else g for g in genes ]
#
#     # removes mistyped and duplicate genes
#     genes = list(set([ g for g in genes if g in genes_dictionary ]))
#
#     return genes
###################################
def convert_genenames_to_ensg(genes):
    # replace gene names with their ensembl ids
    #genes = [ genes_dropdown_options[g] if g in genes_dropdown_options else g for g in genes ]
    genes = [ genes_dropdown_options_inv[g] if g in genes_dropdown_options_inv else g for g in genes ]
    # removes mistyped and duplicate genes
    #genes = list(set([ g for g in genes if g in genes_dictionary ]))
    genes = list(set([ g for g in genes if g in genes_dict ]))
    return genes

###################################
def convert_ensg_to_genenames(genes):
    genes = [ genes_dropdown_options[g] for g in genes ]
    return genes

###################################
###################################
def parse_genes(genes_textbox):
    '''
    sanitizes the list of genes provided in the textbox
    converting it into a python array of valid ensembl genes
    '''

    # split at commas and remove trailing spaces
    genes = [ g.strip() for g in genes_textbox.split(',') ]
    #print("type"+type(genes))
    return convert_genenames_to_ensg(genes)






##################################################################################################################
##### The main window layout
##################################################################################################################

print("\n========= starting server =========\n")

server = flask.Flask(__name__)
app = dash.Dash(
  __name__,
  server=server,
  routes_pathname_prefix='/ccmut/')

#app.config.suppress_callback_exceptions = True ###Note, dangerous -- see if needed. used if components are added after initialization
app.title = "Data viewer: Changchun mutant dataset"
app.layout = html.Div([
    html.Div([
        html.Div([

            # textbox for selecting genes using ensembl id or names; all plots updates are connected to this element
            html.Label('Selected gene:'),
            dcc.Input(
                id='gene-textbox',
                type='text',
                value='',
                #list='list-suggested-inputs',  #Don't suggest EnsemblIDs
                placeholder='Comma-separated list of genes to inspect',
                style={'width': '100%', 'height': '40px'}
            ),
            # gene selection through dropdown; this will add the gene id to the textbox above
            html.Div([
                    dcc.Dropdown(
                    id='genes-dropdown',
                    value ='',
                    options=dropdown_genesymbols_sorted,
                    placeholder='Select a gene using its name',)
            ], id='genes-dropdown-timestamp', n_clicks_timestamp = 1),

            # strain selection through dropdown
            html.Div([
                    dcc.Dropdown(
                    id='strain-dropdown',
                    value ='',
                    options=dropdown_strains_sorted,
                    placeholder='Select a strain using its name',)
            ], id='strains-dropdown-timestamp', n_clicks_timestamp = 1),


            html.Div([html.Label(['X-axis:'])], style = {'display': 'block', 'width': '24%','height': '32px'} ),
            dcc.Dropdown(
                id='xaxis-selected',
                placeholder = 'coordinate',
                options=[{'label': coordinate_list[i], 'value': i} for i in coordinate_list],
                value='fracdel'),

            html.Div([html.Label(['Y-axis:'])], style = {'display': 'block', 'width': '24%','height': '32px'} ),
            dcc.Dropdown(
                id='yaxis-selected',
                placeholder = 'coordinate',
                options=[{'label': coordinate_list[i], 'value': i} for i in coordinate_list],
                value='fracfine'),

            html.Div([
                html.A([
                    html.Img(src=app.get_asset_url('MIMS_logo_blue.svg'), style={
                           'height': '30px',
                           #'float':'right',
                           'padding': '10px 10px'}),
                ], href='http://www.mims.umu.se/')
            ], style={'text-align':'right'})

        ], style={
            'border': 'thin lightgrey solid',
            'backgroundColor': 'rgb(250, 250, 250)',
            'padding': '10px 5px',
            'width':'100%'
        }),


        #####################################################
        ########## Gene information panel ###################
        #####################################################
        html.Div([

            #html.H1(html.Label(id='geneinfo-symbol')),
            html.H4(html.Label(id='geneinfo-longname')),
            html.H4(html.Label(id='straininfo-longname')),
            html.Label(id="geneinfo-alignment"),


            html.Div([
                html.A(['Ensembl'],    id='ensembl-link',   href='', target='_blank')," | ",
                html.A(['UniProt'],    id='uniprot-link',   href='', target='_blank')," | ",
                html.A(['PubMed'],     id='pubmed-link',    href='', target='_blank')," | ",
                html.A(['Genecards'],  id='genecards-link', href='', target='_blank')
            ])

        ], style={
            'margin-top':'50px',
            'border': 'thin lightgrey solid',
            'backgroundColor': 'rgb(250, 250, 250)',
            'padding': '10px 5px',
            'width':'100%'
        }, id='geneinfo-div'),
    ], style={'float':'right','width':'25%', 'padding':'20px'}),
    html.Div([
        dcc.Graph( id='scatter-plot')
    ],style={
      'display': 'inline-block',
      'margin': '0 auto'
    }),
    html.Div([
            html.H4("gene info"),
            html.Label(id="click-data",style={'font-family':'courier','white-space':'pre'}),
    ],style={
      'border':'1px solid black',
      'width':'100%',
      'height':'400px',
      'overflow-y':'scroll',
      'overflow-x':'scroll',
      'display': 'inline-block',
      'margin': '0 auto',
      'white-spacë́':'pre'
    }),
],style={
  'position': 'inline-block',
  'width': '95%',
  'height': '95%',
  'margin': '0 auto', #supposed to center it... but meah.
  'padding':'0',
  'overflow':'hidden'
})


##################################################################################################################
# page callbacks
# test genes = ENSMUSG00000000126,ENSMUSG00000095309,ENSMUSG00000097090
##################################################################################################################




################################################################################
## the genes textbox can be updated manually or using the dropdown
@app.callback(
    Output('gene-textbox', 'value'),
    [Input('genes-dropdown', 'value')])
def update_genes_dropdown(dropdown_value):
    if dropdown_value is None:
        return ''
    else:
        return dropdown_value




##################################################################################################################
# Function: Make the scatter plot for all the genes
##################################################################################################################
@app.callback(
    Output('scatter-plot', 'figure'),
    [Input('gene-textbox',    'value'),
     Input('strain-dropdown', 'value'),
     Input('xaxis-selected', 'value'),
     Input('yaxis-selected', 'value')])
def update_graph(selected_genes, selected_strains, xcoord_name, ycoord_name):

    print("get data")
    ##Load coordinates
    conn = sqlite3.connect("file:data/alignments.sqlite?mode=ro", uri=True)
    coord_data_plot = pd.read_sql_query("SELECT strain, wbid, `"+xcoord_name+"` as x, `"+ycoord_name+"` as y from alignment_sub", conn)
    conn.close()

    ##Check if any genes are selected and should be filtered
    if selected_genes is None or not len(selected_genes):
        selected_genes = []
    else:
        selected_genes = list(selected_genes.split(","))
        #selected_genes = convert_ensg_to_genenames(selected_genes)
    if len(selected_genes)>0:
        print("subset genes "+str(selected_genes))
        coord_data_plot = coord_data_plot[coord_data_plot["wbid"].isin(selected_genes)]

    ##Check if any strains are selected and should be filtered
    if selected_strains is None or not len(selected_strains):
        selected_strains = []
    else:
        selected_strains = list(selected_strains.split(","))
    if len(selected_strains)>0:
        print("subset strains "+str(selected_strains))
        coord_data_plot = coord_data_plot[coord_data_plot["strain"].isin(selected_strains)]

    ##Merge gene symbol
    coord_data_plot = coord_data_plot.merge(data_genelist,left_on = "wbid", right_on = "wbid") #, indicator = True)

    ##Check if there is any data to plot, otherwise return an empty graph
    if coord_data_plot.shape[0]==0:
        fig = go.Figure()
        fig.update_layout( autosize= False, width = 800, height = 800)
        return fig

    print(coord_data_plot.shape)
    print(coord_data_plot.head())

    ##Extract X,Y coordinates of all genes, their feature value (color), and the name
    xaxis = coord_data_plot["x"].values.tolist()
    yaxis = coord_data_plot["y"].values.tolist()
    textvalues = [x[0]+":"+x[1]+":"+x[2] for
	x in zip(
		coord_data_plot["strain"].values.tolist(),
		coord_data_plot["wbid"].values.tolist(),
		coord_data_plot["symbol"].values.tolist())]

    #Create the basic plot
    sgl= go.Scattergl(
            x = xaxis,
            y = yaxis,
            mode = "markers",
            #marker_color = markercolor,   TODO coler by PCR status... separate table? dataframe?
            text = textvalues,
            opacity = 1.0)

    fig = go.Figure(sgl)

    fig.update_layout( autosize= False, width = 800, height = 800, margin={'t':0, 'b':0,'l':0, 'r':0}, clickmode='event+select')
    return fig




@app.callback(
    [Output('click-data', 'children')],
    [Input('scatter-plot', 'clickData')])
def graph_click(clickData):

    #Extract which point
    print(str(json.dumps(clickData, indent=2)))
    if not clickData is None:
        ctext=clickData["points"][0]["text"]
        #print(str(json.dumps(clickData["points"][0], indent=2)))
        vals=ctext.split(":",3)
        strain=vals[0]
        wbid=vals[1]


        conn = sqlite3.connect("file:data/alignments.sqlite?mode=ro", uri=True)
        sqldata = pd.read_sql_query("SELECT strain, wbid, alignment from alignment where strain=? and wbid=?", conn,params=[strain, wbid])
        conn.close()

        al=sqldata["alignment"].tolist()[0]
        print(al)
        al.replace("\n","<br/>")

        print(strain+"   "+wbid)
        print(sqldata)
        return [al]

    return ["foooo"]

##################################################################################################################
##### Callback: Update gene information box ... links
##################################################################################################################
@app.callback(
    [
    Output('geneinfo-div',        'style'),
    #Output('geneinfo-symbol',       'children'),
    Output('geneinfo-longname',   'children'),
    Output('straininfo-longname',     'children'),
    Output('geneinfo-alignment',  'children'),

    Output('ensembl-link',   'href'),
    Output('uniprot-link',   'href'),
    Output('pubmed-link',    'href'),
    Output('genecards-link', 'href')
    ],
    [
      Input('gene-textbox',    'value'),
      Input('strain-dropdown', 'value')
    ])
def update_gene_links(gene,strain):

    ### hide genes links if more than a gene is provided
    selected_genes = parse_genes(gene)
    gene_id = selected_genes
    gene_symbol = convert_ensg_to_genenames(selected_genes)
    geneinfo_longname=""
    straininfo_longname=strain
    geneinfo_alignment=""

    if len(selected_genes) == 1:
        gene_symbol = gene_symbol[0]
        gene_id = gene_id[0]

        ### Make geneinfo box visible
        style= {
            'margin-top':'50px',
            'border': 'thin lightgrey solid',
            'backgroundColor': 'rgb(250, 250, 250)',
            'padding': '10px 5px',
            'width':'100%'
        }
    else:
        gene_symbol=""
        style= {'display': 'none', 'padding':'30px 0 0 30px'}


    ### update gene links
    gene = gene.strip()
    return (
        style,
        geneinfo_longname,
        straininfo_longname,
        geneinfo_alignment,

        'https://www.ensembl.org/Mus_musculus/Gene/Summary?g={}'.format(gene),
        'https://www.uniprot.org/uniprot/?query={}&sort=score'.format(gene),
        'https://pubmed.ncbi.nlm.nih.gov/?term={}'.format(gene_symbol),
        'https://www.genecards.org/cgi-bin/carddisp.pl?gene={}&keywords={}'.format(gene_symbol,gene_symbol)
    )



###################################
@app.callback(Output('genes-dropdown', 'value'),
    [Input('scatter-plot', 'clickData')])
def update_gene_textbox(value):
    # print(value)
    if value is None:
        return ''
    else:
        # print(value)
        pointed_gene = (value['points'])[0]['text']
        # print(pointed_gene)
        # print(str(convert_genenames_to_ensg([pointed_gene])).strip("['']"))
        return str(convert_genenames_to_ensg([pointed_gene])).strip("['']")



######################################
######################################

# run the app on "python app.py";
# default port: 8050
if __name__ == '__main__':
    app.run_server(debug = True)

app = dash.Dash(__name__)
#viewer.show(app)
