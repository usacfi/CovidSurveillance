'''
Date created: Sept 03, 2021
Created by: Jonathan Rico 

Example usage:
>>> import SarsCov2Variants as scv
>>> help(scv)
'''

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from datetime import datetime, date, timedelta
import shutil
from pathlib import Path

import warnings
warnings.filterwarnings('ignore')






def aggregate_divisions(metadata_df, country, division_column='division'):
    '''
    ======================================================================================
    
    Description: Aggregates the metadata divisions. The limitation of this function is that 
    the dictionary of alternative names of divisions varies from country to country.
    
    Returns: DataFrame with new column 'agg_division' 
    
    ======================================================================================
    
    <INPUTS>
    metadata_df: (dataframe)
    division_column: (str) column name 
    country: (str)
    
    ======================================================================================
    
    Example usage:
    >>> import SarsCov2Variants as scv
    >>> new_df, divisions = scv.aggregate_divisions(my_df, division_column='division')
    >>> new_df['agg_division'].value_counts()
    
    ======================================================================================
    '''
    
      
    # Enumerate all possible names of the divisions
    # You may check this with the command below:
    # >>> df[division_column].value_counts()
    philippines_dict = { 
        'NCR'         : ['ncr', 'manila', 'nationalcapitalregion','manilacity','navotascity','pasaycity','quezoncity','taguigcity','caloocancity'], 
        'CAR'         : ['car', 'cordilleraadministrativeregion'],
        'Region I'    : ['regioni', 'ilocosregion','ilocos'],
        'Region II'   : ['regionii', 'cagayanvalley'],
        'Region III'  : ['regioniii', 'centralluzon'],
        'Region IV-A' : ['regioniv-a', 'calabarzon'],
        'Region IV-B' : ['regioniv-b', 'mimaropa'],
        'Region V'    : ['regionv', 'bicolregion','bicol'],
        'Region VI'   : ['regionvi', 'westernvisayas'],
        'Region VII'  : ['regionvii', 'centralvisayas'],
        'Region VIII' : ['regionviii', 'easternvisayas'],
        'Region IX'   : ['regionix', 'zamboangapeninsula','zamboanga'],
        'Region X'    : ['regionx', 'northernmindanao'],
        'Region XI'   : ['regionxi', 'davaoregion','davao'],
        'Region XII'  : ['regionxii', 'soccsksargen'],
        'Caraga'      : ['caraga', 'davaooriental'],
        'BARMM'       : ['barmm', 'armm','bangsamoroautonomousregioninmuslimmindanao'],
        }
                        
                        
    indonesia_dict = {
        'Java'	                    : ['jakarta', 'banten', 'scrofjakarta', 'westjava', 'centraljava', 'srofyogyakarta', 'specialregionofyogyakarta', 'eastjava', 'jawabarat', 'jawatengah', 'jawatimur'],
        'Kalimantan'	            : ['westkalimantan', 'centralkalimantan', 'northkalimantan', 'eastkalimantan', 'southkalimantan', 'kalimantanselatan', 'kalimantantimur', 'kalimantanbarat', 'kalimantantengah'],
        'Maluku Islands'	        : ['northmaluku', 'maluku', 'malukuutara'],
        'Lesser Sunda Islands'	    : ['bali', 'westnusatenggara', 'eastnusatenggara', 'nusatenggaratimur', 'nusatenggarabarat'],
        'Western New Guinea'	    : ['westpapua', 'papua'],
        'Sulawesi'	                : ['northsulawesi', 'gorontalo', 'centralsulawesi', 'westsulawesi', 'southsulawesi', 'southeastsulawesi', 'sulawesitenggara'],
        'Sumatra'	                : ['aceh', 'northsumatra', 'westsumatra', 'riau', 'riauislands', 'jambi', 'bengkulu', 'southsumatra', 'bangkabelitungislands', 'lampung', 'sumatrabarat', 'sumatrautara', 'northsumatera', 'sumaterautara', 'sumaterabarat', 'kepulauanbangkabelitung', 'kepulauanriau'],
        }
                        
                        
    malaysia_dict = {               
        'Central Region'    : ['kualalumpur', 'selangor', 'negerisembilan', 'wilayahpersekutuan', 'dengkil', 'bukitbintang', 'rembau', 'salak'],                       
        'Sarawak'           : ['sarawak'],  
        'Sabah'             : ['sabah'],     
        'Northern Region'   : ['penang', 'perak', 'perlis'], 
        'East Coast'        : ['pahang', 'kelantan'], 
        'Southern Region'   : ['melaka', 'johor'],                             
        }                                   
                                                                         
    brunei_dict =   {'Brunei'    : ['brunei', 'sukang', 'nan']}
                        
    cambodia_dict = {'Phnom Penh'    : ['cambodia', 'phnompenh', 'sihanoukville', 'nan', 'steungtreng']}
                    
    laos_dict = {
        'Northern Region'   : ['vientianecapital', 'vientianeprovince'],
        'Central Region'    : ['savannakhet'],
        'Southern Region'   : ['saravan', 'champasack', 'luangprabang']
        }
    
    # we only considered the 7 states in myanmar and capital                   
    myanmar_dict = { 
        'Nay Pyi Taw': ['naypyitaw', 'yangon', 'mandalay'],
        'Rakhine'   : ['rakhine',], 
        'Kachin'    : ['kachin', 'kalay'],
        'Mon'       : ['mon', 'myeik'], 
        'Shan'      : ['shan', 'tamu'], 
        'Kayin'     : ['kayin'], 
        'Kayah'     : ['kayah'], 
        'Chin'      : ['chin'], 
        
        
        }                    
             
    singapore_dict = {'Singapore' : ['singapore', 'nan']}         
                        
    thailand_dict = {
        'North'     : ['tak', 'chiangrai', 'lamphun', 'phitsanulok', 'phayao', 'chiangmai', 'lampang', 'maehongson', 'nan', 'phayao', 'phrae', 'uttaradit'],
        'Northeast' : ['amnatcharoen', 'buengkan', 'buriram', 'chaiyaphum', 'kalasin', 'khonkaen', 'khonkean', 'loei', 'mahasarakham', 'mukdahan', 'nakhonphanom', 'nakhonratchasima', 'nongbualamphu', 'nongkhai', 'roiet', 'sakonnakhon', 'sisaket', 'surin', 'ubonratchathani', 'udonthani', 'yasothon'],
        'Central'   : ['bangkok', 'pathumthani', 'phatumthani', 'samutsakhon', 'nonthaburi', 'samutprakan', 'ratchaburi', 'angthong', 'lopburi', 'sukhothai', 'phitsanulok', 'phichit', 'kamphaengphet', 'phetchabun', 'nakhonsawan', 'uthaithani', 'nakhonnayok', 'phranakhonsiayutthaya', 'chainat', 'nakhonpathom', 'samutsongkhram', 'saraburi', 'singburi', 'suphanburi', 'kanchanaburi'],
        'East'      : ['chanthaburi', 'trat', 'rayong', 'chonburi', 'chonuri', 'chachoengsao', 'prachinburi', 'sakaeo'],
        'South'     : ['trang', 'suratthani', 'songkhla', 'narathiwat', 'phuket', 'pattani', 'krabi', 'phangnga', 'ranong', 'satun', 'chumphon', 'nakhonsithammarat', 'phatthalung', 'yala', 'phetchaburi', 'prachuapkhirikhan'],
        }
    
    timor_leste_dict = {'Timor-Leste'   : ['easttimor', 'nan']}
                        
    vietnam_dict = {
        'Northwest'             : ['yenbai', 'laocai', 'sonla', 'hoabinh'],
        'Northeast'             : ['thainguyen', 'bacgiang', 'quangninh', 'langson', 'caobang'],
        'Red River Delta'       : ['hanoi', 'vinhphuc', 'haiduong', 'hungyen', 'bacninh', 'haiphong', 'hanam', 'thaibinh', 'namdinh', 'phutho'],
        'North Central Coast'   : ['nghean', 'thanhhoa', 'hatinh'],
        'South Central Coast'   : ['southcentralcoast', 'danang'],
        'Central Highlands'     : ['centralhighlands'],
        'Southeast'             : ['south-easternregion', 'hochiminhcity'],
        'Mekong River Delta'    : ['kiengiang'],
        }
    
    if country.lower().strip() == 'philippines':
        dictionary = philippines_dict
    elif country.lower().strip() == 'indonesia':
        dictionary = indonesia_dict
    elif country.lower().strip() == 'malaysia':
        dictionary = malaysia_dict
    elif country.lower().strip() == 'brunei':
        dictionary = brunei_dict
    elif country.lower().strip() == 'cambodia':
        dictionary = cambodia_dict
    elif country.lower().strip() == 'laos':
        dictionary = laos_dict
    elif country.lower().strip() == 'myanmar':
        dictionary = myanmar_dict
    elif country.lower().strip() == 'singapore':
        dictionary = singapore_dict
    elif country.lower().strip() == 'thailand':
        dictionary = thailand_dict
    elif country.lower().strip() == 'timor-leste':
        dictionary = timor_leste_dict
    elif country.lower().strip() == 'vietnam':
        dictionary = vietnam_dict
    elif country.lower().strip() == 'southeastasia':
        pass
    else:
        raise ValueError('\nError: Unknown country.\n Choose from philippines, indonesia, malaysia, brunei, cambodia,' 
        'laos, myanmar, singapore, thailand, timor-leste, vietnam.')
        
    
    df = metadata_df.copy()
    
    try:
        df['agg_division'] = 'Unknown'
        for division, alt_names in dictionary.items():
            for alt_name in alt_names:
                indices = df[df[division_column].astype(str).str.lower().str.strip().str.replace(" ","") == alt_name].index
                df.loc[indices, 'agg_division'] = division
    except:
        df['agg_division'] = df[division_column]
            
    return df


    























































def area_charts_of_divisions(metadata_df, country, normalized=False, date_column='date', 
                            division_column='agg_division', id_column='gisaid_epi_isl', 
                            variant_column='variant', output_directory='output/11_regions'):
    '''
    ======================================================================================
                                                    
    Description: Plots the area chart of each geographical division
    
    Creates area charts (an area chart per geographical division) in png format 
            and saves to the output_directory.
    ======================================================================================
    
    <INPUTS>
    metadata_df: (dataframe)
    normalized: (bool) If True, the area charts are normalized per aggregation mode (monthly),
                            Default: False
    date_column: (str) column name, Default: 'date'
    division_column: (str) column name, Default: 'agg_division'
    id_column: (str) column name, Default: 'gisaid_epi_isl'    
    variant_column: (str) column name, Default:'variant'
    output_directory: (str) destination folder of the output png files.
                            
    ======================================================================================                                                                                               
    '''
    
    metadata_df[date_column] = pd.to_datetime(metadata_df[date_column])
    metadata_df[date_column] = metadata_df[date_column].dt.to_period('M').dt.start_time
    
    groupby = metadata_df[[ variant_column,
                            date_column,
                            division_column,
                            id_column]].groupby(by=[date_column,
                                                    variant_column,
                                                    division_column], as_index=False).count()  
    plt.style.use('seaborn')  
    countmax = groupby.groupby(by=[date_column,division_column])[id_column].sum().max()                                           
    
    divisions = metadata_df[division_column].unique().tolist()
    divisions = [col for col in divisions if col not in ('Unknown','?','nan')]                                                            
    for i in range(len(divisions)):
        region = divisions[i]
        reg = groupby[groupby[division_column]==region] 
    
        # If there's no sufficient data  
        if len(reg) < 4:  
            df = pd.DataFrame({date_column:['2020-10','2021-01','2021-03'],'Insufficient data':[0,0,0]})
            #df[date_column] = df[date_column].astype('datetime64')
            #df[date_column] = df[date_column].dt.date
            COLORS = 'black'
   
        # If there's sufficient data
        else:
            df = reg.pivot(index=date_column, columns=variant_column)[id_column]
            df = df.replace(np.nan, 0)
            
            # Put the 'Others' column at the end, assumes 'Others' is the last alphabetically
            cols = df.columns.tolist()
            df = df[cols]
            
            # Normalize the occurences
            if normalized:
                countmax = 1
                df = (df.T/df.sum(axis=1)).T
      
            new_index = np.arange((df.index[-1] - df.index[0]).days + 1)
            date_index = pd.date_range( start = df.index[0], 
                                        end = df.index[-1],
                                        freq='D').date

            if normalized:                            
                # Smoothen curves                            
                temp_df = pd.DataFrame()
                for variant in df.columns:
                    func = interp1d((df.index - df.index[0]).days, 
                                    df[variant], 
                                    kind='cubic')
                    temp_df[variant] = func(new_index)

                temp_df[temp_df < 0] = 0     
                temp_df.index = date_index
                df = temp_df
                df = (df.T/df.sum(axis=1)).T  # Ensure the sum equals 1
                
            else:
                temp_df = pd.DataFrame(index=date_index)
                df = df.combine(temp_df, np.maximum, fill_value=0)
                    
                
            variant_color_df = pd.DataFrame({variant_column:df.columns.tolist()})
            variant_color_df = variant_color(variant_color_df)
                
            COLORS = variant_color_df.color.tolist()

        df.plot( kind = 'area', 
                 rot = 0,
                 color = COLORS)

                       
        ax = plt.gca()
        handles, labels = ax.get_legend_handles_labels()
        labels, handles = zip(*sorted(zip(labels, handles), 
                                      key=lambda t: t[0], 
                                      reverse=True))
        ax.legend(handles, labels)
        plt.ylim(0, countmax)
        plt.xlim(metadata_df[date_column].min().date(), metadata_df[date_column].max().date())
        plt.xlabel(date_column)
        year = str(datetime.today().year) 
        month = str(datetime.today().month).zfill(2)
        day = str(datetime.today().day).zfill(2)
        update = year + month + day 
        n = ''
        if normalized:
            plt.ylabel('Count')
            update = f'N_{update}' 
            n = '_N'  
        plt.legend(loc = 'upper left') 
        plt.title('{}'.format(region), 
                        fontsize = 17, 
                        fontstyle = 'oblique') 
                        
        plt.tight_layout()
        Path(f'{output_directory}/{country.lower().strip()}').mkdir(exist_ok=True, parents=True)
        tableau_df = convert_to_tableau_df(df, column_names=['count', 'variant'])
        tableau_df.to_csv(f'{output_directory}/{country}/{country}_{region}{n}.csv', index=True)
        #plt.savefig(f'{output_directory}/{country}/{i}_{region}_{update}.png')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        

def bubble_map_of_country(metadata_df, date_column='date', division_column='agg_division',
                        variant_column='variant', id_column='gisaid_epi_isl', color_column='color',
                        latitude_column='div_latitude', longitude_column='div_longitude',
                        output_filename='output/10_geoplot_variants.html'):
    '''
    ======================================================================================
                        
    Description: Creates a bar chart and a geographical bubble map of the variants per region
    
    ======================================================================================                    
    
    <INPUTS>
    
    metadata_df: (dataframe)
    date_column: (str) column name
    division_column: (str) column name
    variant_column: (str) column name
    id_column: (str) column name
    color_column: (str) column name
    latitude_column: (str) column name
    longitude_column: (str) column name
    output_filename: (str) directory and filename of the output file
    ======================================================================================
    
    Example usage:
    >>> import SarsCov2Variants as scv
    >>> scv.bubble_map_of_country(my_df)      
    
    ======================================================================================                                  
    '''
                        
    meta_df = metadata_df.copy()
      
    # Calculate the count per variant per division
    # | variant | division | case_count | 
    # |'string' |  'str'   |   integer  | 
    casecount_df = meta_df[[ variant_column,
                             division_column,
                             latitude_column,
                             longitude_column,
                             id_column]].groupby([ variant_column,
                                                   division_column,
                                                   latitude_column,
                                                   longitude_column,
                                                  ], as_index=False).count()
                                            
    variants = meta_df[variant_column].unique().tolist()
    colors = meta_df[color_column].unique().tolist()
    color_dict = dict(zip(variants, colors))
    
    fig = make_subplots(rows=1, cols=2, specs=[[{'type':'scattergeo'}, {'type':'xy'}]])

    # Bar Chart 
    # reorder the variants for better view of overlapping bubbles
    variants.sort()
    for variant in variants:
        divisions = casecount_df.query(f'{variant_column} == "{variant}"')[division_column].tolist()
        casecounts = casecount_df.query(f'{variant_column} == "{variant}"')[id_column].tolist()
        
        # Create a temporary dataframe so we can sort the divisions accordingly
        temp_df = pd.DataFrame({'division':divisions,
                                'casecount':casecounts
                                })
        temp_df.sort_values(by='casecount', inplace=True, ascending=False)
        
        fig.add_trace(go.Bar( 
            name = variant, 
            x = temp_df.division, 
            y = temp_df.casecount, 
            legendgroup = f'group_{variant}',
            marker_color = color_dict[variant]),
            row=1, col=2)      



    # Scatter Map
    variants.reverse()
    for variant in variants:
        divisions = casecount_df.query(f'{variant_column} == "{variant}"')[division_column].tolist()
        longitudes = casecount_df.query(f'{variant_column} == "{variant}"')[longitude_column].tolist()
        latitudes = casecount_df.query(f'{variant_column} == "{variant}"')[latitude_column].tolist()
        casecounts = casecount_df.query(f'{variant_column} == "{variant}"')[id_column].tolist()
        
        # Create a temporary dataframe so we can sort the divisions accordingly
        temp_df = pd.DataFrame({'division':divisions,
                                'longitude':longitudes,
                                'latitude':latitudes,
                                'casecount':casecounts
                                })
        temp_df.sort_values(by='casecount', inplace=True, ascending=False)                        
                                
        
        fig.add_trace(go.Scattergeo(
            lon = temp_df.longitude,
            lat = temp_df.latitude,
            text = temp_df.division,
            name = meta_df['country'].iloc[0], 
            showlegend = False,
            customdata = temp_df.casecount,  
            hovertemplate = 'Area: %{text}<br>'
                           f'Variant: {variant}<br>'
                            'Case count: %{customdata}',
            legendgroup = f'group_{variant}',
            marker = dict(
                            size = temp_df.casecount / 10 + 2,  
                            color = color_dict[variant], 
                            line_width = 0,
                            opacity  = 0.5,
                        )),
            row=1, col=1)


    if country.lower().lstrip().rstrip() == 'philippines':
        geo_lonaxis_range = [ 117.0, 127.0 ]
        geo_lataxis_range = [ 5.0, 20.0 ]
    elif country.lower().lstrip().rstrip() == 'indonesia':
        geo_lonaxis_range = [ 94.0, 142.0 ]
        geo_lataxis_range = [ -11.0, 6.50 ]
    else:
        print('\nError: Unknown country.\n')
        exit()

    # Layout design
    fig.update_layout(title = 'Sars-Cov-2 Variants per Region | {} | {} to {} | {} GISAID Sequences'.format(meta_df['country'].iloc[0],
                            meta_df['date'].value_counts().index[-1].date(),
                            meta_df['date'].value_counts().index[0].date(),
                            len(meta_df)),
                      barmode = 'stack',
                      geo_resolution = 50,
                      geo_scope = meta_df.region.iloc[0].lower(),
                      geo_showframe = True,
                      geo_showcoastlines = True,
                      geo_landcolor = "rgb(217, 217, 217)",
                      geo_countrycolor = "black" ,
                      geo_coastlinecolor = "black",
                      geo_projection_type = 'mercator',
                      geo_lonaxis_range = geo_lonaxis_range, 
                      geo_lataxis_range = geo_lataxis_range,
                      legend = dict(
                            yanchor = "top",
                            y = 0.99,
                            xanchor = "left",
                            x = 0.89
                            ),
                      legend_title_text = 'Sars-Cov-2 Variants',
                      )
                            
    fig['layout']['yaxis']['title'] = 'Count'
    fig.write_html(f'{output_filename}')

    























  
    
def calc_metrics(df_, country, country_column='country', submission_date_column='date_submitted', 
    collection_date_column='date', period=30):
    
    cols = [submission_date_column, collection_date_column, country_column]
    df = df_[cols].copy()
    df = df[(~df[country_column].isnull()) & (~df[submission_date_column].isnull()) & (~df[collection_date_column].isnull())]
    df = df[df[country_column]==country]

    # Convert index to datetime dtype
    df[submission_date_column] = pd.to_datetime(df[submission_date_column])
    df[collection_date_column] = pd.to_datetime(df[collection_date_column])

    # Aggregate per submission date
    agg_df = df.groupby(submission_date_column)[country_column].count()
    
    start_date = agg_df.index[0]
    current_date = datetime.now() #date(2022, 4, 17)
    delta = current_date - start_date
        
    date_range = []    
    submission_sum = []    
    submission_freq = []
    submission_cst = []
    for i in range(0, delta.days+1):
        d = start_date + timedelta(days=i) 
        date_range.append(d)
        
        # Filter the dates    
        temp_df = agg_df[(agg_df.index > d - timedelta(days=period)) & (agg_df.index <= d)]
        
        # Calculate the running total (in the past 30 days) of submissions
        submission_sum.append(temp_df.sum())
        
        # Calculate the moving average of submission frequency
        l1 = temp_df.index[:-1].to_frame().reset_index(drop=True)
        l2 = temp_df.index[1:].to_frame().reset_index(drop=True)
        if len(temp_df) < 2:
            submission_freq.append(np.nan)
        else:
            submission_freq.append(np.abs(np.mean(l1-l2).astype('timedelta64[D]').astype(int))[0])
        
        # Calculate average difference between submission and collection dates
        temp2_df = df[(df[submission_date_column] > d - timedelta(days=period)) & (df[submission_date_column] <= d)]
        diff = (temp2_df[submission_date_column] - temp2_df[collection_date_column]).astype('timedelta64[D]')
        submission_cst.append(np.mean(diff))      
    
    # Create output dataframe
    out_df = pd.DataFrame({'date':date_range,
                           'sum':submission_sum,
                           'freq':submission_freq,
                           'cst':submission_cst,
                         })
                         
    out_df['country'] = country
    return out_df





























def convert_to_tableau_df(df, column_names=['values', 'color']):
    '''
    ======================================================================================
                    
    Description: Merges all columns of the dataframe into one column then adds a second column
                whose value is the corresponding original column name.
    
    Returns a dataframe ready to input in Tableau for heatmap or area chart.
    
    ======================================================================================
    <INPUTS>
    df: (dataframe) pandas dataframe with several columns to be merged
    column_names: (list) list of the desired output column names (str) not the input column names
    
    ====================================================================================== 
       
    Example usage:
    >>> import SarsCov2Variants as scv                
    >>> scv.convert_to_tableau_df(my_df, column_names=['values', 'color'])   
    
    my_df:
      |  A  |  B  |  C  | 
    1 | 0.1 | 0.3 | 0.2 |
    2 | 0.5 | 0.4 | 0.7 |
    
    out_df:
      | values | color | 
    1 |  0.1   |   A   | 
    2 |  0.5   |   A   |
    1 |  0.3   |   B   |
    2 |  0.4   |   B   |
    1 |  0.2   |   C   |
    2 |  0.7   |   C   |
    
    ======================================================================================
    '''
    
    cols = df.columns.tolist()
    out_df = df[cols[0]].copy()

    l = [cols[0]]*len(df)
    for col in cols[1:]:
        out_df = pd.concat([out_df, df[col]])
        l.extend([col]*len(df))
        
    out_df = out_df.to_frame()
    out_df[column_names[1]] = l
    out_df.columns = column_names
    
    return out_df   





















































def combine_metadata(directory, meta_cols=['strain','gisaid_epi_isl','date','region','country',
                    'division','age','sex','pangolin_lineage','GISAID_clade','originating_lab',
                    'submitting_lab','date_submitted'], date_column='date', with_phrase='metadata',
                    without_phrase='aggregated'):
    '''
    ======================================================================================
                    
    Description: Searches and combines all 'metadata.tsv' files inside the directory into 
                a single dataframe. The dataframe is sorted by the collection date. 
    
    Returns the combined dataframe.
    
    ======================================================================================
    <INPUTS>
    directory: (str) where the folders containing the metadata.tsv files are located
    meta_cols: (list) list of column names (str)
    date_column: (str) column name by which to sort the combined dataframe. Default:'date'
    
    ====================================================================================== 
       
    Example usage:
    >>> import SarsCov2Variants as scv            
    >>> scv.combine_metadata('/User/Downloads', meta_cols=['virus','division'], date_column='date')   
    
    ====================================================================================== 
    '''
    
    def convert_to_date(df, date_columns):
        # Solve the issues with different date formats
        new_df = df.copy()
        for date_column in date_columns:
            new_df['temp_date'] = pd.to_datetime(new_df[date_column], format='%d/%m/%Y', errors='coerce')  
            mask = new_df['temp_date'].isnull()
            new_df.loc[mask, 'temp_date'] = pd.to_datetime(new_df[date_column], format='%m/%d/%Y', errors='coerce')
            mask = new_df['temp_date'].isnull()
            new_df.loc[mask, 'temp_date'] = pd.to_datetime(new_df[date_column], format='%Y-%m-%d', errors='coerce')    
            mask = new_df['temp_date'].isnull()
            new_df.loc[mask, 'temp_date'] = pd.to_datetime(new_df[date_column], format='%Y-%d-%m', errors='coerce')
            new_df[date_column] = new_df['temp_date'].copy()
            new_df = new_df.drop(columns='temp_date')
        return new_df
    
    
    metadata_df = pd.DataFrame()
    
    count = 0
    for dirname, _, filenames in os.walk(directory):
        for filename in filenames:
            if with_phrase in filename and without_phrase not in filename and 'southeastasia' not in filename:
                print(f'{dirname}/{filename}')
                if filename[-3:] == 'tsv':
                    sep = '\t'
                else:
                    sep = ','
                temp_df = pd.read_csv(f'{dirname}/{filename}', 
                                        sep=sep, 
                                        usecols=meta_cols)
                temp_df = convert_to_date(temp_df, [date_column, 'date_submitted'])
                metadata_df = pd.concat([metadata_df, temp_df], ignore_index=True)
                count = count + 1

    if len(metadata_df) == 0:
        raise ValueError('\nError: No data was found. Check the spelling of filenames and check if the files are not empty.\n')
    else:
        print(f'Found {count} files')
        print(f'Total GISAID sequences: {len(metadata_df)}\n')
        
    metadata_df.sort_values(by=date_column, ascending=False, inplace=True)
    metadata_df = metadata_df.reset_index(drop=True)
    
    return metadata_df























































def combine_output_files(path, use_filename=False, include=[], exclude=[], output_index=False, first_col=None):
    '''
    ======================================================================================
    
    Description: Combines output csv files in a directory
    
    Returns: Combined dataframe
    
    ======================================================================================
    
    <INPUTS>
    path: (directory)
    use_filename: (bool) whether to create another column in the dataframe with the filename 
                    as the value
    include: (list) filters which csv files to combine using a list of keywords
    exclude: (list) filters which csv files to ignore using a list of keywords
    output_index: (bool) whether to save the output csv with index or without index
    first_col: (str) column name for which to rename the first column
    
    ======================================================================================
    
    Example usage:
    >>> import SarsCov2Variants as scv
    >>> scv.combine_output_files('output/12_variants', use_filename=False, include=['aggregated'])
    
    ======================================================================================
    '''
    
    width = os.get_terminal_size().columns
    print('\n\n')
    print('#===================================#'.center(width))
    print('Combining output files'.center(width))
    print('#===================================#'.center(width))
    print('\n\n')
    
    combined_df = pd.DataFrame()
    l = []
    
    for dirname, _, filenames in os.walk(path):
        for filename in filenames:
            if len(include)>0 and len([key for key in include if key in filename]) == 0:
                print('skipped', filename)
                continue
            elif len(exclude)>0 and len([key for key in exclude if key in filename]) > 0:
                print('skipped', filename)
                continue
            
            print(os.path.join(dirname, filename))
            temp_df = pd.read_csv(os.path.join(dirname, filename))    
            combined_df = pd.concat([combined_df, temp_df])
        
            if use_filename:
                l.extend([filename] * len(temp_df))
          
    if first_col != None:
        combined_df = combined_df.rename(columns={combined_df.columns[0]:first_col})
        
    if use_filename:
        combined_df['filename'] = l

    combined_df.to_csv(f'{path}.csv', index=output_index)
    return combined_df























































def geocode_divisions(metadata_df, country, division_column='agg_division'):
    '''
    ======================================================================================
    
    Description: Maps the divisions to the dictionary of latitude and longitude (WSG84) coordinates. 
    
    Returns: DataFrame with two new columns 'div_latitude' and 'div_longitude'
    
    ======================================================================================
    
    <INPUTS>
    metadata_df: (dataframe)
    division_column: (str) column name Default: 'agg_division'
    country: (str)
    
    ======================================================================================
    
    Example usage:
    >>> import SarsCov2Variants as scv
    >>> new_df = scv.geocode_divisions(my_df, division_column='agg_division')
    
    ======================================================================================
    '''
    
    
    philippines_dict =  {
        'NCR'         : [14.5472655, 121],
        'CAR'         : [17.35971005, 121.07525292],
        'Region I'    : [15.855663, 120.1341247],
        'Region II'   : [16.7305959, 121.54067657],
        'Region III'  : [15.5, 121.13104248],
        'Region IV-A' : [14.16052895, 121.24256648],
        'Region IV-B' : [13.072377, 121.3276166],
        'Region V'    : [13.42289175, 123.41103829],
        'Region VI'   : [11, 122.58101501],
        'Region VII'  : [10.6419531, 123.938142],
        'Region VIII' : [11.9, 125],
        'Region IX'   : [8.5, 123.4243754],
        'Region X'    : [8.3, 124.65677618],
        'Region XI'   : [7.0428208, 125.580953],
        'Region XII'  : [6.2965755, 124.9860759],
        'Caraga'      : [9.2471392, 125.85578189],
        'BARMM'       : [5.2, 120.02840346],
        }
    
    indonesia_dict = {
        'Java'	                    : [-6.430323, 107.500769],
        'Sumatra'	                : [-1.094433, 101.977441],
        'Kalimantan'	            : [0.529112, 114.198348],
        'Sulawesi'	                : [-1.833530, 120.287056],
        'Lesser Sunda Islands'	    : [-9.788183, 120.026111],
        'Western New Guinea'	    : [-3.657968, 137.944309],
        'Maluku Islands'	        : [-3.180422, 129.115683],
        }
      
    malaysia_dict = {
        'Central Region'     : [3.128231, 101.839143],
        'Sarawak'            : [2.458302, 113.446666],
        'Sabah'              : [5.253195, 116.874400],
        'Northern Region'    : [5.310380, 101.260808],                           
        'East Coast'         : [3.808925, 103.029524],
        'Southern Region'    : [1.868505, 103.550540],
        }  
    
    brunei_dict =   {'Brunei' : [4.410266, 114.609550]}
    
    cambodia_dict = {'Phnom Penh' : [11.561475, 104.890590]}
    
    laos_dict = {
        'Northern Region'   : [20.076666, 102.596084],
        'Central Region'    : [17.216986, 105.309707],
        'Southern Region'   : [15.291582, 106.531822],
        }

    myanmar_dict = { 
        'Nay Pyi Taw': [20.581489, 96.053544],
        'Rakhine'   : [19.863024, 94.241567], 
        'Kachin'    : [25.509609, 96.818941],
        'Mon'       : [15.985101, 98.097306], 
        'Shan'      : [21.453758, 99.180514], 
        'Kayin'     : [17.805977, 96.727180], 
        'Kayah'     : [19.640226, 97.168350], 
        'Chin'      : [23.337089, 93.952457], 
        }
    
    singapore_dict = {'Singapore' : [1.352, 103.8]}
    
    thailand_dict = {
        'North'     : [18.531524, 99.468971],
        'Northeast' : [16.151894, 103.329354],
        'Central'   : [14.312884, 99.935434],
        'East'      : [13.172319, 101.954093],
        'South'     : [8.469181, 99.311848],
        }
    
    timor_leste_dict = {'Timor-Leste'   : [-8.847402, 125.860190]}
    
    vietnam_dict = {
        'Red River Delta'       : [20.795060, 105.940353],
        'Northwest'             : [21.471386, 103.589279],
        'Northeast'             : [22.286977, 105.918380],
        'North Central Coast'   : [18.122794, 105.830489],
        'Southeast'             : [11.391634, 107.302657],
        'South Central Coast'   : [15.663149, 108.027755],
        'Central Highlands'     : [13.280490, 108.005782],
        'Mekong River Delta'    : [9.944953, 105.544845],
        }
    
    southeastasia_dict = {
        'Brunei'                : [4.410266, 114.609550],
        'Cambodia'              : [11.561475, 104.890590],
        'Indonesia'             : [-6.430323, 107.500769],
        'Laos'                  : [20.076666, 102.596084],
        'Malaysia'              : [3.128231, 101.839143],
        'Myanmar'               : [20.581489, 96.053544],
        'Philippines'           : [14.5472655, 121],
        'Singapore'             : [1.352, 103.8],
        'Thailand'              : [18.531524, 99.468971],
        'Timor-Leste'           : [-8.847402, 125.860190],
        'Vietnam'               : [20.795060, 105.940353],
        }
    
    df = metadata_df.copy()
    df['div_latitude'] = np.nan
    df['div_longitude'] = np.nan
    
    if country.lower().strip() == 'philippines':
        dictionary = philippines_dict
    elif country.lower().strip() == 'indonesia':
        dictionary = indonesia_dict
    elif country.lower().strip() == 'malaysia':
        dictionary = malaysia_dict    
    elif country.lower().strip() == 'brunei':
        dictionary = brunei_dict
    elif country.lower().strip() == 'cambodia':
        dictionary = cambodia_dict
    elif country.lower().strip() == 'laos':
        dictionary = laos_dict
    elif country.lower().strip() == 'myanmar':
        dictionary = myanmar_dict
    elif country.lower().strip() == 'singapore':
        dictionary = singapore_dict
    elif country.lower().strip() == 'thailand':
        dictionary = thailand_dict
    elif country.lower().strip() == 'timor-leste':
        dictionary = timor_leste_dict
    elif country.lower().strip() == 'vietnam':
        dictionary = vietnam_dict
    elif country.lower().strip() == 'southeastasia':
        dictionary = southeastasia_dict
    else:
        raise ValueError('\nError: Unknown country.\n Choose from philippines, indonesia, malaysia, brunei, cambodia,' 
        'laos, myanmar, singapore, thailand, timor-leste, vietnam')
    
    for division, coordinates in dictionary.items():
        latitude = coordinates[0]
        longitude = coordinates[1]
        
        indices = df[df[division_column]==division].index
        df.loc[indices, 'div_latitude'] = latitude
        df.loc[indices, 'div_longitude'] = longitude
    
    return df





































def init_functions(directory, country, lineage_column='pangolin_lineage', 
            division_column='division', date_column='date', variant_column='variant', 
            id_column='gisaid_epi_isl', with_phrase='metadata', without_phrase='aggregated',
            meta_cols=['strain', 'gisaid_epi_isl', 'date', 'region', 
            'country', 'division', 'age', 'sex', 'pangolin_lineage', 'GISAID_clade', 
            'originating_lab', 'submitting_lab', 'date_submitted'],):
            
    '''
    ======================================================================================
    
    Description: Runs all the initial (preprocessing) functions in SarsCov2Variants.
            
    
    Returns a dataframe and saves it as a csv file.     
            
    ======================================================================================   
    '''
    
    output_directory = f'output/12_variants'
    Path(output_directory).mkdir(exist_ok=True, parents=True)
    
    meta_df = combine_metadata(directory, meta_cols, date_column, with_phrase, without_phrase)
    meta_df = lineage_to_variant(meta_df, lineage_column)
    meta_df = aggregate_divisions(meta_df, country, division_column)
    division_column='agg_division'
    meta_df = geocode_divisions(meta_df, country, division_column)
    latitude_column='div_latitude'
    longitude_column='div_longitude'
    meta_df = variant_color(meta_df, variant_column)
    color_column='color'
    
    casecount_df = meta_df[[ variant_column,
                             division_column,
                             latitude_column,
                             longitude_column,
                             id_column]].groupby([ variant_column,
                                                   division_column,
                                                   latitude_column,
                                                   longitude_column,
                                                  ], as_index=False).count() 
    
    casecount_df.columns = ['Variant', 
                            'Location', 
                            'Latitude', 
                            'Longitude',
                            'Case Count']
    
    last_collection_date = meta_df.sort_values(date_column)[date_column].dropna().iloc[-1]
    last_submission_date = meta_df.sort_values('date_submitted')['date_submitted'].dropna().iloc[-1] 
    
    casecount_df['Last Collection Date'] = last_collection_date
    casecount_df['Last Submission Date'] = last_submission_date                                             
    
    meta_df.to_csv(f'{output_directory}/12_variants_in_{country}.csv', index=False)
    casecount_df.to_csv(f'{output_directory}/12_variants_in_{country}_aggregated.csv', index=False)
    
    return meta_df










































def lineage_to_variant(metadata_df, lineage_column='pangolin_lineage'):
    '''
    ======================================================================================
    
    Description: Maps Pangolin lineages into WHO variant labels (check updates below):
    https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/
    
    Returns dataframe
    
    ======================================================================================    
    <INPUTS>
    metadata_df: (object) pandas dataframe
    lineage_column: (str) column name of the Pangolin lineage
    
    ======================================================================================
        
    Example usage:
    >>> import SarsCov2Variants as scv
    >>> scv.lineage_to_variant(my_df, lineage_column='pangolin_lineage')
    
    ======================================================================================
    '''
    
    # ALWAYS CHECK FOR UPDATES IN THE LINK ABOVE!
    variants_dict = {
                'Alpha'  : ['B.1.1.7', 'Q'], 
                'Beta'   : ['B.1.351'], 
                'Gamma'  : ['P.1'], 
                'Delta'  : ['B.1.617.2', 'AY'], 
                'Omicron': ['B.1.1.529', 'BA'],
                'Lambda' : ['C.37'],
                'Mu'     : ['B.1.621'],
                }
    
    df = metadata_df.copy()
    df['variant'] = 'Others'
    
    for variant, lineages in variants_dict.items():
        for lineage in lineages:
            indices = df[df[lineage_column].astype(str).str.startswith(lineage)].index
            df.loc[indices, 'variant'] = variant
            
    return df































def main(countries=['Brunei', 'Cambodia', 'Indonesia', 'Laos', 'Malaysia', 'Myanmar', 'Philippines',
        'Singapore', 'Thailand', 'Timor-Leste', 'Vietnam']):
        
        for country in countries:
                print(country)
                shutil.rmtree(f'output/11_regions/{country}', ignore_errors=True)
                meta_df = init_functions(f'references/Sequences/{country}/Variant_Surveillance', country=country)
                area_charts_of_divisions(meta_df, country, True)
                area_charts_of_divisions(meta_df, country)

        print('Southeastasia')
        southeast_asia('output/12_variants', countries)
        combine_output_files('output/11_regions', 
                            use_filename=True, 
                            exclude=['nan', '2021', 'DS'], 
                            output_index=False, 
                            first_col='Collection Date')
        combine_output_files('output/12_variants', 
                            use_filename=True, 
                            include=['aggregate'], 
                            output_index=False)
            































def southeast_asia(input_directory, countries):
    
    # Create southeast asia csv and southeast asia aggregated csv
    meta_df = init_functions(input_directory, 
                            country='southeastasia', 
                            with_phrase='12_variants_in', 
                            division_column='country')
                           
    # Create area charts for each southeast asian country 
    shutil.rmtree('output/11_regions/southeastasia', ignore_errors=True)                           
    area_charts_of_divisions(meta_df, country='southeastasia')
    area_charts_of_divisions(meta_df, country='southeastasia', normalized=True)
    
    # Create time-series of GISAID submission metrics
    for country in countries:
        print(country)
    
        # Calculate metrics
        temp_df = calc_metrics(meta_df, country=country)
        try:
            df = pd.concat([df, temp_df])
        except:
            df = temp_df.copy()
    
    df.to_csv('output/submission_metrics.csv', index=False)
    return df



































def variant_color(metadata_df, variant_column='variant'):
    '''
    ======================================================================================
    
    Description: Maps each variant to certain colors.

    Returns dataframe with new column 'colors'

    ======================================================================================    
    <INPUTS>
    metadata_df: (object) pandas dataframe
    lineage_column: (str) column name of the WHO variant, Default: 'variant'
    
    ======================================================================================
        
    Example usage:
    >>> import SarsCov2Variants as scv
    >>> scv.variant_color(my_df, variant_column='variant')
    
    ====================================================================================== 
    '''
    
    colors_dict = {
                    'Alpha' : 'maroon',
                    'Beta'  : 'royalblue',
                    'Gamma' : 'lightgreen',
                    'Delta' : 'gold',
                    'Eta'   : 'orange',
                    'Iota'  : 'green',
                    'Omicron': 'pink',
                    'Lambda': 'lightblue',
                    'Mu'    : 'magenta',
                    'Others': 'whitesmoke',        
                }
    
    df = metadata_df.copy()
    df['color'] = 'whitesmoke'
    
    for variant, color in colors_dict.items():
        indices = df[df[variant_column]==variant].index
        df.loc[indices, 'color'] = color
    
    return df














# ================#
#       MAIN      #
# ================#

if __name__ == "__main__":
    main()








