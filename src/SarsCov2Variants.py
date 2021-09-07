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
import warnings

warnings.filterwarnings('ignore')






def aggregate_divisions(metadata_df, division_column='division', country='philippines'):
    '''
    ======================================================================================
    
    Description: Aggregates the metadata divisions. The limitation of this function is that 
    the dictionary of alternative names of divisions varies from country to country.
    
    Returns: DataFrame with new column 'agg_division' 
    
    ======================================================================================
    
    <INPUTS>
    metadata_df: (dataframe)
    division_column: (str) column name 
    country: (str) Default: 'philippines'
    
    ======================================================================================
    
    Example usage:
    >>> from SarsCov2Variants import aggregate_divisions
    >>> new_df, divisions = aggregate_divisions(my_df, division_column='division')
    >>> new_df['agg_division'].value_counts()
    
    ======================================================================================
    '''
    
      
    # Enumerate all possible names of the divisions
    # You may check this with the command below:
    # >>> df[division_column].value_counts()
    philippines_dict = { 
                            'NCR'         : ['ncr', 'manila', 'national capital region','manila city','navotas city','pasay city','quezon city','taguig city','caloocan city'], 
                            'CAR'         : ['car', 'cordillera administrative region'],
                            'Region I'    : ['region i', 'ilocos region','ilocos'],
                            'Region II'   : ['region ii', 'cagayan valley'],
                            'Region III'  : ['region iii', 'central luzon'],
                            'Region IV-A' : ['region iv-a', 'calabarzon'],
                            'Region IV-B' : ['region iv-b', 'mimaropa'],
                            'Region V'    : ['region v', 'bicol region','bicol'],
                            'Region VI'   : ['region vi', 'western visayas'],
                            'Region VII'  : ['region vii', 'central visayas'],
                            'Region VIII' : ['region viii', 'eastern visayas'],
                            'Region IX'   : ['region ix', 'zamboanga peninsula','zamboanga'],
                            'Region X'    : ['region x', 'northern mindanao'],
                            'Region XI'   : ['region xi', 'davao region','davao'],
                            'Region XII'  : ['region xii', 'soccksargen'],
                            'Caraga'      : ['caraga', 'davao oriental'],
                            'BARMM'       : ['barmm', 'armm','bangsamoro autonomous region in muslim mindanao'],
                        }
                        
                        
    indonesia_dict = {
                        'Java'	                : ['banten', 'scr of jakarta', 'west java', 'central java', 'sr of yogyakarta', 'east java'],
                        'Kalimantan'	        : ['west kalimantan', 'central kalimantan', 'north kalimantan', 'east kalimantan', 'south kalimantan'],
                        'Maluku Islands'	    : ['north maluku', 'maluku'],
                        'Lesser Sunda Islands'	: ['bali', 'west nusa tenggara', 'east nusa tenggara'],
                        'Western New Guinea'	: ['west papua', 'papua'],
                        'Sulawesi'	            : ['north sulawesi', 'gorontalo', 'central sulawesi', 'west sulawesi', 'south sulawesi', 'southeast sulawesi'],
                        'Sumatra'	            : ['aceh', 'north sumatra', 'west sumatra, riau', 'riau islands', 'jambi', 'bengkulu', 'south sumatra', 'bangka belitung islands', 'lampung'],
                    }
                        
                        
    if country.lower().lstrip().rstrip() == 'philippines':
        dictionary = philippines_dict
    elif country.lower().lstrip().rstrip() == 'indonesia':
        dictionary = indonesia_dict
    else:
        print('\nError: Unknown country.\n')
        exit()
    
    df = metadata_df.copy()
    df['agg_division'] = 'Unknown'
    
    for division, alt_names in dictionary.items():
        for alt_name in alt_names:
            indices = df[df[division_column].astype(str).str.lower().str.lstrip().str.rstrip() == alt_name].index
            df.loc[indices, 'agg_division'] = division
            
    return df


    























































def area_charts_of_divisions(metadata_df, normalized=False, date_column='date',
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
    normalized: (bool) If True, the area charts are normalized per aggregation mode (weekly),
                            Default: False
    date_column: (str) column name, Default: 'date'
    division_column: (str) column name, Default: 'agg_division'
    id_column: (str) column name, Default: 'gisaid_epi_isl'    
    variant_column: (str) column name, Default:'variant'
    output_directory: (str) destination folder of the output png files.
                            
    ======================================================================================                                                                                               
    '''
    
    groupby = metadata_df[[ variant_column,
                            date_column,
                            division_column,
                            id_column]].groupby(by=[date_column,
                                                    variant_column,
                                                    division_column], as_index=False).count()  
    plt.style.use('seaborn')  
    countmax = groupby.groupby(by=[date_column,division_column])[id_column].sum().max()                                           
    
    divisions = metadata_df[division_column].unique().tolist()
    divisions = [col for col in divisions if col not in ('Unknown')]                    
    random_number = np.random.randint(9999)                                         
    for i in range(len(divisions)):
        region = divisions[i]
        reg = groupby[groupby[division_column]==region] 
    
        # If there's no data  
        if len(reg) == 0:  
            reg_pivot = pd.DataFrame({date_column:['2020-10','2021-01','2021-03'],'No data':[0,0,0]})
            reg_pivot[date_column] = reg_pivot[date_column].astype('datetime64')
            reg_pivot[date_column] = reg_pivot[date_column].dt.date
            COLORS = 'black'
   
    
        # If there's data
        else:
            reg_pivot = reg.pivot(index=reg[date_column], columns=variant_column)[id_column]
            reg_pivot = reg_pivot.replace(np.nan, 0)
            total_df = reg_pivot.sum(axis=1)
      
            # Normalize the occurences
            if normalized:
                countmax = 1
                for variant in reg_pivot.columns:
                    reg_pivot[variant] = reg_pivot[variant]/total_df
               
            # Put the 'Others' column at the end
            cols = reg_pivot.columns.tolist()
            reg_pivot = reg_pivot[cols]
      
            # Smoothen curves if there are more than 3 data points.
            if len(reg) > 3:
                df = pd.DataFrame()
                new_index = np.arange((reg_pivot.index[-1] - reg_pivot.index[0]).days + 1)
                date_index = pd.date_range( start = reg_pivot.index[0], 
                                            end = reg_pivot.index[-1],
                                            freq='D').date

                for variant in reg_pivot.columns:
                    func = interp1d((reg_pivot.index - reg_pivot.index[0]).days, reg_pivot[variant], kind='cubic')
                    df[variant] = func(new_index)
    
                df[df < 0] = 0     
                df.index = date_index
                reg_pivot = df
                
                variant_color_df = pd.DataFrame({variant_column:reg_pivot.columns.tolist()})
                variant_color_df = variant_color(variant_color_df)
                
            COLORS = variant_color_df.color.tolist()

        reg_pivot.plot( kind = 'area', 
                        rot = 0,
                        color = COLORS)

                       
        ax = plt.gca()
        handles, labels = ax.get_legend_handles_labels()
        labels, handles = zip(*sorted(zip(labels, handles), 
                                      key=lambda t: t[0], 
                                      reverse=True))
        ax.legend(handles, labels)
        plt.ylim(0, countmax)
        plt.xlim(metadata_df[date_column].min(), metadata_df[date_column].max())
        plt.xlabel(date_column)
        update = f'{random_number}' 
        if normalized:
            plt.ylabel('Count')
            update = f'N{random_number}'   
        plt.legend(loc = 'upper left') 
        plt.title('{}'.format(region), 
                        fontsize = 17, 
                        fontstyle = 'oblique') 
        plt.savefig(f'{output_directory}/{i}_{region}_{update}.png')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        

def bubble_map_of_country(metadata_df, date_column='date', division_column='agg_division',
                        variant_column='variant', id_column='gisaid_epi_isl', color_column='color',
                        latitude_column='div_latitude', longitude_column='div_longitude',
                        country='philippines',
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
    >>> from SarsCov2Variants import bubble_map_of_country
    >>> bubble_map_of_country(my_df)      
    
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
                            
    fig['layout']['yaxis']['title']='Count'
    fig.write_html(f'{output_filename}')

    

















































def combine_metadata(directory, meta_cols=['strain','gisaid_epi_isl','date','region','country',
                    'division','age','sex','pangolin_lineage','GISAID_clade','originating_lab',
                    'submitting_lab','date_submitted'], date_column='date', agg_per_week=True):
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
    agg_per_week: (bool) if True, aggregates the dates per week, else the dates remain
    
    ====================================================================================== 
       
    Example usage:
    >>> from SarsCov2Variants import combine_metadata                
    >>> combine_metadata('/User/Downloads', meta_cols=['virus','division'], date_column='date')   
    
    ====================================================================================== 
    '''
    
    metadata_df = pd.DataFrame()
    
    count = 0
    for dirname, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename == 'metadata.tsv':
                temp_df = pd.read_csv(f'{dirname}/{filename}', 
                                        sep='\t', 
                                        usecols=meta_cols)
                                                                        
                metadata_df = pd.concat([metadata_df, temp_df], ignore_index=True)
                count = count + 1

    if len(metadata_df) == 0:
        print('\nError: No data was found. Check the spelling of filenames and check if the files are not empty.\n')
        exit()
    else:
        print(f'Found {count} "metadata.tsv" files\n')
        print(f'Total GISAID sequences: {len(metadata_df)}')
        
    metadata_df.sort_values(by=date_column, ascending=False, inplace=True)
    metadata_df = metadata_df.reset_index(drop=True)
    metadata_df[date_column] = metadata_df[date_column].astype('datetime64')
    if agg_per_week:
        metadata_df[date_column] = metadata_df[date_column].dt.to_period('W').dt.start_time

    
    return metadata_df








































def geocode_divisions(metadata_df, division_column='agg_division', country='philippines'):
    '''
    ======================================================================================
    
    Description: Maps the divisions to the dictionary of latitude and longitude (WSG84) coordinates. 
    
    Returns: DataFrame with two new columns 'div_latitude' and 'div_longitude'
    
    ======================================================================================
    
    <INPUTS>
    metadata_df: (dataframe)
    division_column: (str) column name Default: 'agg_division'
    country: (str) Default: 'philippines'
    
    ======================================================================================
    
    Example usage:
    >>> from SarsCov2Variants import geocode_divisions
    >>> new_df = geocode_divisions(my_df, division_column='agg_division')
    
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
                        'Java'	                : [-6.430323, 107.500769],
                        'Kalimantan'	        : [0.529112, 114.198348],
                        'Maluku Islands'	    : [-3.180422, 129.115683],
                        'Lesser Sunda Islands'	: [-9.788183, 120.026111],
                        'Western New Guinea'	: [-3.657968, 137.944309],
                        'Sulawesi'	            : [-1.833530, 120.287056],
                        'Sumatra'	            : [-1.094433, 101.977441],
                    }
    
    
    df = metadata_df.copy()
    df['div_latitude'] = np.nan
    df['div_longitude'] = np.nan
    
    if country.lower().lstrip().rstrip() == 'philippines':
        dictionary = philippines_dict
    elif country.lower().lstrip().rstrip() == 'indonesia':
        dictionary = indonesia_dict
    else:
        print('\nError: Unknown country.\n')
        exit()
    
    for division, coordinates in dictionary.items():
        latitude = coordinates[0]
        longitude = coordinates[1]
        
        indices = df[df[division_column]==division].index
        df.loc[indices, 'div_latitude'] = latitude
        df.loc[indices, 'div_longitude'] = longitude
    
    return df





































def init_functions(directory, country='philippines', lineage_column='pangolin_lineage', 
            division_column='division', date_column='date', variant_column='variant', 
            id_column='gisaid_epi_isl', agg_per_week=True,
            meta_cols=['strain', 'gisaid_epi_isl', 'date', 'region', 
            'country', 'division', 'age', 'sex', 'pangolin_lineage', 'GISAID_clade', 
            'originating_lab', 'submitting_lab', 'date_submitted'],):
            
    '''
    ======================================================================================
    
    Description: Runs all the initial (preprocessing) functions in SarsCov2Variants.
            
    
    Returns a dataframe and saves it as a csv file.     
            
    ======================================================================================   
    '''
    
    output_filename = f'output/12_variants_in_{country}'
    
    meta_df = combine_metadata(directory, meta_cols, date_column, agg_per_week)
    meta_df = lineage_to_variant(meta_df, lineage_column)
    meta_df = aggregate_divisions(meta_df, division_column, country)
    division_column='agg_division'
    meta_df = geocode_divisions(meta_df, division_column, country)
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
                                                   
    casecount_df.columns = ['Variant', 'Location', 'Latitude', 'Longitude', 'Case Count']
    
    meta_df.to_csv(f'{output_filename}.csv', index=False)
    casecount_df.to_csv(f'{output_filename}_aggregated.csv', index=False)
    
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
    >>> from SarsCov2Variants import lineage_to_variant
    >>> lineage_to_variant(my_df, lineage_column='pangolin_lineage')
    
    ======================================================================================
    '''
    
    # ALWAYS CHECK FOR UPDATES IN THE LINK ABOVE!
    variants_dict = {
                'Alpha' : ['B.1.1.7', 'Q'], 
                'Beta'  : ['B.1.351'], 
                'Gamma' : ['P.1'], 
                'Delta' : ['B.1.617.2', 'AY'], 
                'Eta'   : ['B.1.525'], 
                'Iota'  : ['B.1.526'], 
                'Kappa' : ['B.1.617.1'], 
                'Lambda': ['C.37'],
                'Mu'    : ['B.1.621'],
                }
    
    df = metadata_df.copy()
    df['variant'] = 'Others'
    
    for variant, lineages in variants_dict.items():
        for lineage in lineages:
            indices = df[df[lineage_column].astype(str).str.startswith(lineage)].index
            df.loc[indices, 'variant'] = variant
            
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
    >>> from SarsCov2Variants import variant_color
    >>> variant_color(my_df, variant_column='variant')
    
    ====================================================================================== 
    '''
    
    colors_dict = {
                    'Alpha' : 'maroon',
                    'Beta'  : 'royalblue',
                    'Gamma' : 'lightgreen',
                    'Delta' : 'gold',
                    'Eta'   : 'orange',
                    'Iota'  : 'brown',
                    'Kappa' : 'green',
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



























