import pandas as pd
import numpy as np
import geopandas
import matplotlib.pyplot as plt
from matplotlib import cm
from genomic_tools.microhaplotypes import He_from_samples, get_allele_frequencies

def import_genmoz_retrospective_data(data_path = "/home/apujol/isglobal/projects/genmoz/data/retrospective/", \
                                    ucsf_data_path = "/home/apujol/isglobal/projects/genmoz/data/retrospective/ucsf/", \
                                    genmoz_data_path = "/home/apujol/isglobal/projects/genmoz/data/", \
                                    pregmal_data_path = "/home/apujol/isglobal/projects/pregmal/data/", \
                                    summary_meta_filename = "ALL_2015_2018.lowsamples_50percent.metadata.txt", \
                                    all_meta_filename = "all_wt_grc_react_merged_clean.csv", \
                                    all_meta_txt_filename = "mod_all_wt_grc_merged_clean.txt", \
                                    divmetrics_filename = "diversityMeasuresPerTarget.tab.txt", \
                                    v4_amplicon_filename = "v4_calls_WGS.tab.txt", \
                                    codes_WT_ISG2_filename = 'SANGER_shipment_PAU_ALL_june2021_2changes.dta', \
                                    grc1_filename = 'Genotyping_results_March2020_1234-PF-MZ-MAYOR.xlsx', \
                                    grc2_filename = 'Genotyping_results_May_2021_1234-PF-MZ-MAYOR.xlsx', \
                                    ucsf_data_filename = 'EpidemiologiaGenetic_DATA_2018-07-22_0025_V1.csv', \
                                    ucsf_metadata_filename = 'MozambiqueMicrosatellite_Metadata.csv', \
                                    ucsf_barcode2studyid_filename = 'metadata_Arnau.csv', \
                                    index_filename = "index_nida.csv", contacts_filename = "contacts_nida.csv"):
    #It conains metadata of index and contact cases
    index = pd.read_csv(pregmal_data_path + index_filename)
    contacts = pd.read_csv(pregmal_data_path + contacts_filename)
    ##It contains a table connecting nida with external id
    codes_WT_ISG = pd.read_stata(pregmal_data_path + codes_WT_ISG2_filename)
    #GRC data
    grc1 = pd.read_excel(pregmal_data_path + grc1_filename)
    grc2 = pd.read_excel(pregmal_data_path + grc2_filename)

    #Assigning sample ID to contact data
    #Getting table linking nida with sample ID
    codes_WT_ISG1 = pd.merge(codes_WT_ISG[['external_id', 'nida']][codes_WT_ISG['nida'].notnull()], \
         grc1[['Sample Internal ID', 'Sample External ID']], left_on = 'external_id', \
        right_on = 'Sample External ID', how = 'left')
    codes_WT_ISG2 = pd.merge(codes_WT_ISG1, \
         grc2[['Sample Internal ID', 'Sample External ID']], left_on = 'external_id', \
        right_on = 'Sample External ID', how = 'left')
    codes_WT_ISG2['Sample Internal ID_x'][codes_WT_ISG2['Sample Internal ID_x'].isnull()] = codes_WT_ISG2['Sample Internal ID_y'][codes_WT_ISG2['Sample Internal ID_x'].isnull()]
    codes_WT_ISG2['Sample External ID_x'][codes_WT_ISG2['Sample External ID_x'].isnull()] = codes_WT_ISG2['Sample External ID_y'][codes_WT_ISG2['Sample External ID_x'].isnull()]
    codes_WT_ISG2['Sample Internal ID'] = codes_WT_ISG2['Sample Internal ID_x']
    codes_WT_ISG2['Sample External ID'] = codes_WT_ISG2['Sample External ID_x']

    contacts_sampleid = pd.merge(contacts, codes_WT_ISG2[['Sample Internal ID', 'external_id', 'nida']][codes_WT_ISG['nida'].notnull()], \
                             on = 'nida', how = 'left')

    #Loading data
    summary_meta = pd.read_csv(data_path + summary_meta_filename, sep = '\t')
    all_meta = pd.read_csv(data_path + all_meta_filename)
    all_meta = all_meta.rename(columns = {'Sample Internal ID' : 'Sample'})
    divmetrics = pd.read_csv(data_path + divmetrics_filename, sep = '\t')
    v4_amplicon = pd.read_csv(data_path + v4_amplicon_filename, sep = ' ')
    v4_amplicon['s_collection_date'] = pd.to_datetime(v4_amplicon['s_collection_date'])

    #UCSF metadata
    ucsf_data = pd.read_csv(ucsf_data_path + ucsf_data_filename)
    ucsf_data['date'] = pd.to_datetime(ucsf_data['date'])
    ucsf_data['date5'] = pd.to_datetime(ucsf_data['date5'])
    ucsf_data['age'][ucsf_data['age'] == '7Chichango'] = '7'
    ucsf_data['age'] = np.array(ucsf_data['age'], dtype = float)
    ucsf_metadata = pd.read_csv(ucsf_data_path + ucsf_metadata_filename)
    ucsf_metadata = geopandas.GeoDataFrame(ucsf_metadata, geometry=geopandas.points_from_xy(ucsf_metadata['long'], ucsf_metadata['lat']))
    ucsf_metadata['date'] = pd.to_datetime(ucsf_metadata['date'])
    ucsf_metadata['date5'] = pd.to_datetime(ucsf_metadata['date5'])
    ucsf_metadata['age'][ucsf_metadata['age'] == '7Chichango'] = '7'
    ucsf_metadata['age'] = np.array(ucsf_metadata['age'], dtype = float)
    barcode2studyid = pd.read_csv(ucsf_data_path + ucsf_barcode2studyid_filename)

    #Processing data
    #Region from Zambezia was wrongly assigned
    mask = (summary_meta['Province'] == 'Zambezia') | (summary_meta['Province2'] == 'Zambezia')
    summary_meta.loc[mask, 'Region2'] = 'Central'

    #Adding metadata from UCSF samples
    ucsf_data = pd.merge(ucsf_data, barcode2studyid[barcode2studyid['id'].notnull()], \
                        left_on = 'studyid', right_on = 'id', how = 'left')
    ucsf_columns = ['studyarea', 'date', 'clinic', 'age', \
                 'genre', 'province0',  'neighborhood', \
                 'didyouspendnight', 'numbenight', 'datearrival', \
                 'destiny', 'district', 'administrativepost', \
                 'locality', 's_Sample']

    summary_meta = pd.merge(summary_meta, all_meta[['Sample', 'WGS ENA accession', 'study', 'Date']].drop_duplicates(), on = 'Sample', how = 'left')
    summary_meta['sample_name_calls'] = pd.Series(summary_meta['Sample'])
    mask = summary_meta['WGS ENA accession'].notnull()
    summary_meta['sample_name_calls'][mask] = mask = summary_meta['WGS ENA accession']

    summary_meta['source'] = summary_meta['study']
    summary_meta['source'][summary_meta['study'] == 'XMAG15'] = 'X'
    summary_meta['source'][summary_meta['study'] == 'XMAL15'] = 'X'
    summary_meta['source'][summary_meta['study'] == 'XMAG18'] = 'X'
    summary_meta['source'][summary_meta['study'].isnull()] = 'UCSF'

    ucsf_data['s_Sample'] = pd.Series(np.array(np.array(ucsf_data['s_Sample'], \
                                      dtype = int), dtype = str), dtype = object)
    #Avoiding duplicate values of ucsf_data[Sample] == 8034211758
    ucsf_data = ucsf_data.drop(177)
    summary_meta = pd.merge(summary_meta, ucsf_data[ucsf_columns], \
                            left_on = 'Sample', right_on = 's_Sample', \
                            how = 'left')
    return summary_meta, all_meta, divmetrics, v4_amplicon, ucsf_data, ucsf_metadata, barcode2studyid

def import_HFS22_data(data_path = "/home/apujol/isglobal/projects/genmoz/data/HFS/", diversity_filtered = True):
    """
    This method loads the HFS 2022 GenMoz metadata and NextSeq data. 
    
    Parameters:
    -----------
    data_path: str
        The directory where the data is stored. 
    diversity_filtered: bool
        If True, the data selected is the only restricted to diversity loci with applied covering and read filtering. 
    
    Returns:
    --------
        hfs_metadata: pd.DataFrame
            Metadata of the HFS collected data. 
        hfs_run_data: pd.DataFrame
            Data of the NextSeq results
    """
    hfs_excel_filename = data_path + 'HFS_ANO1.xlsx'
    if diversity_filtered:
        hfs_run_filename = data_path + 'HFS22_NextSeq_allele_data_filtered_div.csv'
    else: 
        hfs_run_filename = data_path + 'HFS22_NextSeq_allele_data.csv'
    
    hfs_metadata = pd.read_excel(hfs_excel_filename)
    hfs_metadata = rename_all_xy(hfs_metadata, suf_x = '.x', suf_y = '.y')
    hfs_run_data = pd.read_csv(hfs_run_filename)
    hfs_metadata['date'] = pd.to_datetime(hfs_metadata['date'])
    hfs_metadata['travel_date'] = pd.to_datetime(hfs_metadata['travel_date'])
    hfs_metadata['interview_date'] = pd.to_datetime(hfs_metadata['interview_date'])
    return hfs_metadata, hfs_run_data

def prepare_data4dcifer(v4_amplicon, summary_meta, save = False, \
                        output_path = "/home/isglobal.lan/apujol/isglobal/projects/genmoz/notebooks/testing/data/", \
                       diversity_name = "v4_call_diversity_retrospective_data.csv", \
                       all_amplicons_name = "v4_call_retrospective_data.csv"):
    """
    This method generates data and files for the computation of IBD with Dcifer.

    Parameters:
    -----------
    v4_amplicon: pd.DataFrame
        Amplicon data of samples
    summary_meta: pd.DataFrame
        Metadata of samples
    save: bool
        It specifies if the output data is saved
    output_path: str
        Directory where output data is saved
    diversity_name: str
        Name of file with data with only microhaplotypes of diversity
    all_amplicons_name: str
        Name of file with data containing all the microhaplotypes

    Return:
    -------
    data4dcifer: pd.DataFrame
        Dataframe containing the data with all the microhaplotypes
    data4dcifer_diversity: pd.DataFrame
        Dataframe containing the data the microhaplotypes of diversity
    """
    #Preparing data for Dcifer
    samples4dcifer = summary_meta[['Sample']].drop_duplicates()
    data4dcifer = pd.merge(v4_amplicon, samples4dcifer, left_on = 's_Sample', right_on = 'Sample', how = 'inner')
    diversity_amplicons = pd.Series([i for i in v4_amplicon['amplicon'].unique() if i[-2:] == '1A'])
    diversity_amplicons = pd.DataFrame(diversity_amplicons, columns = ['amplicon'])
    data4dcifer_diversity = pd.merge(data4dcifer, diversity_amplicons, on = 'amplicon', how = 'inner')
    if save:
        data4dcifer_diversity.to_csv(output_path + diversity_name)
        data4dcifer.to_csv(output_path + all_amplicons_name)
    return data4dcifer, data4dcifer_diversity

def get_data_for_glm(ibd_values, dist_values, p_values, min_IBD, max_p, min_dist, max_dist):
    """
    This method filters the data of IBD, distances and p-values for the linear regression
    analysis.

    Parameters:
    -----------
    bd_values: np.array or pd.Series
        Array of IBD values of the pairs.
    dist_values: np.array or pd.Series
        Array of the distance between the pairs.
    p_values: np.array or pd.Series
        Array of p-values of the IBD of the pairs.
    min_IBD: float
        Minimum IBD from which the fraction was calculated.
    max_p: float
        Maximum p-value from which the fraction was calculated.
    min_dist: float
        Minimum distance to include in the bins.
    max_dist: float
        Maximum distance to include in the bins.

    Returns:
    --------
    dist_filt: np.array
        Array of distances included.
    ibd_high: np.array
        Binary array defining the cases that passed the IBD and p-value thresholds.
    """
    if min_dist is None:
        min_dist = np.min(dist_values)
    if max_dist is None:
        max_dist = np.max(dist_values)
    ibd_high = (ibd_values >= min_IBD)&(p_values <= max_p)
    ibd_high = np.array(ibd_high, dtype = int)
    dist_mask = (dist_values >= min_dist)&(dist_values <= max_dist)
    ibd_high = ibd_high[dist_mask]
    dist_filt = dist_values[dist_mask]
    return dist_filt, ibd_high

def rename_all_xy(df, suf_x = '_x', suf_y = '_y'):
    """
    This method identifies the duplicated variables resulting from a merge and removes
    the duplicates. 

    Parameters:
    -----------
    df: pd.DataFrame
        Dataframe to use.
    suf_x: str
        Suffix used to identify one duplication of the variables.
    suf_y: str
        Suffix used to identify the other duplication of the variables.
    
    Returns:
    --------
    df: pd.DataFrame
        Dataframe with the dropped duplicates.
    """
    list_columns = find_all_xy(df, suf_x = suf_x, suf_y = suf_y)
    for name in list_columns:
        df[name+suf_x][df[name+suf_x].isnull()] = df[name+suf_y][df[name+suf_x].isnull()]
        df = df.rename(columns = {name+suf_x:name})
        del df[name+suf_y]
    return df

def find_all_xy(df, suf_x = '_x', suf_y = '_y'):
    """
    This method identifies the columns from a dataframe containing 
    the suffixes suf_x and suf_y. 

    Parameters:
    -----------
    df: pd.DataFrame
        Dataframe to use.
    suf_x: str
        Suffix used to identify one duplication of the variables.
    suf_y: str
        Suffix used to identify the other duplication of the variables.
    
    Returns:
    --------
    list_columns: list
        List of the identified columns.
    """
    list_columns = []
    for col in df:
        if col[-2:] == suf_x:
            if col[:-2] + suf_y not in df.columns:
                print("Warning: " + col[:-2] + suf_y + ' not found')
            list_columns.append(col[:-2])
    return list_columns

list_locs = {
    'Manhica' : [32.80722050, -25.40221980], #From Pau Cisteró excel
    'Maputo' : [32.576388888888889, -25.915277666666],
    'Montepuez' : [38.99972150, -13.12555980], #From Pau Cisteró excel
    'Chokwe' : [33.005166666666666667, -24.5252777777777777],
    'Moatize' : [33.73333040, -16.11666620], #From Pau Cisteró excel
    'Dondo' : [34.75, -19.6166666666666667],
    'Magude' : [32.64216410, -25.02049992], #From Pau Cisteró excel
    'Ilha Josina' : [32.92210000, -25.09330000], #From Pau Cisteró excel
    'Xinavane' : [32.791885, -25.048534],
    'Panjane' : [32.352430, -24.899469],
    'Motaze' : [32.860569, -24.810357],
    'Mapulanguene' : [32.081602, -24.491015],
    'Taninga' : [32.825796, -25.182094],
    'Palmeira' : [32.869766, -25.261457],
    'Massinga' : [35.37405260,-23.32666250], #From Pau Cisteró excel
    'Mopeia' : [35.71338490, -17.97391000], #From Pau Cisteró excel
    'Gaza' : [34.19153, -24.9206],#Chidenguele HF
    'Inhambane' : [35.38, -23.33456],#Massinga
    'Zambezia' : [35.71279, -17.97899],#Mopeia
    'C.Delgado' : [38.99972150, -13.12555980],#Montepuez
    'Cabo-Delgado' : [38.99972150, -13.12555980],#Montepuez
    'Tete' : [33.618156, -16.138187],#Tete
    'Sofala' : [34.846280, -19.833158],#Beira
    }
locations = pd.DataFrame({'location' : [i for i in list_locs], 'longitude': [list_locs[i][0] for i in list_locs], 'latitude': [list_locs[i][1] for i in list_locs]})
locations = geopandas.GeoDataFrame(locations, geometry = geopandas.points_from_xy(locations['longitude'], locations['latitude']))
locations = locations.set_crs(epsg=4326)

ucsf_destiny2province = {
                        1.0 : 'Maputo', \
                        2.0 : 'Gaza', \
                        3.0 : 'Inhambane', \
                        4.0 : 'Sofala', \
                        5.0 : 'Manica', \
                        6.0 : 'Nampula', \
                        7.0 : 'Cabo Delgado', \
                        8.0 : 'Niassa', \
                        9.0 : 'Another Country', \
}

react_destiny2location = {1.0 : 'Magude', \
                        2.0 : 'Xinavane', \
                        3.0 : 'Manhiça', \
                        4.0 : 'Maputo City', \
                        5.0 : 'Maputo Province', \
                        6.0 : 'Gaza', \
                        7.0 : 'Inhambane', \
                        8.0 : 'Sofala', \
                        9.0 : 'Manica', \
                        10.0 : 'Zambezia', \
                        11.0 : 'Tete', \
                        12.0 : 'Nampula', \
                        13.0 : 'Cabo Delgado', \
                        14.0 : 'Niassa', \
                        15.0 : 'South Africa', \
                        16.0 : 'Eswatini', \
                        17.0 : 'Another Country', \
}

def get_overall_exp_He_vs_size(amplicon_data, labels, xlim = None, \
                              ylim = None, locus_name = 'p_name', \
                                allele_name = 'h_popUID', sample_name = 's_Sample', \
                                freq_name = 'c_AveragedFrac'):
    """
    This method plots the relationship between overall
    expected He and sample size from different sub-selections
    of the samples.

    Parameters:
    -----------
    amplicon_data: pd.DataFrame
        Data frame containing the allele frequencies of samples
        and the category label variable.
    lables: list
        List of labels to use to apply different selections.
    xlim: list [float, float]
        Limits of x-axis.
    ylim: list [float, float]
        Limits of y-axis.
    locus_name: str
        Name of the column referring to the locus name. 
    allele_name: str
        Name of the column referring to the allele name. 
    sample_name: str
        Name of the column referring to the sample ID name. 
    freq_name: str
        Name of the column referring to the allele frequency name. 

    Returns:
    --------
    Error bar plot with overall expected He versus sample size.
    """
    #Expected He vs sample size
    for label in labels:
        categories = amplicon_data[label][amplicon_data[label].notnull()].unique()
        He_per_cat = {}
        He_per_cat_err = {}
        for cat in categories:
            mask = amplicon_data[label] == cat
            size = len(amplicon_data[mask][sample_name].unique())
            loci_He, overall_He = He_from_samples(amplicon_data[mask], locus_name = locus_name, \
                                                  allele_name = allele_name, \
                                                  freq_name = freq_name)
            He_per_cat[cat] = overall_He
            He_per_cat_err[cat] = np.std(loci_He['He'])/np.sqrt(len(loci_He['He']))
            plt.errorbar(size, He_per_cat[cat], He_per_cat_err[cat], c = 'tab:blue', marker = 'o')
    plt.grid()
    plt.xlabel('Sample size')
    plt.ylabel('Mean expected He')
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.show()

def plot_He_per_cat(He_per_cat, He_per_cat_err, colours = None):
    """
    This method plots the overall expected He per
    category.

    Parameters:
    -----------
    He_per_cat: pd.DataFrame
        Data frame with the values of overall expected
        He for each category.
    He_per_cat_err: pd.DataFrame
        Data frame with the values of overall expected
        He error for each category.
    colours: list
        List of colours for plotting each category.

    Returns:
    --------
    Error bar plot showing the statistics.
    """
    categories = list(He_per_cat.keys())
    if colours is None:
        colours = [cm.turbo(i/len(categories)) for i in range(len(categories))]
    for i, cat in enumerate(He_per_cat):
        plt.errorbar(i, He_per_cat[cat], He_per_cat_err[cat], marker = 'o', color = colours[i])
    plt.xticks(range(len(He_per_cat)), labels = categories)
    plt.ylabel('Mean expected He')
    plt.show()

def get_He_vs_cat(amplicon_data, label, categories = None, \
                  verbose = True, show = True, ymax = 4, \
                    show_errorbar = True, show_allele_freq = True, \
                    locus_name = 'p_name', allele_name = 'h_popUID', \
                    sample_name = 's_Sample', freq_name = 'c_AveragedFrac'):
    """
    This method calculates the expected He per locus and
    overall for different categories defined by a label.

    Parameters:
    -----------
    amplicon_data: pd.DataFrame
        Data frame containing the allele frequencies of samples
        and the category label variable.
    label: str
        Category name to which amplicon_data is selected.
    categories: list
        List of names from which the label is defined.
    verbose: bool
        If True, some informative text is printed.
    show: bool
        If True, the exp He distributions are plotted.
    ymax: float
        Y axis upper limit to show in the plot.
    show_errorbar: bool
        If True, an errorbar of the overall exp He per
        category is shown.
    show_allele_freq: bool
        If True, an the spectrum of allele frequencies per
        category is shown.
    locus_name: str
        Name of the column referring to the locus name. 
    allele_name: str
        Name of the column referring to the allele name. 
    sample_name: str
        Name of the column referring to the sample ID name. 
    freq_name: str
        Name of the column referring to the allele frequency name. 

    Returns:
    --------
    He_per_cat: pd.DataFrame
        Data frame with the values of overall expected
        He for each category.
    He_per_cat_err: pd.DataFrame
        Data frame with the values of overall expected
        He error for each category.
    """
    if categories is None:
        categories = amplicon_data[label][amplicon_data[label].notnull()].unique()
    colours = [cm.turbo(i/len(categories)) for i in range(len(categories))]
    
    #Calculating exp He per category
    He_per_cat = {}
    He_per_cat_err = {}
    for i, cat in enumerate(categories):
        mask = amplicon_data[label] == cat
        if verbose:
            print("Sample size " + str(cat) + ":" + \
                  str(len(amplicon_data[mask][sample_name].unique())))
        loci_He, overall_He = He_from_samples(amplicon_data[mask], locus_name = locus_name, \
                                              allele_name = allele_name, \
                                              freq_name = freq_name)
        He_per_cat[cat] = overall_He
        He_per_cat_err[cat] = np.std(loci_He['He'])/np.sqrt(len(loci_He['He']))
        if show:
            plt.figure(0)
            #Plotting He
            loci_He['He'].plot.density(bw_method = .1, color = colours[i], alpha = 1)
            plt.vlines(overall_He, 0, ymax, color = colours[i], lw = 2, linestyle = '--')
            plt.annotate(r"Overall H$_e$ in " + cat + " = " + str(round(overall_He,3)), \
                         xy = [overall_He +.01, .95*ymax - i/ymax], color = colours[i])
        if show_allele_freq:
            allele_freq = get_allele_frequencies(amplicon_data[mask], locus_name = locus_name, \
                                                 allele_name = allele_name, \
                                                 freq_name = freq_name)
            plt.figure(1)
            allele_freq['allele_freq'].plot.density(bw_method = .005, color = colours[i], \
                                                    alpha = .5, label = cat)
    if show:
        plt.figure(0)
        plt.ylabel('Probability density')
        plt.xlabel('Expected He')
        plt.xlim(0,1)
        plt.ylim(0, ymax)
    if show_allele_freq:
        plt.figure(1)
        plt.ylabel('Number of alleles')
        plt.xlabel('Allele frequencies')
        plt.yscale('log')
        plt.xlim(0,1)
        plt.yscale('log')
        plt.ylim(.01,100)
        plt.legend()
    if show or show_allele_freq:
        plt.show()
    if show_errorbar:
        plot_He_per_cat(He_per_cat, He_per_cat_err)
    return He_per_cat, He_per_cat_err

def sampleID2nida(dataframe):
    """
    This method creates a variable with the nida values from
    the names of the sampleIDs of the data. 
    
    Parameters:
    -----------
    dataframe: pd.DataFrame
        Data frame with the sample data
    
    Returns:
    --------
    dataframe: pd.DataFrame
        Data frame with the sample data with the new nida variable
    """
    dataframe['nida'] = pd.Series([], dtype = float)
    for i in dataframe.index:
        if dataframe['sampleID'][i][1:8].isdigit():
            if dataframe['sampleID'][i][9].isdigit():
                dataframe['nida'][i] = float(dataframe['sampleID'][i][1:8] + '.' + dataframe['sampleID'][i][9])
            elif dataframe['sampleID'][i][9] == 'S': 
                dataframe['nida'][i] = float(dataframe['sampleID'][i][1:8] + '.0')
            else:
                print("nida error in index", str(i), dataframe['sampleID'][i][1:8] + '.' + dataframe['sampleID'][i][9])
    return dataframe

def import HF_data(data_path = '/home/apujol/isglobal/projects/genmoz/data/', \
                   filename = 'GEN_MOZ_Health_Facilities.xlsx'):
    """
    This method loads and renames the data of heath facilities. 

    Parameters:
    -----------
    data_path: str
        Directory where the database is located.
    filename: str
        Data file name. 
    
    Returns:
    --------
    hf_data: pd.DataFrame
        Dataframe of HF data. 
    """
    hf_filename = data_path + filename
    hf_data = pd.read_excel(hf_filename)
    hf_data = geopandas.GeoDataFrame(hf_data, geometry = geopandas.points_from_xy(hf_data['log'], hf_data['lat']))
    hf_data = hf_data.set_crs(epsg=4326)

    #Renaming HFs
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Lua Lua'] = 'Lua lua'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Posto Campo'] = 'Posto campo'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'CS\xa0Mopeia Sede'] = 'Mopeia sede'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Guro Chivuli'] = 'C.S. Chivuli'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'C.S. Mabil'] = 'C.S. Mabil'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'C.S. Chichuco'] = 'C.S. Chicuque' #right?
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Chemba Catulene'] = 'C.S. Catulene' #right?
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Chemba Mulima'] = 'C.S. Mulima' #right?
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Chemba Chiramba'] = 'C.S. Chiramba' #right?
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Chemba Cado'] = 'C.S. Cado' #right?
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Chemba Senhabuzua'] = 'C.S. Senhabuzua' #right?
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'CHEMBA Chemba Sede'] = 'C.S. Chemba Sede'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Guro CS Guro Sede'] = 'C.S. Guro' #right?
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'C.S. Magude Sede'] = 'CS Magude'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Manhica CS Malavela'] = 'CS Malavela'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Nhlamanculu Chamanculo'] = 'CS Chamanculo' #right?
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'BOANE CS Boane'] = 'C.S. Boane' #there are more Boanes
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'BOANE CS Massaca II'] = 'C.S. de Massaca' #right? Or Manaca II?
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'KaMavota Romao'] = 'CS Romão' #right?
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Nlhamankulo CS Jose Macamo'] = 'CS José Macamo' #right?
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Nhlamanculu Xipamanine'] = 'CS Xipamanine' #right?
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Gondola H.D de Gondola'] = 'H.D. Gondola'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'KaMaxakeni CS 1 de Maio'] = 'CS 1 de Maio' #right?
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Gondola CS Francismo Manhanga'] = 'F.Manhanga'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Guro CS Guro Sede'] = 'C.S. Guro sede'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'CUAMBA Adine III'] = 'C.S. Adine 3'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Cuamba CS Etatara'] = 'C.S. Cuamba (CB)' #probably not this Cuamba
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'CUAMBA CS Mitucué'] = 'C.S. Mitucue'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Cuamba TITIMANE'] = 'C.S. Titimane'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Cuamba Mepessene'] = 'C.S. Mepessene'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Cuamba Etatara'] = 'C.S. Etatara'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'C.S. Ratorga'] = 'C.S. Ratorga'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Gondola Josina Machel'] = 'C.S. Josina Machel'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'Guro Nhansana'] = 'C.S. Nhansana'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'C.S. Rio das Pedras'] = 'C.S. Rio de Pedras'
    hf_data['Health facility'].loc[hf_data['Health facility'] == 'KaMavota 1 de Junho'] = 'CS 1 de Junho'

    return hf_data