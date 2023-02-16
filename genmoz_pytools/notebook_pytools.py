import pandas as pd
import numpy as np
import geopandas

def import_genmoz_retrospective_data(data_path = "/home/isglobal.lan/apujol/isglobal/projects/genmoz/data/retrospective/", \
                                    ucsf_data_path = "/home/isglobal.lan/apujol/isglobal/projects/genmoz/data/retrospective/ucsf/", \
                                    genmoz_data_path = "/home/isglobal.lan/apujol/isglobal/projects/genmoz/data/", \
                                    pregmal_data_path = "/home/isglobal.lan/apujol/isglobal/projects/pregmal/data/", \
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
                                    ucsf_barcode2studyid_filename = 'msbarcodeToStudyID.tab.txt', \
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
    barcode2studyid = pd.read_csv(ucsf_data_path + ucsf_barcode2studyid_filename, delimiter = '\t')

    #Processing data
    #Region from Zambezia was wrongly assigned
    mask = (summary_meta['Province'] == 'Zambezia') | (summary_meta['Province2'] == 'Zambezia')
    summary_meta.loc[mask, 'Region2'] = 'Central'

    #Adding metadata
    summary_meta = pd.merge(summary_meta, all_meta[['Sample', 'WGS ENA accession', 'study', 'Date']].drop_duplicates(), on = 'Sample', how = 'left')
    summary_meta['sample_name_calls'] = pd.Series(summary_meta['Sample'])
    mask = summary_meta['WGS ENA accession'].notnull()
    summary_meta['sample_name_calls'][mask] = mask = summary_meta['WGS ENA accession']

    summary_meta['source'] = summary_meta['study']
    summary_meta['source'][summary_meta['study'] == 'XMAG15'] = 'X'
    summary_meta['source'][summary_meta['study'] == 'XMAL15'] = 'X'
    summary_meta['source'][summary_meta['study'] == 'XMAG18'] = 'X'
    summary_meta['source'][summary_meta['study'].isnull()] = 'UCSF'
    return summary_meta, all_meta, divmetrics, v4_amplicon, ucsf_data, ucsf_metadata, barcode2studyid
