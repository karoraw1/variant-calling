import pandas as pd
import numpy as np
import os
from Bio import SeqIO

gb_file = "/Volumes/KeithSSD/variant-calling/refs/GCF_000171795.2_ASM17179v2_genomic.gbff"

gb_records = {i.name:i for i in list(SeqIO.parse(open(gb_file,"r"), "genbank"))}

gb_subdfs = []
for chrom, gb_record in gb_records.items():
    feat_idx = list(range(len(gb_record.features)))
    exp_cols = set(['molecule', 'type', 'strand', 'start', 'end', 'sequence'])
    exp_cols.update([k for i in feat_idx for k in gb_record.features[i].qualifiers.keys()])
    gb_subdf = pd.DataFrame(index=feat_idx, columns=exp_cols)
    for i in gb_subdf.index:
        feature_i = gb_record.features[i]
        gb_subdf.at[i, 'molecule'] = chrom
        gb_subdf.at[i, 'type'] = feature_i.type
        gb_subdf.at[i, 'strand'] = feature_i.location.strand
        gb_subdf.at[i, 'start'] = int(feature_i.location.start)
        gb_subdf.at[i, 'end'] = int(feature_i.location.end)
        # This already does reverse complementing of crick strand features
        gb_subdf.at[i, 'sequence'] = str(feature_i.extract(gb_record.seq))
        for qk, qv in feature_i.qualifiers.items():
            if len(qv) != 1:
                qv = ", ".join(qv)
            gb_subdf.at[i, qk] = qv[0]
    gb_subdfs.append(gb_subdf.copy())
    print("Parsed {} features from {} to a dataframe size {}".format(len(feat_idx), chrom, gb_subdf.shape))
    
gb_df = pd.concat(gb_subdfs, axis=0, ignore_index=True)
print("Full genbank dataframe sized {}".format(gb_df.shape))

var_file = '/Volumes/KeithSSD/variant-calling/misc_data/Crofts_Gene_Variants.xlsx'
var_df = pd.read_excel(var_file, sheet_name='Spec_Variants', index_col=None)
var_df['Variant Position'] = var_df['Variant Position'].astype(str).str.replace("\.\.", ',')
var_df['Variant Position'] = var_df['Variant Position'].str.replace("^", ',')

full_genome = gb_df.loc[0, 'sequence']

lociii = []
for var_i in var_df.index:
    locii = str(var_df.loc[var_i, 'Variant Position']).replace("*", "").split(",")
    locii = tuple([int(i) for i in locii])
    lociii.append(locii)
    detectable_muts = var_df.loc[var_i, 'Variant Type'] in ['Deletion', 'SNV']
    if len(locii) == 1 and detectable_muts:
        assert full_genome[locii[0]-1] == var_df.at[var_i, 'Reference Sequence']
    elif len(locii) == 2 and detectable_muts:
        assert full_genome[locii[0]-1:locii[1]] == var_df.at[var_i, 'Reference Sequence']

var_df['var_locii'] = pd.Series(index=var_df.index, data=lociii)

var_df['Variant Location'] = var_df['Variant Location'].apply(lambda x: x.strip().upper())

# fix formatting
var_df['affected'] = var_df["Affected CG8421 Gene(s)"].str.replace(", Promoter of ", " & ")
var_df['affected'] = var_df['affected'].apply(lambda x: x.split(" & "))

# manually verified = [12, 13, 14, 25, 26, 29, 30, 41, 42]
updated_tags = {"CJ8421_RS03020": "CJ8421_RS08905", 
                "CJ8421_RS03305": "CJ8421_RS08915", 
                "CJ8421_RS06480": "CJ8421_RS08860",
                'CJ8421_RS06580': 'CJ8421_RS08870',
                'CJ8421_RS06585': 'CJ8421_RS08870',
                "CJ8421_RS07090": "CJ8421_RS07085"}
def fix_tags(x):
    y = []
    for ip in x:
        i = ip.strip()
        if not i in updated_tags.keys():
            y.append(i)
        else:
            y.append(updated_tags[i])
    return list(set(y))

var_df['affected_2'] = var_df['affected'].apply(fix_tags)
changes_made = var_df[['affected', 'affected_2']].apply(lambda x: len(set(x[0]) - set(x[1])), axis=1).sum()
print("{} locus tags were either changed or removed".format(changes_made))

for var in ['start', 'end', 'strand']:
    ngb_df[var] = ngb_df[var].astype(int)

anchor_dicts = {}
anch_counter = 0 
# location window
for var_i in var_df.index:
     mut_pos = var_df.loc[var_i, 'var_locii'][0]
     mut_us = mut_pos + 150
     mut_ds = mut_pos - 150
     
     coding_bool = ngb_df['type'] == 'CDS'
     within_bool = (ngb_df['start'] <= mut_pos) & (ngb_df['end'] >= mut_pos) & (coding_bool)
     downstream = (ngb_df['start'] > mut_pos) & (ngb_df['start'] <= mut_us) & (coding_bool)
     upstream = (ngb_df['end'] < mut_pos) & (ngb_df['end'] >= mut_ds) & (coding_bool)
     sensor_bools = [(downstream) & (ngb_df['strand'] == 1),
                     (upstream) & (ngb_df['strand'] == 1), 
                     (downstream) & (ngb_df['strand'] == -1),
                     (upstream) & (ngb_df['strand'] == -1), 
                     within_bool]
     
     sensor_types = ['upstream', 'downstream', 'upstream', 'downstream', 'in_CDS']
     keys_needed = ['locus_tag', 'start', 'end', 'strand', 'product', 'sequence', 'pseudo']
     all_pinged = within_bool.sum()+downstream.sum()+upstream.sum()

     print("{}. {} possible genes affected by mutation at {}".format(var_i, all_pinged, mut_pos))
     print("\t= {}".format(var_df.loc[var_i, 'affected_2']))
     print("{}, {}".format(mut_pos, ngb_df.loc[within_bool | upstream | downstream, keys_needed].T))
     input()
     
     for sensor, senseType in zip(sensor_bools, sensor_types):
         if sensor.sum() > 0:
             sensee = ngb_df.loc[sensor, keys_needed]
             for hit_i in sensee.index:
                 pos_dict = {'variation_index':var_i, 'mut_start':mut_pos, 'var_location': senseType}
                 for kn in keys_needed:
                     pos_dict[kn] = sensee.loc[hit_i, kn]
                 if pos_dict['strand'] == 1:
                     pos_dict['upstream_seq'] = full_genome[pos_dict['start']-20-1:pos_dict['start']+180]
                     pos_dict['downstream_seq'] = full_genome[pos_dict['end']-20-1:pos_dict['end']+180]
                     pos_dict['upstream_region'] = (pos_dict['start']-20, pos_dict['start']+180)
                     pos_dict['downstream_region'] = (pos_dict['end']-20, pos_dict['end']+180)
                 else:
                     pos_dict['upstream_seq'] = full_genome[pos_dict['end']-20-1:pos_dict['end']+180]
                     pos_dict['downstream_seq'] = full_genome[pos_dict['start']-20-1:pos_dict['start']+180]
                     pos_dict['upstream_region'] = (pos_dict['end']-20, pos_dict['end']+180)
                     pos_dict['downstream_region'] = (pos_dict['start']-20, pos_dict['start']+180)
                 anchor_dicts[anch_counter] = pos_dict
                 anch_counter += 1




# pull out start, end, locus_tag, strand, and sequence of any gene hit 
# add whether it was in one of the three categories





