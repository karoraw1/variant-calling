anchor_cols = ['primer_type', 'mutation_rows', 'support_rows', 'variant_position', 
               'support_positions', 'target_start', 'target_end', 'nucleotides', 'amino_acids']

var_df['midint'] = var_df['midpoint'].apply(int)

type_1_anchor = []
for midx in var_df['midint'].unique():
    varrows = list(var_df[var_df['midint'] == midx].index)
    varpos = [var_df.at[i, 'var_locii'] for i in varrows]
    nt_seq = full_genome[int(midx)-1-75:int(midx)-1+75]
    type_1_anchor.append(('exact', varrows, np.nan, varpos, np.nan,
                          int(midx)-75, int(midx)+75, nt_seq, np.nan))

print("{} nucleotide primers for exact mutation locations were created".format(len(type_1_anchor)))

target_locii = [j for i in list(var_df['affected_3'].dropna()) for j in i ]
support_locii = [j for i in list(var_df['support_features'].dropna()) for j in i ]
all_possible = set(target_locii).union(set(support_locii))

print("{} primers for full genes will be created".format(len(all_possible)))

type_2_anchor = []
for tf in all_possible:
    mut_rows = [i for i in var_df.index if tf in var_df.loc[i, 'affected_3']]
    sup_rows = [i for i in var_df.index if tf in var_df.loc[i, 'support_features']]
    varpos = [var_df.at[i, 'var_locii'] for i in mut_rows]
    suppos = [var_df.at[i, 'var_locii'] for i in sup_rows]
    feat_idx = list(ngb_df[ngb_df['locus_tag'] == tf].index)
    assert len(feat_idx) == 1
    feat_deets = list(ngb_df.loc[feat_idx[0], ['start', 'end', 'sequence', 'translation']])
    type_2_anchor.append(tuple([tf, mut_rows, sup_rows, varpos, suppos] + feat_deets))
    
anchors = type_1_anchor + type_2_anchor
print("A total of {} primers were created".format(len(anchors)))

anchor_df = pd.DataFrame(index=range(len(anchors)), columns=anchor_cols, data=anchors)
anchor_df['primer_idx'] = pd.Series(index=anchor_df.index, data=range(len(anchors)))
anchor_df

#############################################################################################

data_dir = "/Volumes/KeithSSD/variant-calling/variant_homologs"
if not os.path.exists(data_dir):
    os.mkdir(data_dir)
    
prot_file = os.path.join(data_dir, 'prot_primers.fa')
nt_file = os.path.join(data_dir, 'nt_primers.fa')

header_fxn = lambda x: ">"+"_".join([str(f) for f in x])
anchor_df['headers'] = anchor_df[['primer_type', 'primer_idx']].apply(header_fxn, axis=1)

make_records = lambda x: "\n".join(x)
for rec_type, file_handle in zip(['amino_acids', 'nucleotides'], [prot_file, nt_file]):
    rec_list = anchor_df[['headers', rec_type]].dropna().apply(make_records, axis=1)
    print("{} {} records written to {}".format(len(rec_list), rec_type, file_handle))
    full_str = "\n".join(list(rec_list))
    with open(file_handle, 'w') as fh:
        fh.write(full_str)

#############################################################################################

