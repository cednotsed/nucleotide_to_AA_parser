import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO, Seq
from pathlib import Path
cwd = Path.cwd() / 'avengers2'  # Edit paths as appropriate
dnds = cwd / 'dnds'
df = pd.read_csv(dnds / 'gisaid_cov2020_sequences.QC.human.nextstrain_filter.QC.NSmask.subset.noambig.aln.BASECOUNTS.csv')
df = df.loc[:, df.any(0)]

meta = pd.read_csv(dnds / 'orf_positions.csv')

# Remove UTRs
meta = meta.loc[meta['name'] != '5UTR', :]
meta = meta.loc[meta['name'] != '3UTR', :]

# ONE BASED INDEXING
df.index = df.index + 1

# Get separate ORFs for translation
orf_dict = {}

for i in meta.index:
    orf_dict[meta.loc[i, "name"]] = df.loc[meta.loc[i, "start"]:meta.loc[i, "end"], :]

to_save = pd.DataFrame(columns=['bp', 'codon_number', 'pos_in_codon', 'snp_count', 'Ref', 'snp',
                                'ref_AA', 'variant', 'NS/S', 'ORF'])

for orf_key in orf_dict.keys():
    orf = orf_dict[orf_key]
    # Get ref genome pos
    bp = pd.Series(orf.index, name='genome_pos')
    bp.index = bp.index + 1

    # Get codon position
    import math
    codon_number = pd.Series(bp.index / 3)
    codon_number = codon_number.apply(math.ceil)

    pos_in_codon = pd.Series(bp.index % 3).replace({0: 3})

    # Get reference protein sequence
    seq = ''.join(list(orf['Ref']))
    seq = seq.upper()
    seq = Seq.Seq(seq)
    prot_seq = seq.translate()
    prot_seq = pd.Series(list(prot_seq))
    prot_seq.index = prot_seq.index + 1

    # match AA to codon number
    ref_AA = codon_number.apply(func=lambda x: prot_seq[x])  # -1 to match 0-based indexing

    to_concat = pd.DataFrame({'bp': bp.reset_index(drop=True), 'codon_number': codon_number, 'pos_in_codon': pos_in_codon, 'ref_AA': ref_AA})
    orf_df = pd.concat([to_concat, orf.reset_index(drop=True)], axis=1)

    # Get variant positions
    def is_variant(x):
        reference_base = x['Ref']

        if x[reference_base] + x['n'] != 12626:
            return 1
        else:
            return 0

    orf_df.insert(loc=0, column='is_variant', value=orf_df.apply(is_variant, axis=1))

    # Get dict of codons
    gb = orf_df.groupby('codon_number')
    codon_dict = {}

    for key in list(gb.groups.keys()):
        codon = gb.get_group(key)['Ref']
        codon = ''.join(list(codon))
        codon_dict[key] = codon

    # Retain variant positions
    orf_df = orf_df.loc[orf_df['is_variant'] == 1, :]

    # Expand df to include triallelics
    x_df = pd.DataFrame(columns=['bp', 'codon_number', 'pos_in_codon', 'snp_count', 'Ref', 'snp', 'ref_AA'])

    for row in list(orf_df.index):
        morsel = orf_df.loc[row, :].to_frame().transpose()
        ref_base = morsel['Ref'].iloc[0]
        allele_list = ['a', 't', 'g', 'c']
        allele_list.remove(ref_base)
        for a in allele_list:
            if morsel[a].iloc[0] != 0:
                crumb = morsel.copy()
                crumb['snp'] = a
                crumb['snp_count'] = morsel[a].iloc[0]
                crumb = crumb[['bp', 'codon_number', 'pos_in_codon', 'snp_count', 'Ref', 'snp', 'ref_AA']]
                x_df = x_df.append(crumb)

    x_df = x_df.reset_index()

# x_df.loc[x_df['codon_number'] == 3606, :]
    def get_sav(x):
        codon_no = x['codon_number']
        codon_pos = x['pos_in_codon']
        snp = x['snp']

        codon = codon_dict[codon_no]
        codon = list(codon)
        codon[codon_pos - 1] = snp
        codon = ''.join(codon)

        sav = str(Seq.Seq(codon).translate())

        return sav


    x_df.insert(loc=7, column='variant', value=x_df.apply(get_sav, axis=1))

    def get_syn(x):
        if x['ref_AA'] == x['variant']:
            return 'S'
        else:
            return 'NS'

    x_df.insert(loc=8, column='NS_S', value=x_df.apply(get_syn, axis=1))
    x_df.insert(loc=9, column='ORF', value=orf_key)

    to_save = to_save.append(x_df)

to_save = to_save[['snp_count', 'Ref', 'bp', 'snp', 'ref_AA', 'codon_number', 'variant', 'pos_in_codon', 'NS_S', 'ORF']]
to_save.columns = ['snp_count', 'ref', 'bp', 'snp', 'ref_AA', 'codon_number', 'var_AA', 'pos_in_codon', 'NS_S', 'ORF']
to_save['ref'] = to_save['ref'].str.upper()
to_save['snp'] = to_save['snp'].str.upper()

to_save.to_csv(dnds / 'parsed_savs_12626_210520.csv', header=True, index=False)

