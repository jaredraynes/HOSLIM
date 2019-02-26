import argparse, textwrap
from argparse import RawDescriptionHelpFormatter
import pandas as pd
import numpy as np
import Bio as bio
import itertools as it
from io import TextIOWrapper
from Bio import ExPASy
from Bio import SeqIO

#       __  __  _____   __      __  __    __
#      / / / / / ___/  / /     /_/ /  \  /  \
#     / /_/ / / /___  / /     __  / /\ \/ /\ \
#    / __  / /___  / / /     / / / /  \__/  \ \
#   / / / / ____/ / / /___  / / / /          \ \
#  /_/ /_/ /_____/ /_____/ /_/ /_/            \_\

#--------------------------------------------------

def fasta_download(accession_number, NresA, NresB):
    '''
    Input the protein accession number, starting and ending residue number and return the amino acid sequence.

    >>> fasta_download("Q16082", 0, 65)
    ('MSGRSVPHAHPATAEYEFANPSRLGEQRFGEGLLPEEILTPTLYHGYYVRPRAAPAGEGSRAGAS', 'Q16082', 'HSPB2', 'Human')
    '''

    #retrieve protein sequence and save as a string
    with ExPASy.get_sprot_raw(accession_number) as handle:
        seq_record = SeqIO.read(handle, "swiss")
        if NresA == -1:
            seq = str(seq_record.seq)
        else:
            seq = str(seq_record.seq[NresA:NresB])
    if 'gene_name' in seq_record.annotations:
        if 'csn1s2' in seq_record.annotations['gene_name'].lower():
            protein_name = 'as2-casein'
        elif 'csn1s1' in seq_record.annotations['gene_name'].lower():
            protein_name = 'as1-casein'
        elif 'csn2' in seq_record.annotations['gene_name'].lower():
            protein_name = 'beta-casein'
        elif 'csn3' in seq_record.annotations['gene_name'].lower():
            protein_name = 'kappa-casein'
        else:
            protein_name = seq_record.annotations['gene_name'][5:(seq_record.annotations['gene_name'].find(';'))]
    else:
        protein_name = seq_record.name
        print('NO GENE NAME FOUND IN UNIPROT RECORD, SEE UNIPROT WEBSITE FOR MORE INFOMATION')
    
    return(seq, seq_record.id, protein_name.translate({ord(i):None for i in '{!@#$}|.:'}), seq_record.annotations["organism"][(seq_record.annotations['organism'].find('(')+1):seq_record.annotations['organism'].find(')')])

def scales(hydropathy_scale):
    '''
    Input either K_D or Guy to choose the hydropathy scale to use

    >>> scale("K_D")
    {'A':0.40, 'R':0.00, 'N':0.11, 'D':0.11, 'C':0.78, 'Q':0.11, 'E':0.11, 'G':0.46, 'H':0.14, 'I':1.00, 'L':0.92, 'K':0.07, 'M':0.71, 'F':0.81, 'P':0.32, 'S':0.41, 'T':0.42, 'W':0.40, 'Y':0.36, 'V':0.97}
    '''

    #selector for which hydropathy scale diction of values for each amino acid the program uses
    if hydropathy_scale == "K_D":
        scale = {'A':0.40, 'R':0.00, 'N':0.11, 'D':0.11, 'C':0.78, 'Q':0.11, 'E':0.11, 'G':0.46, 'H':0.14, 'I':1.00, 'L':0.92, 'K':0.07, 'M':0.71, 'F':0.81, 'P':0.32, 'S':0.41, 'T':0.42, 'W':0.40, 'Y':0.36, 'V':0.97}
    else:
        scale = {'A':0.45, 'R':0.00, 'N':0.35, 'D':0.28, 'C':0.83, 'Q':0.27, 'E':0.24, 'G':0.39, 'H':0.60, 'I':0.75, 'L':0.77, 'K':0.13, 'M':0.87, 'F':1.00, 'P':0.29, 'S':0.34, 'T':0.46, 'W':0.60, 'Y':0.53, 'V':0.79}
    return(scale)

def order(order_scale):
    '''
    Only input is "IDP". Background program to get the order scale into the program

    >>> order("IDP")
    {'A':0.39, 'R':0.27, 'N':0.41, 'D':0.23, 'C':0.65, 'Q':0.18, 'E':0.12, 'G':0.27, 'H':0.42, 'I':0.69, 'L':0.62, 'K':0.10, 'M':0.44, 'F':0.67, 'P':0.00, 'S':0.14, 'T':0.35, 'W':1.00, 'Y':0.72, 'V':0.59}
    '''

    #selector for order scale
    if order_scale == "IDP":
        IDP = {'A':0.39, 'R':0.27, 'N':0.41, 'D':0.23, 'C':0.65, 'Q':0.18, 'E':0.12, 'G':0.27, 'H':0.42, 'I':0.69, 'L':0.62, 'K':0.10, 'M':0.44, 'F':0.67, 'P':0.00, 'S':0.14, 'T':0.35, 'W':1.00, 'Y':0.72, 'V':0.59}
    return(IDP)

def norm_hydro_order(seq, seq_record_id, protein_name, seq_record_organism, scale, IDP, output_file_location, save):
    '''
    Input amino acid sequence, hydropathy_scale library and order scale library, return fo, fp, fn, ho, hp, hn and h, do, dp, dn and d

    >>> norm_hydro("MSGRSVPHAHPATAEYEFANPSRLGEQRFGEGLLPEEILTPTLYHGYYVRPRAAPAGEGSRAGAS", {'A':0.40, 'R':0.00, 'N':0.11, 'D':0.11, 'C':0.78, 'Q':0.11, 'E':0.11, 'G':0.46, 'H':0.14, 'I':1.00, 'L':0.92, 'K':0.07, 'M':0.71, 'F':0.81, 'P':0.32, 'S':0.41, 'T':0.42, 'W':0.40, 'Y':0.36, 'V':0.97}, {'A':0.39, 'R':0.27, 'N':0.41, 'D':0.23, 'C':0.65, 'Q':0.18, 'E':0.12, 'G':0.27, 'H':0.42, 'I':0.69, 'L':0.62, 'K':0.10, 'M':0.44, 'F':0.67, 'P':0.00, 'S':0.14, 'T':0.35, 'W':1.00, 'Y':0.72, 'V':0.59})
    fo = 0.8
    ho = 0.477
    do = 0.363
    fp = 0.092
    hp = 0.0
    dp = 0.27
    fn = 0.108
    hn = 0.11
    dn = 0.12
    h =  0.393
    d =  0.329
    '''

    scale_seq = []
    for a in seq:
        if a in scale.keys():
            scale_seq.append((scale[a]))

    scale_seq_neu = []
    for a in seq:
        if a == 'K':
            continue
        if a == 'R':
            continue
        if a == 'D':
            continue
        if a == 'E':
            continue
        if a in scale.keys():
            scale_seq_neu.append((scale[a]))

    scale_seq_pos = []
    for a in seq:
        if a == 'K':
            scale_seq_pos.append((scale[a]))
        elif a == 'R':
            scale_seq_pos.append((scale[a]))
        else:
            continue

    scale_seq_neg = []
    for a in seq:
        if a == 'D':
            scale_seq_neg.append((scale[a]))
        elif a == 'E':
            scale_seq_neg.append((scale[a]))
        else:
            continue

    IDP_seq = []
    for a in seq:
        if a in IDP.keys():
            IDP_seq.append((IDP[a]))

    IDP_seq_neu = []
    for a in seq:
        if a == 'K':
            continue
        if a == 'R':
            continue
        if a == 'D':
            continue
        if a == 'E':
            continue
        if a in IDP.keys():
            IDP_seq_neu.append((IDP[a]))

    IDP_seq_pos = []
    for a in seq:
        if a == 'K':
            IDP_seq_pos.append((IDP[a]))
        elif a == 'R':
            IDP_seq_pos.append((IDP[a]))
        else:
            continue

    IDP_seq_neg = []
    for a in seq:
        if a == 'D':
            IDP_seq_neg.append((IDP[a]))
        elif a == 'E':
            IDP_seq_neg.append((IDP[a]))
        else:
            continue

    fo = round(len(scale_seq_neu) / len(seq), 3)
    ho = round(sum(scale_seq_neu) / len(scale_seq_neu), 3)
    do = round(sum(IDP_seq_neu) / len(IDP_seq_neu), 3)
    fp = round(len(scale_seq_pos) / len(seq), 3)
    hp = round(sum(scale_seq_pos) / len(scale_seq_pos), 3)
    dp = round(sum(IDP_seq_pos) / len(IDP_seq_pos), 3)
    fn = round(len(scale_seq_neg) / len(seq), 3)
    hn = round(sum(scale_seq_neg) / len(scale_seq_neg), 3)
    dn = round(sum(IDP_seq_neg) / len(IDP_seq_neg), 3)

    h = round((fo*ho) + (fp*hp) + (fn*hn), 3)
    d = round((fo*do) + (fp*dp) + (fn*dn), 3)

    data = [[protein_name, seq_record_id, seq_record_organism, len(seq), fn, fp, h, ho, d, do]]

    df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Nres', 'f-', 'f+', 'h', 'ho', 'd', 'do'])

    if save == 'yes':
        df.to_csv(output_file_location + seq_record_id + '_' + protein_name + '_' + seq_record_organism + '_parameters' + '.csv', index = False)
    print(df.to_string(index=False))
    return(df.to_string(index=False))

def h_slim(seq, seq_record_id, protein_name, seq_record_organism, scale, IDP, motif_length, output_file_location, save_Sw):
    '''
    Input amino acid sequence, hydropathy_scale library and order scale library, return Sw which is a parameter for the hydrophobic and order-promoting short-linear motif (H-SLiM)

    >>> h_slim("MSGRSVPHAHPATAEYEFANPSRLGEQRFGEGLLPEEILTPTLYHGYYVRPRAAPAGEGSRAGAS", {'A':0.40, 'R':0.00, 'N':0.11, 'D':0.11, 'C':0.78, 'Q':0.11, 'E':0.11, 'G':0.46, 'H':0.14, 'I':1.00, 'L':0.92, 'K':0.07, 'M':0.71, 'F':0.81, 'P':0.32, 'S':0.41, 'T':0.42, 'W':0.40, 'Y':0.36, 'V':0.97}, {'A':0.39, 'R':0.27, 'N':0.41, 'D':0.23, 'C':0.65, 'Q':0.18, 'E':0.12, 'G':0.27, 'H':0.42, 'I':0.69, 'L':0.62, 'K':0.10, 'M':0.44, 'F':0.67, 'P':0.00, 'S':0.14, 'T':0.35, 'W':1.00, 'Y':0.72, 'V':0.59})
    For sequence 47 YYV 49 Sw = 6.297
    '''

    #Convert the amino acid sequence to a list of hydropathy values by accessing the scale dictionary
    scale_seq = []
    for a in seq:
        if a in scale.keys():
            scale_seq.append((scale[a]))
    #print(scale_seq)

    #Convert only neutral amino acids to list of neutral hydropathy values by accessing the scale dictionary
    scale_seq_neu = []
    for a in seq:
        if a == 'K':
            continue
        if a == 'R':
            continue
        if a == 'D':
            continue
        if a == 'E':
            continue
        if a in scale.keys():
            scale_seq_neu.append((scale[a]))

    #Convert the amino acid sequence to a list of order propensity values by accessing the order dictionary
    IDP_seq = []
    for a in seq:
        if a in IDP.keys():
            IDP_seq.append((IDP[a]))
    #print(IDP_seq)

    #Convert only neutral amino acids to a list of neutral order values by accessing the order dictionary
    IDP_seq_neu = []
    for a in seq:
        if a == 'K':
            continue
        if a == 'R':
            continue
        if a == 'D':
            continue
        if a == 'E':
            continue
        if a in IDP.keys():
            IDP_seq_neu.append((IDP[a]))


    #Go through the amino acid sequence and return Sw for H-SLiMs
    data = []
    start = 0
    s = 3
    e = motif_length
    motif = ['', '', '', 0, 0, 0, '', 0]

    for a, b, c in it.zip_longest(seq, scale_seq, IDP_seq):
        if e < (len(seq) + motif_length): #struggling with range to get to go right up to end of sequence!
            for i in range(s, e):
                seq1 = ''
                h = []
                d = []
                seq1 += seq[start:i]
                #screen out positively charged amino acids
                if 'K' in seq1:
                    continue
                if 'R' in seq1:
                    continue
                #screen out negatively charged amino acids
                if 'D' in seq1:
                    continue
                if 'E' in seq1:
                    continue
                #screen out all amino acids that have an hd less than that for Arg
                if 'P' in seq1:
                    continue
                if 'G' in seq1:
                    continue
                if 'N' in seq1:
                    continue
                if 'S' in seq1:
                    continue
                if 'Q' in seq1:
                    continue
                if 'T' in seq1:
                    continue
                #this selector determines whether to bypass 'H' if it is in the K_D scale or keep it in for the Hug scale
                if scale['H'] == 0.14:
                    if 'H' in seq1:
                        continue
                h += scale_seq[start:i]
                d += IDP_seq[start:i]
                ho = np.mean(scale_seq_neu)
                do = np.mean(IDP_seq_neu)
                num = 0

                for i in range(len(h)):
                    num += (h[i] * d[i])
                    Sw = round(num / (ho * do), 1)

                if len(seq1) > len(motif[6]):
                    motif = [protein_name, seq_record_id, seq_record_organism, Sw, len(seq1), (s - 2), seq1, ((s - 3) + len(seq1))]
                    if len(data) == 0:
                        data.append(motif)
                else:
                    data.append(motif)
                    motif = [protein_name, seq_record_id, seq_record_organism, Sw, len(seq1), (s - 2), seq1, ((s - 3) + len(seq1))]
                    if len(seq1) == len(motif[6]):
                        data.append(motif)
                        motif = [protein_name, seq_record_id, seq_record_organism, Sw, len(seq1), (s - 2), seq1, ((s - 3) + len(seq1))]

        s += 1
        e += 1
        start += 1

    df = pd.DataFrame(data, columns = ['Protein name', 'Uniprot Code', 'Organism', 'Sw', 'Nres', 'Start', 'Sequence', 'End'])
    df.drop_duplicates(subset = 'Start', keep = 'last', inplace = True) #could make this an argument
    df.drop_duplicates(subset = 'End', keep = 'first', inplace = True) #could make this an argument
    df = df[df.Nres >= 3]
    df[['Protein name', 'Uniprot Code', 'Organism']] = df[['Protein name', 'Uniprot Code', 'Organism']].where(df[['Protein name', 'Uniprot Code', 'Organism']].apply(lambda x: x != x.shift()), '')
    
    print(df.to_string(index=False))
    
    if save_Sw == 'yes':
        df.to_csv(output_file_location + seq_record_id + '_' + protein_name + '_' + seq_record_organism + '_Sw' + '.csv', index = False)
    
    return(df.to_string(index=False))

def main():
    input = argparse.ArgumentParser(prog='HSLiM', formatter_class = argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
         -------------------------------------------------------
         For the investigation of protein sequences for H-SLiMs.
         See the published paper: Holt, Raynes and Carver (2019)'''))
    input.add_argument("accession_number", help = 'Input your protein accession number such as from UniProt e.g. "P02662".', type = str)
    input.add_argument("-a", "--NresA", help = 'Input the start amino acid residue. If not used the full protein sequence will be used.', default = -1, type = int)
    input.add_argument("-b", "--NresB", help = 'Input the end amino acid residue. If not used the full protein sequence will be used.', default = -1, type = int)
    input.add_argument("-m", "--motif_length", help = 'Input the maximum length of the H-SLiM you are looking for. The default maximum is 11 amino acids. Python indexing means add 1 to the sequence length you would like.', default = 12, type = int)
    input.add_argument("-s", "--hydropathy_scale", help = 'Choose to use either the Kyte & Doolittle or Guy hydropathy scale for normalized amino acid residue hydropathy. The default is K_D.', default = "K_D", type = str)
    input.add_argument("-p", "--save_parameters", help = 'Type "yes" if you would like to save the parameters data frame. Make sure to input your absolute file location in the output_file_location option.', default = "no_save", type = str,)
    input.add_argument("-q", "--save_Sw", help = 'Type "yes" if you would like to save the H-SLiMs data frame. Make sure to input your absolute file location in the output_file_location option.', default = "no_save", type = str,)
    input.add_argument("-r", "--output_file_location", help = 'Type your absolute file location to save a .csv of the parameters table. The default will save it to the same file location', default = "", type = str)
    args = input.parse_args()
    seq = fasta_download(args.accession_number, args.NresA, args.NresB)
    scale = scales(args.hydropathy_scale)
    IDP = order("IDP")
    norm_hydro_order(seq[0], seq[1], seq[2], seq[3], scale, IDP, args.output_file_location, args.save_parameters)
    h_slim(seq[0], seq[1], seq[2], seq[3], scale, IDP, args.motif_length, args.output_file_location, args.save_Sw)

if __name__ == "__main__":
	main()
