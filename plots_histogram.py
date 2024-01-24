import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import plotly.graph_objects as go
import plotly.express as px

BGAL_ECOLI_CONTAMINANT = "sp|P00722|BGAL_ECOLI"


def fasta_peptide_dict(raw, proteomes_fasta_list, folder_path):
    """
    creates a dictionary of {fasta file name: peptide's array}
    """
    dictionary = {}
    for fas in proteomes_fasta_list:
        if "peptide.tsv" in os.listdir(f'{folder_path}\\{fas}\\{raw}'):
            df = pd.read_csv(f"{folder_path}\\{fas}\\{raw}\\peptide.tsv", sep="\t")
            col_name = 'Peptide'
            col_array = df[col_name].values
            dictionary[fas] = col_array
    return dictionary


def fasta_human_protein_dict(raw, proteomes_fasta_list, folder_path, contam_notation):
    """
    creates dictionary(fasta_name: list of human proteins (length = number of peptides found))
    """
    dictionary = {}
    for fas in proteomes_fasta_list:
        if "peptide.tsv" in os.listdir(f'{folder_path}\\{fas}\\{raw}'):
            df = pd.read_csv(f"{folder_path}\\{fas}\\{raw}\\peptide.tsv", sep="\t")
            col_name = 'Protein'
            col_array = df[col_name].values
            sp_array = []
            for val in col_array:
                if val.split("_")[-1] in contam_notation:
                    sp_array.append(val)

            dictionary[fas] = sp_array
    return dictionary


def fasta_pandas_df_dict(raw, proteomes_fasta_list, folder_path):
    """
    dict {fasta file name: pandas file of peptide}
    """
    dictionary = {}
    for fas in proteomes_fasta_list:
        if "peptide.tsv" in os.listdir(f'{folder_path}\\{fas}\\{raw}'):
            df = pd.read_csv(f"{folder_path}\\{fas}\\{raw}\\peptide.tsv", sep="\t")
            dictionary[fas] = df
    return dictionary


def fasta_non_human_protein_dict(peptides_dictionary, df_dictionary, contam_notation):
    """
    creates dict of number of peptides.
    values are lists of proteins for each fasta file (key) (length = number of peptides found)
    """
    contaminants_cleaned_dict = {}
    iterable_keys = np.array(list(peptides_dictionary.keys()), dtype=object)
    for fas in iterable_keys:
        # if (len(human_proteins_dictionary[fas]) >= 200) and (len(human_proteins_dictionary[fas]) <= 240):
        contaminants_cleaned_dict[fas] = []
        for val in df_dictionary[fas]["Protein"].values:
            if val.split("_")[-1] not in contam_notation and val != BGAL_ECOLI_CONTAMINANT:
                contaminants_cleaned_dict[fas].append(val)

    return contaminants_cleaned_dict


def create_histogram(folder_path, dictionary, raw, x_title, y_title, histogram_title, n_bins, steps):
    peptides_arrays = np.array(list(dictionary.values()), dtype=object)
    lens = []
    for array in peptides_arrays:
        lens.append(len(array))

    std_value = round(np.std(lens), 3)
    mean = round(np.mean(lens), 3)

    fig = go.Figure()
    fig.add_trace(go.Histogram(x=lens, nbinsx=n_bins, marker_color='blue', opacity=0.7, histnorm='probability density'))

    # Customize layout
    fig.update_layout(
        xaxis_title=x_title,
        yaxis_title=y_title,
        title=f'{histogram_title} mean:{mean} std:{std_value}',
        xaxis=dict(tickvals=np.arange(0, max(lens) + 1, step=steps))
    )

    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='lightgrey')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='lightgrey')

    fig.update_traces(
        marker_line_color='black',
        marker_line_width=1
    )
    fig.write_html(f"{folder_path}\\{histogram_title}_{raw}.html")
    # Show the plot
    fig.show()
