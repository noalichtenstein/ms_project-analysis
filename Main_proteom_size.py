import os.path
import numpy as np
import pandas as pd
import plots_histogram
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import kaleido
import plotly


def convert_file_to_list(folder, organism_name):
    fasta_list = []
    with open(f"{folder}\\{organism_name}.txt", 'r') as file:
        line = file.readline()

        while line:
            fasta_list.append(line.strip())
            line = file.readline()

        file.close()
    return fasta_list


def convert_list_to_file(folder, file_name, points_list):
    with open(f"{folder}\\{file_name}.txt", 'w') as file:
        for point in points_list:
            file.write(f"{point[0]} {point[1]}\n")
        file.close()


def calculate_proteomes_sizes(all_proteomes_folder, fas):
    """
    return number of protein in fasta file
    """
    prot_num = 0
    with open(f"{all_proteomes_folder}\\{fas}\\{fas}.fasta") as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                prot_num += 1
    return prot_num


if __name__ == '__main__':
    """
    for each fasta:
    find the number of proteins in the file
    array of fastas sizes
    array of number of peptides
    """
    raw_list = [
        # "2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2",
        # "2022_05_25_07_PreNatal3_Sample2_S1_SCX_p2",
        "2022_05_25_05_PreNatal3_Sample3_S1_SCX_p2",
        # "2022_05_25_03_PreNatal3_Sample4_S1_SCX_p2",
        # "2022_05_25_01_PreNatal3_Sample5_S1_SCX_p2",
    ]
    all_proteomes_folder = r'C:\Users\owner\Microbiome\Analyses\Klebsiella_aerogenes_analasys_human_ms'
    organism_name = 'Klebsiella+aerogenes'
    results_folder = (f"C:\\Users\\owner\\Microbiome\\Analyses\\proteome_size_vs_peptides\\"
                      f"{organism_name.replace('+', '_')}")
    # os.makedirs(results_folder)
    contam_notation = ['HUMAN', 'PIG', 'BOVIN', 'SHEEP', 'YEAST']
    proteomes_fasta_list = convert_file_to_list(all_proteomes_folder, organism_name)

    proteomes_sizes = []
    number_of_peptides = []
    proteomes_sizes_std = []
    number_of_peptides_std = []
    ref_prot_peptides = 0
    ref_prot_proteins = 0

    for chosen_raw in raw_list:
        peptides_dictionary = plots_histogram.fasta_peptide_dict(chosen_raw, proteomes_fasta_list,
                                                                 all_proteomes_folder)
        human_proteins_dictionary = plots_histogram.fasta_human_protein_dict(chosen_raw, proteomes_fasta_list,
                                                                             all_proteomes_folder, contam_notation)
        df_dictionary = plots_histogram.fasta_pandas_df_dict(chosen_raw, proteomes_fasta_list, all_proteomes_folder)
        non_human_dict = plots_histogram.fasta_non_human_protein_dict(peptides_dictionary, df_dictionary, contam_notation)

        for dictionary_type in [human_proteins_dictionary, non_human_dict]:
            points_in_plot = []

            for fas in proteomes_fasta_list:
                if fas in dictionary_type.keys():
                    proteomes_sizes.append(calculate_proteomes_sizes(all_proteomes_folder, fas))
                    number_of_peptides.append(len(dictionary_type[fas]))
                    if fas == "uniprotkb_proteome_UP000000818_2023-11-12":
                        ref_prot_peptides = len(dictionary_type[fas])
                        ref_prot_proteins = calculate_proteomes_sizes(all_proteomes_folder, fas)
                    points_in_plot.append(
                        [calculate_proteomes_sizes(all_proteomes_folder, fas), len(dictionary_type[fas])])
                else:
                    print(f"not found {fas}")

            std_value = round(np.std(number_of_peptides), 3)
            mean = round(np.mean(number_of_peptides), 3)

            # removing outliers
            for fas in proteomes_fasta_list:
                if (fas in dictionary_type.keys() and
                        (mean - 3 * std_value) < len(dictionary_type[fas]) < (mean + 3 * std_value)):
                    proteomes_sizes_std.append(calculate_proteomes_sizes(all_proteomes_folder, fas))
                    number_of_peptides_std.append(len(dictionary_type[fas]))

            final_std_value = round(np.std(number_of_peptides_std), 3)
            final_mean = round(np.mean(number_of_peptides_std), 3)
            if final_mean != 0:
                coef_of_variation = final_std_value / final_mean * 100
            else:
                coef_of_variation = 0

            if dictionary_type == human_proteins_dictionary:
                name_human_or_bacteria = 'human'
            else:
                name_human_or_bacteria = 'bacteria'
            print(f"Coefficient of variation in {name_human_or_bacteria}:", coef_of_variation)

            convert_list_to_file(results_folder, f'points_{name_human_or_bacteria}_{chosen_raw}', points_in_plot)
            # Create a scatter plot with a linear regression line
            raw_number = chosen_raw.split('_')[5]
            scatter_plot = px.scatter(x=proteomes_sizes, y=number_of_peptides, trendline='ols',
                                      labels={'x': 'Number of proteins in fasta file', 'y': 'Number of peptides'},
                                      trendline_color_override="blue",
                                      height=1000,
                                      width=1000
                                      )
            scatter_plot.update_traces(marker=dict(color="black", size=10))
            scatter_plot.update_layout(plot_bgcolor="white",
                title={
                    'text': (
                        f'Connection Between Number Of Peptides And Number Of Proteins <br>in'
                        f' {organism_name.replace("+", " ")}, sample {raw_number[-1]} in '
                        f'{name_human_or_bacteria} proteins'
                        # f'<br>std:{final_std_value} mean:{final_mean}\n'
                        # f'reference proteome: proteins-{ref_prot_proteins}, peptides-{ref_prot_peptides}'

                    ),
                    'font': {'size': 26},  # Font properties for the main title
                })
            scatter_plot.update_layout(
                xaxis=dict(
                    title='<b>Number of proteins in fasta file</b>',
                    tickfont=dict(size=24),
                    title_font=dict(size=24),
                    range=[1, max(proteomes_sizes) + 100],
                    linecolor='black',  # Set the color of the x-axis line
                    linewidth=2,  # Set the width of the x-axis line
                ),

                yaxis=dict(
                    title='<b>Number of peptides</b>',
                    title_font=dict(size=24),
                    tickfont=dict(size=24),
                    range=[-2, max(number_of_peptides) + 100],
                    linecolor='black',  # Set the color of the x-axis line
                    linewidth=2,  # Set the width of the x-axis line
                ),
                margin = dict(t=100)
            )
            # scatter_plot.show()
            scatter_plot.write_image(
                f"{results_folder}\\{organism_name.replace('+', '_')}_{raw_number}_{name_human_or_bacteria}.png")
            scatter_plot.write_html(
                f"{results_folder}\\{organism_name.replace('+', '_')}_{raw_number}_{name_human_or_bacteria}.html")

            points_in_plot = []
            number_of_peptides = []
            proteomes_sizes = []
            proteomes_sizes_std = []
            number_of_peptides_std = []
