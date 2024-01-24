import os.path
import numpy as np
import pandas as pd
import plots_histogram
import matplotlib.pyplot as plt


def convert_file_to_list(folder, organism_name):
    fasta_list = []
    with open(f"{folder}\\{organism_name}.txt", 'r') as file:
        line = file.readline()

        while line:
            fasta_list.append(line.strip())
            line = file.readline()

        file.close()
    return fasta_list


if __name__ == '__main__':
    results_folder = r"C:\Users\owner\Microbiome\Analyses\compare_100_human\plots\Klebsiella_michiganensis"
    raw_list = [
        "2022_05_25_09_PreNatal3_Sample1_S1_SCX_p2", "2022_05_25_07_PreNatal3_Sample2_S1_SCX_p2",
        "2022_05_25_05_PreNatal3_Sample3_S1_SCX_p2",
        "2022_05_25_03_PreNatal3_Sample4_S1_SCX_p2",
        "2022_05_25_01_PreNatal3_Sample5_S1_SCX_p2",
    ]
    list_of_fasta_directory = r'C:\Users\owner\Microbiome\Analyses\compare_100_human\Klebsiella_michiganensis'
    organism_name = 'Klebsiella+michiganensis'
    contam_notation = ['HUMAN', 'PIG', 'BOVIN', 'SHEEP', 'YEAST']
    proteomes_fasta_list = convert_file_to_list(list_of_fasta_directory, organism_name)

    for chosen_raw in raw_list:
        peptides_dictionary = plots_histogram.fasta_peptide_dict(chosen_raw, proteomes_fasta_list,
                                                                 list_of_fasta_directory)

        # this plot shows the distribution of peptides length
        plots_histogram.create_histogram(results_folder, peptides_dictionary, chosen_raw, 'Peptides Number',
                                         'Proteomes', f"All Peptides Histogram {chosen_raw}", 50, steps=100)

        # this plot shows the distribution of peptides of human
        human_proteins_dictionary = plots_histogram.fasta_human_protein_dict(chosen_raw, proteomes_fasta_list,
                                                                             list_of_fasta_directory, contam_notation)
        plots_histogram.create_histogram(results_folder, human_proteins_dictionary, chosen_raw, 'Human Peptides number',
                                         'Proteomes',
                                         f"Human Peptides Histogram {chosen_raw}", 50, steps=10)

        # this plot shows the distribution of non human peptides (only from peptides files with sp proteins in range 200-240)
        df_dictionary = plots_histogram.fasta_pandas_df_dict(chosen_raw, proteomes_fasta_list, list_of_fasta_directory)
        non_human_dict = plots_histogram.fasta_non_human_protein_dict(peptides_dictionary, df_dictionary, contam_notation)
        plots_histogram.create_histogram(results_folder, non_human_dict, chosen_raw, 'Non Human Peptides number',
                                         'Proteomes', f"Non Human Peptides Histogram {chosen_raw}", 50, steps=100)

        # reference proteome for: Enterococcus faecalis
        # print(f"\nfor raw file: {chosen_raw} in the reference proteome")
        # print("human:", len(human_proteins_dictionary["uniprotkb_proteome_UP000001415_2023-11-20"]))
        # print("bacteria:", len(non_human_dict["uniprotkb_proteome_UP000001415_2023-11-20"]))

        # reference proteome for: Clostridium_perfringens
        # print(f"\nfor raw file: {chosen_raw} in the reference proteome")
        # print("human:", len(human_proteins_dictionary["uniprotkb_proteome_UP000000818_2023-11-12"]))
        # print("bacteria:", len(non_human_dict["uniprotkb_proteome_UP000000818_2023-11-12"]))

        # # proteomes of Klebsiella+aerogenes
        # print(f"\nfor raw file: {chosen_raw} in the proteome: uniparc_upid_UP000011170_2023-11-14")
        # print("human:", len(human_proteins_dictionary["uniparc_upid_UP000011170_2023-11-14"]))
        # print("bacteria:", len(non_human_dict["uniparc_upid_UP000011170_2023-11-14"]))
        #
        # print(f"\nfor raw file: {chosen_raw} in the proteome: uniparc_upid_UP000033582_2023-11-14")
        # print("human:", len(human_proteins_dictionary["uniparc_upid_UP000033582_2023-11-14"]))
        # print("bacteria:", len(non_human_dict["uniparc_upid_UP000033582_2023-11-14"]))

        # proteomes of Klebsiella+michiganensis
        print(f"\nfor raw file: {chosen_raw} in the proteome: uniprotkb_proteome_UP000020202_2023_11_22")
        print("human:", len(human_proteins_dictionary["uniprotkb_proteome_UP000020202_2023_11_22"]))
        print("bacteria:", len(non_human_dict["uniprotkb_proteome_UP000020202_2023_11_22"]))

        print(f"\nfor raw file: {chosen_raw} in the proteome: uniprotkb_proteome_UP000077713_2023_11_22")
        print("human:", len(human_proteins_dictionary["uniprotkb_proteome_UP000077713_2023_11_22"]))
        print("bacteria:", len(non_human_dict["uniprotkb_proteome_UP000077713_2023_11_22"]))
