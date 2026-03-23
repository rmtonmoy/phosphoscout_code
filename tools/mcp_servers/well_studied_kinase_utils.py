

from pathlib import Path
import pandas as pd
from Bio import SeqIO
import random

PROJECT_ROOT = Path(__file__).resolve().parents[2]

cantley_kinase_path = str(PROJECT_ROOT / 'data' / 'kinase_data' / '41586_2022_5575_MOESM3_ESM.xlsx')

curated_dom_seq_path = str(PROJECT_ROOT / 'data' / 'kinase_data' / 'pkfold_hs_curated.fa')

def get_well_studied_kinase_ids():
    data = pd.read_excel(cantley_kinase_path, sheet_name = 1, engine='openpyxl')
    return data['Uniprot id'].to_list()

def get_well_studied_human_kinase_ids():
    ids = get_well_studied_kinase_ids()
    cancelled = ["O14874", "Q15118", "Q16654", "Q3UTQ8", "E0W1I1", "Q63285"]
    return [id for id in ids if id not in cancelled]



def get_domain_seq_from_id(uniprot_id, fasta_file = curated_dom_seq_path):

    domain_seqs = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == uniprot_id or record.id == f"{uniprot_id}|1" or record.id == f"{uniprot_id}|2":
            domain_seqs.append(str(record.seq))

    if(len(domain_seqs) == 0):
        raise AssertionError("Domain sequence not found")

    return domain_seqs


if __name__ == "__main__":
    #print((get_well_studied_human_kinase_ids()))
    print(get_domain_seq_from_id("Q2M2I8"))
