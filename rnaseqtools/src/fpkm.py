import numpy as np


def get_FPKM(gene_reads_dict, gene_length_dict):
    if set(list(gene_reads_dict.keys())) != set(list(gene_length_dict.keys())):
        raise ValueError("input should have same keys!")

    total_mapped_reads = sum([gene_reads_dict[i] for i in gene_reads_dict])
    mm = total_mapped_reads/1000000

    fpkm_dict = {}
    for g_id in gene_reads_dict:
        c = gene_reads_dict[g_id]
        l = gene_length_dict[g_id]
        kl = l/1000
        fpkm = c/kl/mm
        fpkm_dict[g_id] = fpkm

    return fpkm_dict


def get_TPM(gene_reads_dict, gene_length_dict):
    if set(list(gene_reads_dict.keys())) != set(list(gene_length_dict.keys())):
        raise ValueError("input should have same keys!")

    rpk_dict = {}
    for g_id in gene_reads_dict:
        c = gene_reads_dict[g_id]
        l = gene_length_dict[g_id]
        kl = l/1000
        rpk_dict[g_id] = c/kl

    srpkm = sum([rpk_dict[i] for i in rpk_dict])/1000000

    tpm_dict = {}
    for g_id in rpk_dict:
        rpk = rpk_dict[g_id]
        tpm_dict[g_id] = rpk/srpkm

    return tpm_dict


def get_FPKM_matrix(gene_count_matrix, gene_length_list):
    mm = np.array([np.sum(gene_count_matrix, axis=0)])/1000000
    gene_length_array = np.array([gene_length_list])
    rpk = gene_count_matrix / gene_length_array.T * 1000
    fpkm = rpk/mm
    return fpkm


def get_TPM_matrix(gene_count_matrix, gene_length_list):
    gene_length_array = np.array([gene_length_list])
    rpk = gene_count_matrix / gene_length_array.T * 1000
    mm = np.array([np.sum(rpk, axis=0)])/1000000
    tpm = rpk/mm
    return tpm


if __name__ == "__main__":

    """
    gene count matrix
        s1  s2  s3  s4      
    g1  0   0   0   0
    g2  10  10  10  10
    g3  20  18  22  23
    g4  5   0   5   5
    g5  20  22  50  20

    gene length
    g1  1000
    g2  2000
    g3  5000
    g4  3000
    g5  1000
    """

    # for a matrix
    gene_count_matrix = np.array([[0,  0,  0,  0],
                                  [10, 10, 10, 10],
                                  [20, 18, 22, 23],
                                  [5,  0,  5,  5],
                                  [20, 22, 50, 20]])

    gene_length_list = [1000, 2000, 5000, 3000, 1000]

    fpkm_matrix = get_FPKM_matrix(gene_count_matrix, gene_length_list)

    tpm_matrix = get_TPM_matrix(gene_count_matrix, gene_length_list)

    gene_reads_dict = {
        'g1': 0,
        'g2': 10,
        'g3': 20,
        'g4': 5,
        'g5': 20,
    }

    gene_length_dict = {
        'g1': 1000,
        'g2': 2000,
        'g3': 5000,
        'g4': 3000,
        'g5': 1000,
    }
