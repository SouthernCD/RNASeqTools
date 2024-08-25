from yxutil import read_matrix_file, cmd_run
from yxseq import read_gff_file, Gene
from .fpkm import get_TPM_matrix, get_FPKM_matrix


def Count2TMM_main(args):
    count_matrix, sample_list, gene_list = read_matrix_file(args.count_matrix)

    gene_length_dict = {}
    with open(args.gene_length_file, 'r') as f:
        for i in f:
            i.strip()
            g_id, length = i.split()
            gene_length_dict[g_id] = float(length)

    gene_length_list = [gene_length_dict[i] for i in gene_list]

    fpkm = get_FPKM_matrix(count_matrix, gene_length_list)
    tpm = get_TPM_matrix(count_matrix, gene_length_list)

    # output
    with open(args.output_prefix+".fpkm.matrix", 'w') as f:
        f.write("Gene\t" + "\t".join(sample_list) + "\n")

        for g in range(len(gene_list)):
            c_list = []
            for s in range(len(sample_list)):
                c = fpkm[g][s]
                c_list.append(str(c))
            "\t".join(c_list)

            g_id = gene_list[g]
            f.write(g_id + "\t" + "\t".join(c_list) + "\n")

    # output
    with open(args.output_prefix+".tpm.matrix", 'w') as f:
        f.write("Gene\t" + "\t".join(sample_list) + "\n")

        for g in range(len(gene_list)):
            c_list = []
            for s in range(len(sample_list)):
                c = tpm[g][s]
                c_list.append(str(c))
            "\t".join(c_list)

            g_id = gene_list[g]
            f.write(g_id + "\t" + "\t".join(c_list) + "\n")

    tpm_file = args.output_prefix+".tpm.matrix"
    tmm_file = args.output_prefix+".TMM.matrix"

    cmd_string = "run_TMM_scale_matrix.pl --matrix %s > %s" % (
        tpm_file, tmm_file)
    cmd_run(cmd_string)


def get_gene_length_dict(gff_file):
    gff_dict = read_gff_file(gff_file)

    gene_length_dict = {}
    for i in gff_dict:
        for j in gff_dict[i]:
            gene_tmp = Gene(from_gf=gff_dict[i][j])
            gene_max_length = 0
            for mRNA_tmp in gene_tmp.sub_features:
                mRNA_tmp.sgf_len()
                if 'exon' in mRNA_tmp.sgf_len_dir:
                    gene_length = mRNA_tmp.sgf_len_dir['exon']
                else:
                    gene_length = mRNA_tmp.sgf_len_dir['CDS']
                if gene_length > gene_max_length:
                    gene_max_length = gene_length
            gene_length_dict[j] = gene_max_length

    return gene_length_dict


def GetGeneLength_main(args):
    gene_length_dict = get_gene_length_dict(args.gff_file)

    with open(args.gene_length_file, 'w') as f:
        for i in gene_length_dict:
            f.write("%s\t%d\n" % (i, gene_length_dict[i]))
