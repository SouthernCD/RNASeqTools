from yxutil import cmd_run, mkdir, ln_file, have_file, tsv_file_dict_parse_big
from yxalign import read_psl_file, outfmt6_read_big, hit_CIP, hit_CALP
from yxseq import read_fasta_by_faidx, read_gff_file, Genome
from interlap import InterLap


def run_blat(ref_genome_fasta, de_novo_fasta, work_dir):
    if not have_file(work_dir + "/blat.psl"):
        ln_file(ref_genome_fasta, work_dir + "/ref.fa")
        ln_file(de_novo_fasta, work_dir + "/de_novo.fa")
        blat_build_ref_cmd_string = "blat ref.fa /dev/null /dev/null -tileSize=11 -makeOoc=11.ooc -repMatch=1024"
        cmd_run(blat_build_ref_cmd_string, work_dir)
        blat_run_cmd_string = "blat ref.fa de_novo.fa -q=rna -dots=100 -maxIntron=500000 -out=psl -ooc=11.ooc blat.psl"
        cmd_run(blat_run_cmd_string, work_dir)
    return work_dir + "/blat.psl"


def read_gene_trans_map(gene_trans_map_file):
    all_tg2tt_dict = {}
    for tmp_id, tmp_dict in tsv_file_dict_parse_big(gene_trans_map_file, fieldnames=['tg', 'tt'], key_col='tt'):
        tg = tmp_dict['tg']
        all_tg2tt_dict.setdefault(tg, []).append(tmp_dict['tt'])
    return all_tg2tt_dict


def read_rsem_gene_results(rsem_gene_results):
    gene_count_dict = {}
    for tmp_id, tmp_dict in tsv_file_dict_parse_big(rsem_gene_results, key_col='gene_id'):
        gene_count_dict[tmp_id] = round(
            float(tmp_dict['expected_count']))
    return gene_count_dict


def get_unigene_hits_dict_from_blat_results(blat_results_file, tt2tg_dict, coverage_threshold=0.6):
    psl_gf_dict = read_psl_file(blat_results_file)

    unigene_hits_dict = {}
    for hit_id in psl_gf_dict:
        hit_gf = psl_gf_dict[hit_id]
        coverage = hit_gf.qualifiers['match']/hit_gf.qualifiers['q_size']
        if coverage >= coverage_threshold:
            iso_id = hit_gf.qualifiers['q_name'].split(".")[0]
            gene_id = tt2tg_dict[iso_id]
            unigene_hits_dict.setdefault(gene_id, []).append(hit_gf)

    return unigene_hits_dict


def hit_sum_hsp_aln(hit):
    aln_list = []
    for hsp in hit.hsp:
        aln_list.append(hsp.Hsp_align_len)
    return sum(aln_list)


def get_ref_gene_to_unigene_map_by_blast(ref_genome_fasta, ref_genome_gff, de_novo_fasta, tt2tg_dict, work_dir, ref_cDNA_fasta=None):
    # read de_novo_fasta
    seq_dict = read_fasta_by_faidx(de_novo_fasta)
    len_dict = {i: seq_dict[i].len() for i in seq_dict}

    # read ref genome and get cDNA sequences

    cDNA_fasta_file = work_dir + "/ref_cDNA.fa"
    if not have_file(cDNA_fasta_file):
        if ref_cDNA_fasta:
            ln_file(ref_cDNA_fasta, cDNA_fasta_file)
        else:
            ref_genome = Genome(genome_file=ref_genome_fasta,
                                gff_file=ref_genome_gff)
            ref_genome.genome_feature_parse()
            ref_genome.build_gene_sequence()

            with open(cDNA_fasta_file, 'w') as OUT:
                for gene_id in ref_genome.feature_dict['gene']:
                    gf = ref_genome.feature_dict['gene'][gene_id]
                    OUT.write(">{}\n{}\n".format(
                        gene_id, gf.model_cDNA_seq))

    ref_gene_dict = read_fasta_by_faidx(cDNA_fasta_file)

    # blast
    bls_out_file = work_dir + "/trans_vs_cds.bls"
    if not have_file(bls_out_file):
        cmd_string = "makeblastdb -in %s -dbtype nucl" % cDNA_fasta_file
        cmd_run(cmd_string, cwd=work_dir)
        cmd_string = "blastn -query %s -db %s -out %s -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -num_threads 1" % (
            de_novo_fasta, cDNA_fasta_file, bls_out_file)
        cmd_run(cmd_string, cwd=work_dir)

    # read blast
    t2g_dict = {}
    for query in outfmt6_read_big(bls_out_file):
        if query.qID not in tt2tg_dict:
            continue

        q_id = tt2tg_dict[query.qID]

        t2g_dict[q_id] = []
        query.qLen = len_dict[query.qID]
        for hit in query.hit:
            cip = hit_CIP(hit)
            calp = hit_CALP(hit)
            sum_aln = hit_sum_hsp_aln(hit)

            if cip > 0.6 and calp > 0.5 and sum_aln > 50:
                t2g_dict[q_id].append((hit.Hit_id, (cip * calp * sum_aln)))

        if len(t2g_dict[q_id]) != 0:
            t2g_dict[q_id] = sorted(
                t2g_dict[q_id], key=lambda x: x[1], reverse=True)[0][0]
        else:
            del t2g_dict[q_id]

    g2t_dict = {}
    for t in t2g_dict:
        g = t2g_dict[t]
        if g not in g2t_dict:
            g2t_dict[g] = []
        g2t_dict[g].append(t)

    for g_id in ref_gene_dict:
        if g_id not in g2t_dict:
            g2t_dict[g_id] = []

    return g2t_dict


def get_ref_gene_to_unigene_map_by_blat(ref_genome_fasta, ref_genome_gff, de_novo_fasta, tt2tg_dict, work_dir, coverage_threshold=0.6):

    # read ref genome
    gene_gff_dict = read_gff_file(ref_genome_gff)
    gene_gff_dict = gene_gff_dict['gene']

    # run blat
    blat_results_file = run_blat(ref_genome_fasta, de_novo_fasta, work_dir)

    # read blat results
    unigene_hits_dict = get_unigene_hits_dict_from_blat_results(
        blat_results_file, tt2tg_dict, coverage_threshold=coverage_threshold)
    unigene_top_hits_dict = {}
    for gene_id in unigene_hits_dict:
        unigene_top_hits_dict[gene_id] = sorted(
            unigene_hits_dict[gene_id], key=lambda x: x.qualifiers['match'], reverse=True)[0]

    # get gene to unigene map
    gene_interlap_dict = {}
    for gene_id in gene_gff_dict:
        gene_gf = gene_gff_dict[gene_id]
        gene_interlap_dict.setdefault(gene_gf.chr_id, InterLap())
        gene_interlap_dict[gene_gf.chr_id].add(
            (gene_gf.start, gene_gf.end, gene_id))

    ug2g_dict = {}
    for ug_id in unigene_top_hits_dict:
        ug_gf = unigene_top_hits_dict[ug_id]
        ug2g_dict[ug_id] = []
        if ug_gf.chr_id in gene_interlap_dict:
            for g_s, g_e, g_id in gene_interlap_dict[ug_gf.chr_id].find((ug_gf.start, ug_gf.end)):
                ug2g_dict[ug_id].append(g_id)

    g2ug_dict = {}
    for ug_id in ug2g_dict:
        for g_id in ug2g_dict[ug_id]:
            g2ug_dict.setdefault(g_id, []).append(ug_id)

    for g_id in gene_gff_dict:
        if not g_id in g2ug_dict:
            g2ug_dict[g_id] = []

    return g2ug_dict


def DenovoCount2RefCount_main(args):

    mkdir(args.work_dir, True)

    # read de novo info
    tg2tt_dict = read_gene_trans_map(args.de_novo_trans_gene_map)
    tt2tg_dict = {}
    for tg in tg2tt_dict:
        for tt in tg2tt_dict[tg]:
            tt2tg_dict[tt] = tg

    # read rsem gene results
    unigene_count_dict = read_rsem_gene_results(args.rsem_gene_results)

    # get gene to unigene map
    if args.map_program == 'blat':
        g2ug_dict = get_ref_gene_to_unigene_map_by_blat(
            args.ref_genome_fasta, args.ref_genome_gff, args.de_novo_fasta, tt2tg_dict, args.work_dir)
    elif args.map_program == 'blast':
        g2ug_dict = get_ref_gene_to_unigene_map_by_blast(
            args.ref_genome_fasta, args.ref_genome_gff, args.de_novo_fasta, tt2tg_dict, args.work_dir, ref_cDNA_fasta=args.ref_cDNA_fasta)

    ug2g_dict = {}
    for g_id in g2ug_dict:
        for ug_id in g2ug_dict[g_id]:
            ug2g_dict.setdefault(ug_id, []).append(g_id)

    # get gene count
    gene_count_dict = {}
    for g_id in g2ug_dict:
        gene_count_dict[g_id] = 0
        for ug_id in g2ug_dict[g_id]:
            gene_count_dict[g_id] += unigene_count_dict[ug_id]

    for i in g2ug_dict:
        if not i in gene_count_dict:
            gene_count_dict[i] = 0

    gene_count_file = args.work_dir + "/gene_count.tsv"
    with open(gene_count_file, 'w') as OUT:
        OUT.write("gene_id\tcount\n")
        for g_id in gene_count_dict:
            OUT.write("{}\t{}\n".format(g_id, gene_count_dict[g_id]))

    # get statistics
    unigene_number = len(
        [i for i in unigene_count_dict if unigene_count_dict[i] > 0])
    unigene_count = sum([unigene_count_dict[i] for i in unigene_count_dict])
    unigene_number_with_hits = len(
        [i for i in ug2g_dict if len(ug2g_dict[i]) > 0])
    unigene_count_with_hits = sum([unigene_count_dict[i] for i in ug2g_dict if len(
        ug2g_dict[i]) > 0 and i in unigene_count_dict and unigene_count_dict[i] > 0])
    gene_number = len(g2ug_dict)
    gene_number_with_count = len(
        [i for i in gene_count_dict if gene_count_dict[i] > 0])
    gene_count = sum([gene_count_dict[i] for i in gene_count_dict])

    with open(args.work_dir + "/statistics.txt", 'w') as OUT:
        OUT.write("unigene_number\t{}\n".format(unigene_number))
        OUT.write("unigene_count\t{}\n".format(unigene_count))
        OUT.write("unigene_number_with_hits\t{}\n".format(
            unigene_number_with_hits))
        OUT.write("unigene_count_with_hits\t{}\n".format(
            unigene_count_with_hits))
        OUT.write("gene_number\t{}\n".format(gene_number))
        OUT.write("gene_number_with_count\t{}\n".format(
            gene_number_with_count))
        OUT.write("gene_count\t{}\n".format(gene_count))

# def DenovoCount2RefCount_main(args):

#     args.work_dir = os.path.abspath(args.work_dir)

#     # build env
#     mkdir(args.work_dir, False)

#     gene_model_cds_file = args.work_dir + "/gene_model_cds.fna"
#     tran_seq_file = args.work_dir + "/trans.fna"
#     gene_trans_map_file = args.work_dir + "/gene_trans_map"
#     os.symlink(os.path.abspath(args.gene_model_cds_file), gene_model_cds_file)
#     os.symlink(os.path.abspath(args.tran_fasta_file), tran_seq_file)
#     os.symlink(os.path.abspath(args.gene_trans_map), gene_trans_map_file)

#     if args.tran_count_fof:
#         tran_count_fof = args.work_dir + "/trans.counts.fof"
#         os.symlink(os.path.abspath(args.tran_count_fof), tran_count_fof)
#     else:
#         tran_count_fof = None

#     if args.tran_count_matrix:
#         tran_count_matrix = args.work_dir + "/tran_count_matrix"
#         os.symlink(os.path.abspath(args.tran_count_matrix), tran_count_matrix)
#     else:
#         tran_count_matrix = None

#     # blast
#     cmd_string = "makeblastdb -in %s -dbtype nucl" % gene_model_cds_file
#     cmd_run(cmd_string, cwd=args.work_dir)

#     bls_out_file = args.work_dir + "/trans_vs_cds.bls"
#     cmd_string = "blastn -query %s -db %s -out %s -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -num_threads %d" % (
#         tran_seq_file, gene_model_cds_file, bls_out_file, args.threads)
#     cmd_run(cmd_string, cwd=args.work_dir)

#     # id map
#     seq_dict = read_fasta_by_faidx(tran_seq_file)
#     len_dict = {i: seq_dict[i].len() for i in seq_dict}

#     all_tg2tt_dict = {}
#     for tmp_id, tmp_dict in tsv_file_dict_parse_big(gene_trans_map_file, fieldnames=['tg', 'tt'], key_col='tt'):
#         tg = tmp_dict['tg']
#         all_tg2tt_dict.setdefault(tg, []).append(tmp_dict['tt'])

#     tg2tt_dict = {}
#     for tg in all_tg2tt_dict:
#         tg2tt_dict[tg] = sorted(
#             all_tg2tt_dict[tg], key=lambda x: len_dict[x], reverse=True)[0]

#     tt2tg_dict = {tg2tt_dict[tg]: tg for tg in tg2tt_dict}

#     # read blast
#     t2g_dict = {}
#     for query in outfmt6_read_big(bls_out_file):
#         if query.qID not in tt2tg_dict:
#             continue

#         q_id = tt2tg_dict[query.qID]

#         t2g_dict[q_id] = []
#         query.qLen = len_dict[query.qID]
#         for hit in query.hit:
#             cip = hit_CIP(hit)
#             calp = hit_CALP(hit)
#             sum_aln = hit_sum_hsp_aln(hit)

#             if cip > 0.6 and calp > 0.5 and sum_aln > 50:
#                 t2g_dict[q_id].append((hit.Hit_id, (cip * calp * sum_aln)))

#         if len(t2g_dict[q_id]) != 0:
#             t2g_dict[q_id] = sorted(
#                 t2g_dict[q_id], key=lambda x: x[1], reverse=True)[0][0]
#         else:
#             del t2g_dict[q_id]

#     g2t_dict = {}
#     for t in t2g_dict:
#         g = t2g_dict[t]
#         if g not in g2t_dict:
#             g2t_dict[g] = []
#         g2t_dict[g].append(t)

#     # read tran count
#     if tran_count_fof:
#         rsem_file_dict = {}
#         sample_list = []
#         with open(tran_count_fof, 'r') as f:
#             for l in f:
#                 l.strip()
#                 s, r_file = l.split()
#                 rsem_file_dict[s] = r_file
#                 sample_list.append(s)

#         tran_count_dict = {}
#         for s in rsem_file_dict:
#             r_file = rsem_file_dict[s]
#             for tmp_id, tmp_dict in tsv_file_dict_parse_big(r_file, key_col='gene_id'):
#                 tran_count_dict.setdefault(
#                     tmp_id, {sd: 0 for sd in sample_list})
#                 tran_count_dict[tmp_id][s] = round(
#                     float(tmp_dict['expected_count']))

#     elif tran_count_matrix:
#         ttran_count_dict = {}
#         for tmp_id, tmp_dict in tsv_file_dict_parse_big(tran_count_matrix, key_col=''):
#             sample_list = sorted([i for i in tmp_dict if i != ''])
#             ttran_count_dict[tmp_id] = {
#                 s: round(float(tmp_dict[s])) for s in sample_list}

#         tran_count_dict = {}
#         for tg in all_tg2tt_dict:
#             tran_count_dict[tg] = {s: 0 for s in sample_list}
#             for tt in all_tg2tt_dict[tg]:
#                 if tt in ttran_count_dict:
#                     for s in ttran_count_dict[tt]:
#                         tran_count_dict[tg][s] += ttran_count_dict[tt][s]

#     # map count to gene
#     gene_count_dict = {}
#     for g in g2t_dict:
#         gene_count_dict[g] = {i: 0 for i in sample_list}
#         for t in g2t_dict[g]:
#             sample_count_dict = tran_count_dict[t]
#             for s in sample_count_dict:
#                 gene_count_dict[g][s] += tran_count_dict[t][s]

#     # output
#     with open(args.output_prefix+".counts.matrix", 'w') as f:
#         f.write("Gene\t" + "\t".join(sample_list) + "\n")

#         for g in gene_count_dict:
#             c_list = []
#             for s in sample_list:
#                 c = gene_count_dict[g][s]
#                 c_list.append(str(c))
#             "\t".join(c_list)

#             f.write(g + "\t" + "\t".join(c_list) + "\n")

#     # get TPM
#     gene_list = sorted(list(gene_count_dict.keys()))
#     count_matrix = []
#     for g in gene_list:
#         g_r = []
#         for s in sample_list:
#             g_r.append(gene_count_dict[g][s])
#         count_matrix.append(g_r)
#     count_matrix = np.array(count_matrix)

#     seq_dict = read_fasta_by_faidx(gene_model_cds_file)
#     cds_len_list = [seq_dict[i].len() for i in gene_list]

#     tpm_matrix = get_TPM_matrix(count_matrix, cds_len_list)

#     # output
#     with open(args.output_prefix+".tpm.matrix", 'w') as f:
#         f.write("Gene\t" + "\t".join(sample_list) + "\n")

#         for g in range(len(gene_list)):
#             c_list = []
#             for s in range(len(sample_list)):
#                 c = tpm_matrix[g][s]
#                 c_list.append(str(c))
#             "\t".join(c_list)

#             g_id = gene_list[g]
#             f.write(g_id + "\t" + "\t".join(c_list) + "\n")

#     # get FPKM
#     fpkm_matrix = get_FPKM_matrix(count_matrix, cds_len_list)

#     # output
#     with open(args.output_prefix+".fpkm.matrix", 'w') as f:
#         f.write("Gene\t" + "\t".join(sample_list) + "\n")

#         for g in range(len(gene_list)):
#             c_list = []
#             for s in range(len(sample_list)):
#                 c = fpkm_matrix[g][s]
#                 c_list.append(str(c))
#             "\t".join(c_list)

#             g_id = gene_list[g]
#             f.write(g_id + "\t" + "\t".join(c_list) + "\n")
