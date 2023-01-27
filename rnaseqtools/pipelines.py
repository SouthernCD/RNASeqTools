from BCBio import GFF
from pyfaidx import Fasta
from toolbiox.api.common.genome.blast import outfmt6_read_big, hit_CIP, hit_CALP
from toolbiox.lib.common.fileIO import tsv_file_dict_parse, tsv_file_dict_parse_big, read_matrix_file
from toolbiox.lib.common.genome.seq_base import read_fasta_big, read_fasta_by_faidx
from toolbiox.lib.common.genome.genome_feature2 import read_gff_file, Gene
from toolbiox.lib.common.os import cmd_run, mkdir, get_file_dir
from toolbiox.api.xuyuxing.transcriptome.FPKM import get_TPM_matrix, get_FPKM_matrix
from toolbiox.lib.xuyuxing.evolution.taxonomy import read_tax_record_dict_db
import numpy as np
import os
import re


def extract_longest_isofrom_from_gene_in_trinity(trinity_output):
    longest_dict = {}
    for record in read_fasta_big(trinity_output):
        cluster_id, gene_id, isoform_id = re.match(r'^(\S+c\d+)_(g\d+)_(i\d+)$',
                                                   record.seqname_short()).groups()
        gene_name = cluster_id + "_" + gene_id
        if not gene_name in longest_dict:
            longest_dict[gene_name] = record
        elif len(longest_dict[gene_name].seq) < len(record.seq):
            longest_dict[gene_name] = record
    return longest_dict


def ExtractTrinityGene_main(args):
    Trinity_output = args.Trinity_output
    Trinity_gene = args.Trinity_gene
    output = extract_longest_isofrom_from_gene_in_trinity(Trinity_output)
    with open(Trinity_gene, "w") as f:
        for i in output:
            i = output[i]
            i.wrap(60)
            f.write(">%s\n%s" % (i.seqname, i.seq))


def GenerateTrinityGeneTransMap_main(args):
    Trinity_output = args.Trinity_output
    Trinity_gene = args.gene_map_file
    with open(Trinity_output, 'r') as f:
        with open(Trinity_gene, 'w') as w:
            for each_line in f:
                each_line = each_line.strip()
                record_head = re.match(r"^>", each_line)
                if record_head:
                    seqname = re.sub(r'^>', '', each_line)
                    name_short = re.search('^(\S+)', seqname).group(1)
                    gene_name = re.sub(r'\_i\d+$', '', name_short)
                    w.write("%s\t%s\n" % (gene_name, name_short))


def CheckTrinotateStatus_main(args):
    work_dir = os.path.dirname(args.transcripts_file)

    check_list_dir = {
        "TransDecoder": {
            "bed_file": "uncheck",
            "cds_file": "uncheck",
            "gff3_file": "uncheck",
            "pep_file": "uncheck"
        },
        "blast": {
            "blastx": "uncheck",
            "blastp": "uncheck"
        },
        "pfam": {
            "pfam_log": "uncheck",
            "pfam_tab": "uncheck"
        },
        "signalp": {
            "signalp": "uncheck"
        },

        "tmhmm": {
            "tmhmm": "uncheck"
        },
        "rnammer": {
            "rnammer_out": "uncheck",
            "rnammer_gff": "uncheck",
        },
        "trinotate": {
            "trinotate_tab": "uncheck",
            "trinotate_go": "uncheck",
        }
    }

    # TransDecoder
    transdecoder_files = {
        "bed_file": args.transcripts_file + ".transdecoder.bed",
        "cds_file": args.transcripts_file + ".transdecoder.cds",
        "gff3_file": args.transcripts_file + ".transdecoder.gff3",
        "pep_file": args.transcripts_file + ".transdecoder.pep"
    }

    for i in transdecoder_files:
        file_name = transdecoder_files[i]
        if os.path.exists(file_name):
            if i == 'bed_file':
                tran_id_list = []
                with open(file_name, 'r') as f:
                    for each_line in f:
                        if re.match(r'^track name=', each_line):
                            continue
                        each_line = re.sub('\n', '', each_line)
                        tran_id = each_line.split('\t')[0]
                        tran_id_list.append(tran_id)
                check_list_dir["TransDecoder"][i] = len(
                    list(set(tran_id_list)))
            elif i == "gff3_file":
                num = 0
                with open(file_name, 'r') as f:
                    for rec in GFF.parse(f):
                        num = num + 1
                check_list_dir["TransDecoder"][i] = num
            else:
                check_list_dir["TransDecoder"][i] = len(
                    Fasta(file_name).keys())
        else:
            check_list_dir["TransDecoder"][i] = 'failed'

    # Blast/diamond
    blast_files = {
        "blastx": work_dir + "/swissprot.blastx.outfmt6",
        "blastp": work_dir + "/swissprot.blastp.outfmt6"
    }

    for i in blast_files:
        file_name = blast_files[i]
        if os.path.exists(file_name):
            tran_id_list = []
            with open(file_name, 'r') as f:
                for each_line in f:
                    each_line = re.sub('\n', '', each_line)
                    tran_id = each_line.split('\t')[0]
                    tran_id_list.append(tran_id)
            check_list_dir["blast"][i] = len(list(set(tran_id_list)))
        else:
            check_list_dir["blast"][i] = 'failed'

    # pfam
    pfam_files = {
        "pfam_log": work_dir + "/pfam.log",
        "pfam_tab": work_dir + "/TrinotatePFAM.out"
    }

    for i in pfam_files:
        file_name = pfam_files[i]
        if os.path.exists(file_name):
            if i == "pfam_log":
                tran_id_list = []
                with open(file_name, 'r') as f:
                    for each_line in f:
                        each_line = re.sub('\n', '', each_line)
                        match_list = re.findall(
                            r'^Query:\s+(\S+)\s+\S+$', each_line)
                        if len(match_list) != 0:
                            tran_id = match_list[0]
                            tran_id_list.append(tran_id)
                check_list_dir["pfam"][i] = len(list(set(tran_id_list)))
            elif i == "pfam_tab":
                tran_id_list = []
                with open(file_name, 'r') as f:
                    for each_line in f:
                        if re.match(r'^#', each_line):
                            continue
                        each_line = re.sub('\n', '', each_line)
                        tran_id = each_line.split()[3]
                        tran_id_list.append(tran_id)
                check_list_dir["pfam"][i] = len(list(set(tran_id_list)))
        else:
            check_list_dir["pfam"][i] = 'failed'

    # SIGNALP
    signalp_files = {
        "signalp": work_dir + "/signalp.out_summary.signalp5"
    }

    for i in signalp_files:
        file_name = signalp_files[i]
        if os.path.exists(file_name):
            tran_id_list = []
            with open(file_name, 'r') as f:
                for each_line in f:
                    if re.match(r'^#', each_line):
                        continue
                    each_line = re.sub('\n', '', each_line)
                    tran_id = each_line.split()[0]
                    tran_id_list.append(tran_id)
            check_list_dir["signalp"][i] = len(list(set(tran_id_list)))
        else:
            check_list_dir["signalp"][i] = 'failed'

    # TMHMM
    tmhmm_files = {
        "tmhmm": work_dir + "/tmhmm.out"
    }

    for i in tmhmm_files:
        file_name = tmhmm_files[i]
        if os.path.exists(file_name):
            tran_id_list = []
            with open(file_name, 'r') as f:
                for each_line in f:
                    if re.match(r'^#', each_line):
                        continue
                    each_line = re.sub('\n', '', each_line)
                    tran_id = each_line.split()[0]
                    tran_id_list.append(tran_id)
            check_list_dir["tmhmm"][i] = len(list(set(tran_id_list)))
        else:
            check_list_dir["tmhmm"][i] = 'failed'

    # RNAMMER
    rnammer_files = {
        "rnammer_out": work_dir + "/tmp.superscaff.rnammer.gff",
        "rnammer_gff": args.transcripts_file + ".rnammer.gff",
    }

    for i in rnammer_files:
        file_name = rnammer_files[i]
        if os.path.exists(file_name):
            tran_id_list = []
            num = 0
            with open(file_name, 'r') as f:
                for rec in GFF.parse(f):
                    for gene in rec.features:
                        num = num + 1
            check_list_dir["rnammer"][i] = num
        else:
            check_list_dir["rnammer"][i] = 'failed'

    # TRINOTATE report
    trinotate_files = {
        "trinotate_tab": work_dir + "/Trinotate.xls",
        "trinotate_go": work_dir + "/Trinotate.xls.gene_ontology"
    }

    col_name = ['sprot_Top_BLASTX_hit',
                'RNAMMER',
                'prot_id',
                'prot_coords',
                'sprot_Top_BLASTP_hit',
                'Pfam',
                'SignalP',
                'TmHMM',
                'eggnog',
                'Kegg',
                'gene_ontology_blast',
                'gene_ontology_pfam',
                'transcript',
                'peptide']

    for i in trinotate_files:
        file_name = trinotate_files[i]
        if os.path.exists(file_name):
            if i == "trinotate_tab":

                col_hit_num_dir = {}

                file_dict = tsv_file_dict_parse(file_name)

                for j in col_name:
                    gene_num = len(
                        list(set([file_dict[k]['#gene_id'] for k in file_dict if file_dict[k][j] != '.'])))
                    transcript_num = len(
                        list(set([file_dict[k]['transcript_id'] for k in file_dict if file_dict[k][j] != '.'])))
                    col_hit_num_dir[j] = (gene_num, transcript_num)

                check_list_dir["trinotate"][i] = col_hit_num_dir
            elif i == "trinotate_go":
                tran_id_list = []
                with open(file_name, 'r') as f:
                    for each_line in f:
                        if re.match(r'^#', each_line):
                            continue
                        each_line = re.sub('\n', '', each_line)
                        tran_id = each_line.split()[0]
                        tran_id_list.append(tran_id)
                check_list_dir["trinotate"][i] = len(tran_id_list)
        else:
            check_list_dir["trinotate"][i] = 'failed'

    if args.output_file is not None:
        import json

        with open(args.output_file, 'w') as f:
            json.dump(check_list_dir, f)

    # make output
    if check_list_dir["trinotate"]["trinotate_tab"] == 'failed':
        print("FAILED: complete failed")
    else:
        num = 0
        for i in col_name:
            if check_list_dir["trinotate"]["trinotate_tab"][i] == (0, 0):
                num = num + 1
        if num > 4:
            print("FAILED: too many field failed")
        else:
            print("OK: with %d bad" % num)


def hit_sum_hsp_aln(hit):
    aln_list = []
    for hsp in hit.hsp:
        aln_list.append(hsp.Hsp_align_len)
    return sum(aln_list)


def DenovoCount2RefCount_main(args):

    args.work_dir = os.path.abspath(args.work_dir)

    # build env
    mkdir(args.work_dir, False)

    gene_model_cds_file = args.work_dir + "/gene_model_cds.fna"
    tran_seq_file = args.work_dir + "/trans.fna"
    gene_trans_map_file = args.work_dir + "/gene_trans_map"
    os.symlink(os.path.abspath(args.gene_model_cds_file), gene_model_cds_file)
    os.symlink(os.path.abspath(args.tran_fasta_file), tran_seq_file)
    os.symlink(os.path.abspath(args.gene_trans_map), gene_trans_map_file)

    if args.tran_count_fof:
        tran_count_fof = args.work_dir + "/trans.counts.fof"
        os.symlink(os.path.abspath(args.tran_count_fof), tran_count_fof)
    else:
        tran_count_fof = None

    if args.tran_count_matrix:
        tran_count_matrix = args.work_dir + "/tran_count_matrix"
        os.symlink(os.path.abspath(args.tran_count_matrix), tran_count_matrix)
    else:
        tran_count_matrix = None

    # blast
    cmd_string = "makeblastdb -in %s -dbtype nucl" % gene_model_cds_file
    cmd_run(cmd_string, cwd=args.work_dir)

    bls_out_file = args.work_dir + "/trans_vs_cds.bls"
    cmd_string = "blastn -query %s -db %s -out %s -outfmt 6 -evalue 1e-5 -max_target_seqs 10 -num_threads %d" % (
        tran_seq_file, gene_model_cds_file, bls_out_file, args.threads)
    cmd_run(cmd_string, cwd=args.work_dir)

    # id map
    seq_dict = read_fasta_by_faidx(tran_seq_file)
    len_dict = {i: seq_dict[i].len() for i in seq_dict}

    all_tg2tt_dict = {}
    for tmp_id, tmp_dict in tsv_file_dict_parse_big(gene_trans_map_file, fieldnames=['tg', 'tt'], key_col='tt'):
        tg = tmp_dict['tg']
        all_tg2tt_dict.setdefault(tg, []).append(tmp_dict['tt'])

    tg2tt_dict = {}
    for tg in all_tg2tt_dict:
        tg2tt_dict[tg] = sorted(
            all_tg2tt_dict[tg], key=lambda x: len_dict[x], reverse=True)[0]

    tt2tg_dict = {tg2tt_dict[tg]: tg for tg in tg2tt_dict}

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

    # read tran count
    if tran_count_fof:
        rsem_file_dict = {}
        sample_list = []
        with open(tran_count_fof, 'r') as f:
            for l in f:
                l.strip()
                s, r_file = l.split()
                rsem_file_dict[s] = r_file
                sample_list.append(s)

        tran_count_dict = {}
        for s in rsem_file_dict:
            r_file = rsem_file_dict[s]
            for tmp_id, tmp_dict in tsv_file_dict_parse_big(r_file, key_col='gene_id'):
                tran_count_dict.setdefault(
                    tmp_id, {sd: 0 for sd in sample_list})
                tran_count_dict[tmp_id][s] = round(
                    float(tmp_dict['expected_count']))

    elif tran_count_matrix:
        ttran_count_dict = {}
        for tmp_id, tmp_dict in tsv_file_dict_parse_big(tran_count_matrix, key_col=''):
            sample_list = sorted([i for i in tmp_dict if i != ''])
            ttran_count_dict[tmp_id] = {
                s: round(float(tmp_dict[s])) for s in sample_list}

        tran_count_dict = {}
        for tg in all_tg2tt_dict:
            tran_count_dict[tg] = {s: 0 for s in sample_list}
            for tt in all_tg2tt_dict[tg]:
                if tt in ttran_count_dict:
                    for s in ttran_count_dict[tt]:
                        tran_count_dict[tg][s] += ttran_count_dict[tt][s]

    # map count to gene
    gene_count_dict = {}
    for g in g2t_dict:
        gene_count_dict[g] = {i: 0 for i in sample_list}
        for t in g2t_dict[g]:
            sample_count_dict = tran_count_dict[t]
            for s in sample_count_dict:
                gene_count_dict[g][s] += tran_count_dict[t][s]

    # output
    with open(args.output_prefix+".counts.matrix", 'w') as f:
        f.write("Gene\t" + "\t".join(sample_list) + "\n")

        for g in gene_count_dict:
            c_list = []
            for s in sample_list:
                c = gene_count_dict[g][s]
                c_list.append(str(c))
            "\t".join(c_list)

            f.write(g + "\t" + "\t".join(c_list) + "\n")

    # get TPM
    gene_list = sorted(list(gene_count_dict.keys()))
    count_matrix = []
    for g in gene_list:
        g_r = []
        for s in sample_list:
            g_r.append(gene_count_dict[g][s])
        count_matrix.append(g_r)
    count_matrix = np.array(count_matrix)

    seq_dict = read_fasta_by_faidx(gene_model_cds_file)
    cds_len_list = [seq_dict[i].len() for i in gene_list]

    tpm_matrix = get_TPM_matrix(count_matrix, cds_len_list)

    # output
    with open(args.output_prefix+".tpm.matrix", 'w') as f:
        f.write("Gene\t" + "\t".join(sample_list) + "\n")

        for g in range(len(gene_list)):
            c_list = []
            for s in range(len(sample_list)):
                c = tpm_matrix[g][s]
                c_list.append(str(c))
            "\t".join(c_list)

            g_id = gene_list[g]
            f.write(g_id + "\t" + "\t".join(c_list) + "\n")

    # get FPKM
    fpkm_matrix = get_FPKM_matrix(count_matrix, cds_len_list)

    # output
    with open(args.output_prefix+".fpkm.matrix", 'w') as f:
        f.write("Gene\t" + "\t".join(sample_list) + "\n")

        for g in range(len(gene_list)):
            c_list = []
            for s in range(len(sample_list)):
                c = fpkm_matrix[g][s]
                c_list.append(str(c))
            "\t".join(c_list)

            g_id = gene_list[g]
            f.write(g_id + "\t" + "\t".join(c_list) + "\n")


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


def ContaminationDetector_main(args):
    """
    running diamond
    diamond blastp --query Trinity.model.faa --max-target-seqs 10 --db /lustre/home/xuyuxing/Database/NCBI/nr/2020/nr.taxon.dmnd --evalue 1e-5 --out Trinity.model.faa.bls --outfmt 6 qseqid sseqid staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads 56

    give me target species taxon

    bls_results_file = '/lustre/home/xuyuxing/Work/Orobanchaceae/Trans/clean_data/Tve/Trinity.model.faa.bls'
    taxon_db_file = '/lustre/home/xuyuxing/Database/NCBI/taxonomy/tax_xyx.db'

    target_taxon = 'Lamiales'
    """

    bls_results_file, taxon_db_file, target_taxon = args.bls_results_file, args.taxon_db_file, args.target_taxon

    taxon_dict = read_tax_record_dict_db(taxon_db_file)

    try:
        target_taxon = str(int(target_taxon))
    except:
        t = [i for i in taxon_dict if taxon_dict[i].sci_name == target_taxon]
        if len(t) == 0:
            target_taxon = None
        else:
            target_taxon = t[0]

    contaminate_flag_dict = {}

    for query in outfmt6_read_big(bls_results_file, fieldname=["qseqid", "sseqid", "staxids", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"], gzip_flag=False):
        high_sim_taxon = []
        for hit in query.hit:
            hit_taxon_list = hit.Hit_taxon_id
            for hsp in hit.hsp:
                if hsp.Hsp_identical_ratio > 0.98 and hsp.Hsp_align_len > 50:
                    high_sim_taxon.extend(hit_taxon_list)
        high_sim_taxon = list(set(high_sim_taxon))
        if "" in high_sim_taxon:
            high_sim_taxon.remove("")

        contaminate_flag = False
        if len(high_sim_taxon) > 0:
            contaminate_flag = True
            for t_id in high_sim_taxon:
                if t_id in taxon_dict:
                    t = taxon_dict[t_id]
                    if target_taxon in set(i[0] for i in t.get_lineage):
                        contaminate_flag = False

        contaminate_flag_dict[query.qDef] = contaminate_flag

    with open(get_file_dir(bls_results_file)+"/contamination.id", 'w') as f:
        for i in contaminate_flag_dict:
            if contaminate_flag_dict[i]:
                f.write(i+"\n")

    with open(get_file_dir(bls_results_file)+"/passed.id", 'w') as f:
        for i in contaminate_flag_dict:
            if contaminate_flag_dict[i] is False:
                f.write(i+"\n")

def GetGeneLength_main(args):
    gff_dict = read_gff_file(args.gff_file)

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

    with open(args.gene_length_file, 'w') as f:
        for i in gene_length_dict:
            f.write("%s\t%d\n" % (i, gene_length_dict[i]))
