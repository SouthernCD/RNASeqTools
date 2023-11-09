from BCBio import GFF
from pyfaidx import Fasta
from toolbiox.lib.common.fileIO import tsv_file_dict_parse
from toolbiox.lib.common.genome.seq_base import read_fasta_big
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
