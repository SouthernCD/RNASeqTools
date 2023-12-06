from toolbiox.api.common.mapping.blast import outfmt6_read_big
from toolbiox.lib.common.os import get_file_dir
from toolbiox.lib.xuyuxing.evolution.taxonomy import read_tax_record_dict_db


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
