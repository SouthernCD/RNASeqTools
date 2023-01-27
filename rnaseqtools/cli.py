import argparse
from rnaseqtools.pipelines import *


def args_init():

    # argument parse
    parser = argparse.ArgumentParser(
        prog='TaxonTools',
    )

    subparsers = parser.add_subparsers(
        title='subcommands', dest="subcommand_name")

    # argparse for ExtractTrinityGene
    parser_a = subparsers.add_parser('ExtractTrinityGene',
                                     help='Extract longest isoform from a gene in trinity output\n',
                                     description='Extract longest isoform from a gene in trinity output')

    parser_a.add_argument('Trinity_output', type=str,
                          help='a fasta file from trinity output')
    parser_a.add_argument('Trinity_gene', type=str,
                          help='a fasta file from this script')

    # argparse for CheckTrinotateStatus
    parser_a = subparsers.add_parser('CheckTrinotateStatus',
                                     help='test if a trinotate pipeline is good')

    parser_a.add_argument('transcripts_file', type=str,
                          help='path of transcripts file')
    parser_a.add_argument('-o', "--output_file",
                          type=str, help='output file path')

    # argparse for GenerateTrinityGeneTransMap
    parser_a = subparsers.add_parser('GenerateTrinityGeneTransMap',
                                     help='Generate Trinity Gene Trans Map file from trinity fasta\n'
                                     )

    parser_a.add_argument('Trinity_output', type=str,
                          help='a fasta file from trinity output')
    parser_a.add_argument('gene_map_file', type=str,
                          help='a tab file for gene_map')

    # argparse for DenovoCount2RefCount
    parser_a = subparsers.add_parser('DenovoCount2RefCount',
                                     help='Mapping the count of denovo-assembled unigene to the gene of the reference genome\n'
                                     )

    parser_a.add_argument('gene_model_cds_file', type=str,
                          help='cds file from reference genome')
    parser_a.add_argument('tran_fasta_file', type=str,
                          help='a fasta file from trinity output')
    parser_a.add_argument('gene_trans_map', type=str,
                          help='a gene_trans_map file from trinity output')
    parser_a.add_argument('-f', '--tran_count_fof', type=str, default=None,
                          help='a list of rsem genes.results file: (sampleid rsem_count_file)')
    parser_a.add_argument('-c', '--tran_count_matrix', type=str, default=None,
                          help='a count matrix')
    parser_a.add_argument('-w', '--work_dir', type=str, default='tmp',
                          help='tmp work dir')
    parser_a.add_argument('-o', '--output_prefix', type=str, default='Gene',
                          help='output prefix')
    parser_a.add_argument('-t', '--threads', type=int, default=56,
                          help='number of threads: defaults 56')

    # argparse for GetGeneLength
    parser_a = subparsers.add_parser('GetGeneLength',
                                     help='get gene model (longest mRNA) lenght from gff file\n'
                                     )

    parser_a.add_argument('gff_file', type=str,
                          help='a gff file')
    parser_a.add_argument('gene_length_file', type=str,
                          help='gene length file (output)')

    # argparse for Count2TMM
    parser_a = subparsers.add_parser('Count2TMM',
                                     help='from count matrix 2 TMM matrix, FPKM and TPM will give to\n'
                                     )

    parser_a.add_argument('count_matrix', type=str,
                          help='a count matrix')
    parser_a.add_argument('gene_length_file', type=str,
                          help='gene length file')
    parser_a.add_argument('output_prefix', type=str,
                          help='output prefix')

    # argparse for ContaminationDetector
    parser_a = subparsers.add_parser('ContaminationDetector',
                                     help='find contamination from diamond output: diamond blastp --query Trinity.model.faa --max-target-seqs 10 --db /lustre/home/xuyuxing/Database/NCBI/nr/2020/nr.taxon.dmnd --evalue 1e-5 --out Trinity.model.faa.bls --outfmt 6 qseqid sseqid staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore --threads 56\n'
                                     )

    parser_a.add_argument('bls_results_file', type=str,
                          help='diamond blast output')
    parser_a.add_argument('taxon_db_file', type=str,
                          help='taxon_db_file, from ncbi taxonomy parsed by TaxonTools')
    parser_a.add_argument('target_taxon', type=str,
                          help='can be taxon sciname or taxon id')

    args = parser.parse_args()

    return args


def main():

    args = args_init()
    args_dict = vars(args)

    if args_dict["subcommand_name"] == "ExtractTrinityGene":
        ExtractTrinityGene_main(args)

    elif args_dict["subcommand_name"] == "CheckTrinotateStatus":
        CheckTrinotateStatus_main(args)

    elif args_dict["subcommand_name"] == "GenerateTrinityGeneTransMap":
        GenerateTrinityGeneTransMap_main(args)

    elif args_dict["subcommand_name"] == "DenovoCount2RefCount":
        DenovoCount2RefCount_main(args)

    elif args_dict["subcommand_name"] == "Count2TMM":
        Count2TMM_main(args)

    elif args_dict["subcommand_name"] == "ContaminationDetector":
        ContaminationDetector_main(args)

    elif args_dict["subcommand_name"] == "GetGeneLength":
        GetGeneLength_main(args)

if __name__ == '__main__':
    main()
