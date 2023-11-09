import argparse


def args_init():

    # argument parse
    parser = argparse.ArgumentParser(
        prog='RNASeqTools',
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
                                     help='Test if a trinotate pipeline is good')

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

    # argparse for GenerateSpadesGeneTransMap
    parser_a = subparsers.add_parser('GenerateSpadesGeneTransMap',
                                     help='Generate Spades Gene Trans Map file from Spades fasta\n'
                                     )

    parser_a.add_argument('Spades_transcripts_fasta', type=str,
                          help='transcripts.fasta from rnaspades')
    parser_a.add_argument('gene_map_file', type=str,
                          help='a tab file for gene_map')

    # argparse for DenovoCount2RefCount
    parser_a = subparsers.add_parser('DenovoCount2RefCount',
                                     help='Mapping the count of denovo-assembled unigene to the gene of the reference genome\n'
                                     )

    parser_a.add_argument('de_novo_fasta', type=str,
                          help='a fasta file from trinity output')
    parser_a.add_argument('de_novo_trans_gene_map', type=str,
                          help='a gene_trans_map file from trinity output')
    parser_a.add_argument('rsem_gene_results', type=str,
                          help='a rsem gene results file')
    parser_a.add_argument('-f', '--ref_genome_fasta', type=str, default=None,
                          help='refernce genome fasta')
    parser_a.add_argument('-g', '--ref_genome_gff', type=str, default=None,
                          help='refernce genome gff')
    parser_a.add_argument('-c', '--ref_cDNA_fasta', type=str, default=None,
                          help='refernce cDNA fasta, if not given, will use ref_genome_fasta and ref_genome_gff, only work for blast')
    parser_a.add_argument('-o', '--work_dir', type=str, default='DenovoCount2RefCount_out',
                          help='work dir')
    parser_a.add_argument('-m', '--map_program', type=str, default='blat',
                          help='map program, blat or blast')

    # argparse for GetGeneLength
    parser_a = subparsers.add_parser('GetGeneLength',
                                     help='get gene model (longest mRNA) length from gff file\n'
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
        from rnaseqtools.src.trinity_tools import ExtractTrinityGene_main
        ExtractTrinityGene_main(args)

    elif args_dict["subcommand_name"] == "CheckTrinotateStatus":
        from rnaseqtools.src.trinity_tools import CheckTrinotateStatus_main
        CheckTrinotateStatus_main(args)

    elif args_dict["subcommand_name"] == "GenerateTrinityGeneTransMap":
        from rnaseqtools.src.trinity_tools import GenerateTrinityGeneTransMap_main
        GenerateTrinityGeneTransMap_main(args)

    elif args_dict["subcommand_name"] == "DenovoCount2RefCount":
        from rnaseqtools.src.DenovoCount2RefCount import DenovoCount2RefCount_main
        DenovoCount2RefCount_main(args)

    elif args_dict["subcommand_name"] == "Count2TMM":
        from rnaseqtools.src.small_tools import Count2TMM_main
        Count2TMM_main(args)

    elif args_dict["subcommand_name"] == "ContaminationDetector":
        from rnaseqtools.src.ContaminationDetector import ContaminationDetector_main
        ContaminationDetector_main(args)

    elif args_dict["subcommand_name"] == "GetGeneLength":
        from rnaseqtools.src.small_tools import GetGeneLength_main
        GetGeneLength_main(args)


if __name__ == '__main__':
    main()
