# RNASeqTools
Some small tools to help with RNA-seq analysis

## Installation
```
pip install git+ssh://git@github.com/SouthernCD/RNASeqTools.git
```

## Usage

**DenovoCount2RefCount**

Sometimes we want to map the expression information obtained from de novo transcriptome assembly data to a reference genome of a closely related species, and DenovoCount2RefCount helps you do that!

```
RNASeqTools DenovoCount2RefCount [-h] [-f REF_GENOME_FASTA] [-g REF_GENOME_GFF] [-c REF_CDNA_FASTA] [-o WORK_DIR] [-m MAP_PROGRAM] de_novo_fasta de_novo_trans_gene_map rsem_gene_results

positional arguments:
  de_novo_fasta         a fasta file from trinity output
  de_novo_trans_gene_map
                        a gene_trans_map file from trinity output
  rsem_gene_results     a rsem gene results file

optional arguments:
  -h, --help            show this help message and exit
  -f REF_GENOME_FASTA, --ref_genome_fasta REF_GENOME_FASTA
                        refernce genome fasta
  -g REF_GENOME_GFF, --ref_genome_gff REF_GENOME_GFF
                        refernce genome gff
  -c REF_CDNA_FASTA, --ref_cDNA_fasta REF_CDNA_FASTA
                        refernce cDNA fasta, if not given, will use ref_genome_fasta and ref_genome_gff, only work for blast
  -o WORK_DIR, --work_dir WORK_DIR
                        work dir
  -m MAP_PROGRAM, --map_program MAP_PROGRAM
                        map program, blat or blast
```

## Example
### DenovoCount2RefCount
Mapping trinity results and rsem results to a reference genome using blat
```
RNASeqTools DenovoCount2RefCount -f ref_genome.fna -g ref_genome.gff -o map_out -m blat Trinity.fasta Trinity.fasta.gene_trans_map rsem.genes.results
```
Mapping trinity results and rsem results to a reference genome cDNA using blast
```
RNASeqTools DenovoCount2RefCount -c ref_cDNA.fna -o map_out -m blast Trinity.fasta Trinity.fasta.gene_trans_map rsem.genes.results
```