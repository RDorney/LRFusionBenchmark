-------------
COSMIC Fusion
-------------
 All gene fusion mutation data from the current release in a tab separated table. This file includes all the tested samples, with and without fusion detected.  [ Cosmic_Fusion_v102_GRCh38.tsv.gz ]

  File Description

[column number:label] Heading                             Description
--------------------------------------------------------------------------------------------------------
[1:A]                 COSMIC_SAMPLE_ID                    A unique COSMIC sample identifier (COSS) is used to identify a sample. This identifier can be used to retrieve additional Sample information from the Cosmic_Sample file.
[2:B]                 SAMPLE_NAME                         The sample name can be derived from a number of sources. In many cases it originates from the cell line name. Other sources include names assigned by the annotators, or an incremented number assigned during an anonymization process.
[3:C]                 COSMIC_PHENOTYPE_ID                 A unique COSMIC identifier (COSO) for the classification. This identifier can be used to retrieve tissue and histology information from the classification file.
[4:D]                 COSMIC_FUSION_ID                    A fusion mutation identifier (COSF). This identifier can be null for samples tested but where no fusion was detected.
[5:E]                 FUSION_SYNTAX                       Syntax describing the portions of mRNA present (in HGVS 'r.' format) from each gene (allows representation of UTR sequences).
[6:F]                 FIVE_PRIME_CHROMOSOME               Chromosome of 5' gene.
[7:G]                 FIVE_PRIME_STRAND                   Positive or negative of the 5' gene (+/-).
[8:H]                 FIVE_PRIME_TRANSCRIPT_ID            The Ensembl Transcript identifier (ENST) of the 5' gene. This identifier can be used to retrieve additional Transcript information from the Cosmic_Transcripts file.
[9:I]                 FIVE_PRIME_GENE_SYMBOL              Gene symbol for the 5' gene fusion partner for which the data has been curated in COSMIC. In most cases this is the accepted HGNC identifier.
[10:J]                FIVE_PRIME_LAST_OBSERVE_EXON        Last observed exon number of the 5' gene fusion partner.
[11:K]                FIVE_PRIME_GENOME_START_FROM        The genomic coordinate of the start (+ strand)/breakpoint (- strand) of the 5' fusion gene as described in the fusion syntax.
[12:L]                FIVE_PRIME_GENOME_START_TO          The range of genomic coordinates of the start (+ strand)/breakpoint (- strand) of the 5' fusion gene if it is an unknown base position.
[13:M]                FIVE_PRIME_GENOME_STOP_FROM         The genomic coordinate of the breakpoint (+ strand)/start (- strand) of the 5' fusion gene as described in the Translocation Name.
[14:N]                FIVE_PRIME_GENOME_STOP_TO           The range of genomic coordinates of the breakpoint (+ strand)/start (- strand) of the 5' fusion gene if it is an unknown base position.
[15:O]                THREE_PRIME_CHROMOSOME              Chromosome of 3' gene.
[16:P]                THREE_PRIME_STRAND                  Positive or negative of the 3' gene (+/-).
[17:Q]                THREE_PRIME_TRANSCRIPT_ID           The Ensembl Transcript identifier (ENST) of the 3' gene. This identifier can be used to retrieve additional Transcript information from the Cosmic_Transcripts file.
[18:R]                THREE_PRIME_GENE_SYMBOL             Gene symbol for the 3' gene fusion partner  for which the data has been curated in COSMIC. In most cases this is the accepted HGNC identifier.
[19:S]                THREE_PRIME_FIRST_OBSERVE_EXON      First observed exon number of the 3' gene fusion partner.
[20:T]                THREE_PRIME_GENOME_START_FROM       The genomic coordinate of the breakpoint (+ strand)/stop (- strand) of the 3' fusion gene as described in the Translocation Name.
[21:U]                THREE_PRIME_GENOME_START_TO         The range of genomic coordinates of the breakpoint (+ strand)/stop (- strand) of the 3' fusion gene if it is an unknown base position.
[22:V]                THREE_PRIME_GENOME_STOP_FROM        The genomic coordinate of the stop (+ strand)/breakpoint (- strand) of the 3' fusion gene as described in the Translocation Name.
[23:W]                THREE_PRIME_GENOME_STOP_TO          The range of genomic coordinates of the stop (+ strand)/breakpoint (- strand) of the 3' fusion gene if it is an unknown base position.
[24:X]                FUSION_TYPE                         Type of mutation.
[25:Y]                PUBMED_PMID                         The PUBMED ID for the paper that the sample was noted in.
