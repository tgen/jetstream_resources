# contig map

This table describes the contigs found in various GRCh38/hg38 reference genome 
fastas. It was anchored by the GRCh38.p12 standard release. All contigs, 
including patch, alt, unplaced/unlocalized, from that release are found in 
the first column "refseq". The rest of the columns are the identifiers other
organizations use to identify that same contig. 


# GRC

The GRC actually publishes two separate releases of the genome:

The first, I'll call the "standard release", is updated with patches periodically. 
On their FTP, you can find GRCh38 or GRCh38.p12 for example. The .p12 indicates 
that this is patch version 12. These patch releases are also called minor 
releases. As new sequences are identified, gaps filled, etc., the GRC creates 
incremental updates that _do_not_ affect the coordinates of the other contigs. 
They are added as additional seqeunces in the fasta files. 

Additionaly, GRC has published a "pipeline release" specifically for data 
analysis pipelines. This version has different sequence names that match the 
UCSC-style naming convention. It also does not include any patches, and is not 
updated with minor releases.

Both of these are described below:


# GRC - Standard release 

Every contig is named by its RefSeq accession number, the standard builds will
also include patch contigs, notice "FIX PATCH" in some of the descriptions.

## Fasta

Source: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.27_GRCh38.p12/GCA_000001405.27_GRCh38.p12_genomic.fna.gz

Excerpt:

```
>CM000663.2 Homo sapiens chromosome 1, GRCh38 reference primary assembly
>KI270706.1 Homo sapiens chromosome 1 unlocalized genomic contig, GRCh38 reference primary assembly
>KI270707.1 Homo sapiens chromosome 1 unlocalized genomic contig, GRCh38 reference primary assembly
...
>KI270740.1 Homo sapiens chromosome Y unlocalized genomic contig, GRCh38 reference primary assembly
>KI270302.1 Homo sapiens unplaced genomic contig, GRCh38 reference primary assembly
>KI270304.1 Homo sapiens unplaced genomic contig, GRCh38 reference primary assembly
...
>KZ208924.1 Homo sapiens chromosome Y genomic contig HG1535_PATCH, GRC reference assembly FIX PATCH for GRCh38
>KN196487.1 Homo sapiens chromosome Y genomic contig HG2062_PATCH, GRC reference assembly FIX PATCH for GRCh38
>KI270762.1 Homo sapiens chromosome 1 genomic contig, GRCh38 reference assembly alternate locus group ALT_REF_LOCI_1
...
>KI270933.1 Homo sapiens chromosome 19 genomic contig, GRCh38 reference assembly alternate locus group ALT_REF_LOCI_34
>GL000209.2 Homo sapiens chromosome 19 genomic contig, GRCh38 reference assembly alternate locus group ALT_REF_LOCI_35
>J01415.2 Homo sapiens mitochondrion, complete genome
```


# GRC - Pipeline release

With the GRCh38 primary assembly release, the GRC added a folder to their 
FTP [seqs_for_alignment_pipelines.ucsc_ids](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids)
presumably because researchers would prefer to see chr1, chrX, etc. instead of 
the raw accession numbers. Some portions of the genome have been hard masked to
prevent erroneous alignments. In addition, they added several extra "decoy" sequences.
These sequences are not part of the human genome, but improve the alignment of short
read sequencing data.

## Fasta 

Source: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/eqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz

There are several fastas available in that same FTP directory. This is from one with the 
most content. The other files will omit the alt or decoy sequences for example.

Excerpt: 

```
>chr1  AC:CM000663.2  gi:568336023  LN:248956422  rl:Chromosome  M5:6aef897c3d6ff0c78aff06ac189178dd  AS:GRCh38
>chr2  AC:CM000664.2  gi:568336022  LN:242193529  rl:Chromosome  M5:f98db672eb0993dcfdabafe2a882905c  AS:GRCh38
>chr3  AC:CM000665.2  gi:568336021  LN:198295559  rl:Chromosome  M5:76635a41ea913a405ded820447d067b0  AS:GRCh38
...
>chrX  AC:CM000685.2  gi:568336001  LN:156040895  rl:Chromosome  M5:2b3a55ff7f58eb308420c8a9b11cac50  AS:GRCh38
>chrY  AC:CM000686.2  gi:568336000  LN:57227415  rl:Chromosome  M5:ce3e31103314a704255f3cd90369ecce  AS:GRCh38  hm:10001-2781479,56887903-57217415
>chrM  AC:J01415.2  gi:113200490  LN:16569  rl:Mitochondrion  M5:c68f52674c9fb33aef52dcf399755519  AS:GRCh38  tp:circular
...
>chrUn_GL000216v2  AC:GL000216.2  gi:568335254  LN:176608  rl:unplaced  M5:725009a7e3f5b78752b68afa922c090c  AS:GRCh38
>chrUn_GL000218v1  AC:GL000218.1  gi:224183305  LN:161147  rl:unplaced  M5:1d708b54644c26c7e01c2dad5426d38c  AS:GRCh38
>chr1_KI270762v1_alt  AC:KI270762.1  gi:568335926  LN:354444  rg:chr1:2448811-2791270  rl:alt-scaffold  M5:b0397179e5a9
...
>chr19_GL949753v2_alt  AC:GL949753.2  gi:568335996  LN:796479  rg:chr19:54025634-55084318  rl:alt-scaffold  M5:19162055ca3e800f14797b6cd37c1d4c  AS:GRCh38
>chr19_KI270938v1_alt  AC:KI270938.1  gi:568335998  LN:1066800  rg:chr19:54025634-55084318  rl:alt-scaffold  M5:9363b56f7b34fb35ab3400b1093f431a  AS:GRCh38
>chrEBV  AC:AJ507799.2  gi:86261677  LN:171823  rl:decoy  M5:6743bd63b3ff2b5b8985d8933c53290a  SP:Human_herpesvirus_4  tp:circular
...
>chrUn_JTFH01001996v1_decoy  AC:JTFH01001996.1  gi:725020270  LN:2009  rl:decoy  M5:a6503ea36c1b69162d3eda87c4f33922  AS:hs38d1
>chrUn_JTFH01001997v1_decoy  AC:JTFH01001997.1  gi:725020269  LN:2003  rl:decoy  M5:d3a1bed41244634725882235662d11a6  AS:hs38d1
>chrUn_JTFH01001998v1_decoy  AC:JTFH01001998.1  gi:725020268  LN:2001  rl:decoy  M5:35916d4135a2a9db7bc0749955d7161a  AS:hs38d1
```


# UCSC

Prefixes every contig name with "chr" including unplaced, unlocalized, patches, alts,
etc. The tricky part is that they also change the name of _some_ of the contigs.
Don't fall into the trap of thinking you can convert from UCSC to Ensembl IDs by just
removing the "chr". Even if you choose to ignore the alt contigs, those GL.. KI.. 
contigs will still be wrong. Here are examples:

    Ensembl     UCSC
    KI270706.1	chr1_KI270706v1_random
    KI270707.1	chr1_KI270707v1_random
    KI270708.1	chr1_KI270708v1_random
    GL000009.2	chr14_GL000009v2_random
    GL000225.1	chr14_GL000225v1_random

## Fasta

Source: ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.fa.gz

There are several fastas available in that same FTP directory. This is 
from one with the most content. The other files will omit the alt or 
unplaced contigs for example. Each contig can also be downloaded 
individually from the "chromosomes" directory found here: 
ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes

Excerpt:

```
>chr1
>chr10
>chr11
>chr11_KI270721v1_random
...
>chr19_KI270938v1_alt
>chrM
>chrUn_KI270302v1
...
>chrX
>chrY
>chrY_KI270740v1_random
...
>chr17_KZ559114v1_alt
>chr18_KZ559116v1_alt
>chr18_KZ559115v1_fix
```


# Ensembl

## Fasta

Source: ftp://ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

Ensembl uses single integers to represent the autosomes. For unplaced, 
unlocalized, and alts, they generally just use the RefSeq accession numbers. 
But, the patch contigs are not consistent. The names found in data downloads 
from Ensembl FTP are different than what you'll find while browsing their 
website or other sources like BioMart. Most contigs match, but for patches, 
the name in their downloads includes `CHR_<patchname>`, whereas the website 
lists these contigs with just the patch name. The difference might be very 
small, but it will break some simple match operations that researchers often
use.

Excerpt:

```
1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
2 dna:chromosome chromosome:GRCh38:2:1:242193529:1 REF
3 dna:chromosome chromosome:GRCh38:3:1:198295559:1 REF
...
X dna:chromosome chromosome:GRCh38:X:1:156040895:1 REF
Y dna:chromosome chromosome:GRCh38:Y:2781480:56887902:1 REF
MT dna:chromosome chromosome:GRCh38:MT:1:16569:1 REF
CHR_HG76_PATCH dna:chromosome chromosome:GRCh38:CHR_HG76_PATCH:1:144938989:1 PATCH_FIX
CHR_HSCHR15_4_CTG8 dna:chromosome chromosome:GRCh38:CHR_HSCHR15_4_CTG8:1:102071387:1 HAP
CHR_HSCHR6_MHC_SSTO_CTG1 dna:chromosome chromosome:GRCh38:CHR_HSCHR6_MHC_SSTO_CTG1:1:170946136:1 HAP
...
KI270423.1 dna:scaffold scaffold:GRCh38:KI270423.1:1:981:1 REF
KI270392.1 dna:scaffold scaffold:GRCh38:KI270392.1:1:971:1 REF
KI270394.1 dna:scaffold scaffold:GRCh38:KI270394.1:1:970:1 REF
```


## GTF

Source: ftp://ftp.ensembl.org/pub/release-95/gtf/homo_sapiens/Homo_sapiens.GRCh38.95.gtf.gz

Looking at the first column in this file, some lines redacted for more concise 
output. Note that the autosomes don't match UCSC or RefSeq names, the unplaced 
contigs match RefSeq names, and the patches don't match anyone, including 
their own database. They use their own ID system for exons, genes, transcripts 
etc. But, they also integrate HAVANA ids when applicable.

```
1
10
11
...
CHR_HG107_PATCH
CHR_HG126_PATCH
CHR_HG1311_PATCH
...
KI270744.1
KI270750.1
MT
...
```

# Gencode

Gencode is nice. They tried to add Ensembl contig names along with their 
own names. But, due to the problems with Ensembl discussed above, they don't 
exactly match the names that you'll find in files downloaded from Ensembl.

## Fasta

Source: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.p12.genome.fa.gz

Excerpt:

```
>chr1 1
>chr2 2
>chr3 3
...
>chrX X
>chrY Y
>chrM MT
>GL000008.2 GL000008.2
...
>GL000226.1 GL000226.1
>KQ759759.1 HG107_PATCH
>KN538364.1 HG126_PATCH
...
>KV766199.1 HSCHRX_3_CTG7
>KI270302.1 KI270302.1
>KI270303.1 KI270303.1
...
>KI270755.1 KI270755.1
>KI270756.1 KI270756.1
>KI270757.1 KI270757.1
```

## GTF

Source: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.chr_patch_hapl_scaff.annotation.gtf.gz

Looking at the first column in this file, some lines redacted for more concise 
output. The `gene_id`, `transcript_id`, and `exon_id` fields are Ensembl IDs. 
They also include HAVANA ids when applicable.

```
chr1
chr10
chr11
...
chrM
chrX
chrY
...
GL000009.2
GL000194.1
GL000195.1
...
KZ559114.1
KZ559115.1
KZ559116.1
```

# Other info 

## GRC vs UCSC

It looks like UCSC was publishing its own genome builds for hg1-hg8, but 
later adopted the NCBI builds this table 
[here](https://genome.ucsc.edu/FAQ/FAQreleases.html). In December 2001, with 
the release of NCBI Build 28, UCSC lists the NCBI build IDs associated with the
their own IDs hg10-19. Finally, upon the release of build 38, they harmonized 
their build numbers with the GRC and called the latest release hg38.

## Sequence descriptions

Often times the sequence descriptions will hold useful information. There
is actually a standard published by the NCBI but nobody really uses it as far 
as I can tell:

From wikipedia FASTA_format:,

    NCBI identifiers
    The NCBI defined a standard for the unique identifier used for the sequence 
    (SeqID) in the header line. This allows a sequence that was obtained from a 
    database to be labelled with a reference to its database record. The database 
    identifier format is understood by the NCBI tools like makeblastdb and 
    table2asn. The following list describes the NCBI FASTA defined format for 
    sequence identifiers.

Type	Format(s)	Example(s)
local (i.e. no database reference)	lcl|integer lcl|string lcl|123 lcl|hmm271
GenInfo backbone seqid	bbs|integer	bbs|123
GenInfo backbone moltype	bbm|integer	bbm|123
GenInfo import ID	gim|integer	gim|123
GenBank	gb|accession|locus	gb|M73307|AGMA13GT
EMBL	emb|accession|locus	emb|CAM43271.1|
PIR	pir|accession|name	pir||G36364
SWISS-PROT	sp|accession|name	sp|P01013|OVAX_CHICK
patent	pat|country|patent|sequence-number	pat|US|RE33188|1
pre-grant patent	pgp|country|application-number|sequence-number	pgp|EP|0238993|7
RefSeq	ref|accession|name	ref|NM_010450.1|
general database reference (a reference to a database that's not in this list)	gnl|database|integer gnl|database|string gnl|taxon|9606 gnl|PID|e1632
GenInfo integrated database	gi|integer	gi|21434723
DDBJ	dbj|accession|locus	dbj|BAC85684.1|
PRF	prf|accession|name	prf||0806162C
PDB	pdb|entry|chain	pdb|1I4L|D
third-party GenBank	tpg|accession|name	tpg|BK003456|
third-party EMBL	tpe|accession|name	tpe|BN000123|
third-party DDBJ	tpd|accession|name	tpd|FAA00017|
TrEMBL	tr|accession|name	tr|Q90RT2|Q90RT2_9HIV1

