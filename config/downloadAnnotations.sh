#!/usr/bin/env bash

# Modify FASTA header to include only first name (ENST)
curl http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz \
| zcat | sed s'/|.*$//' | gzip > gencode.v38.transcripts-modHeader.fa.gz

# Remove 'chr' prefix from chromosome names (e.g. chr1 -> 1)
curl http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gff3.gz \
| zcat | sed 's/^chr//' | gzip > gencode.v38.annotation-noChr.gff3.gz

# Download Hallmark gene set from MSigDB
wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.4/h.all.v7.4.symbols.gmt
