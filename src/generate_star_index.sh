#!/bin/bash

# Paths
GENOME_DIR="/Users/isidorabeslic/Workspace/egfrKD-AstroRNAseq/data/genomes"
FASTA_GZ="${GENOME_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
GTF_GZ="${GENOME_DIR}/Mus_musculus.GRCm39.104.gtf.gz"

# Paths for unzipped versions
FASTA="${GENOME_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
GTF="${GENOME_DIR}/Mus_musculus.GRCm39.104.gtf"

# STAR index output directory
STAR_INDEX_DIR="${GENOME_DIR}/star_index"

# Create output directory if it doesn't exist
mkdir -p "${STAR_INDEX_DIR}"

# Decompress if needed
if [ ! -f "$FASTA" ]; then
    echo "Decompressing FASTA..."
    gunzip -c "$FASTA_GZ" > "$FASTA"
fi

if [ ! -f "$GTF" ]; then
    echo "Decompressing GTF..."
    gunzip -c "$GTF_GZ" > "$GTF"
fi

# Run STAR genomeGenerate
/Users/isidorabeslic/STAR-2.7.11b/bin/MacOSX_x86_64/STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir $STAR_INDEX_DIR \
--genomeFastaFiles $FASTA \
--sjdbGTFfile $GTF \
--sjdbOverhang 74
