#!/usr/bin/env python
import pysam, os
from snakemake import input, output, params

# Open input BAM
inbam = pysam.AlignmentFile(input.bam_all, "rb")

# Prepare output BAM for filtered reads
os.makedirs(params.outdir, exist_ok=True)
outbam = pysam.AlignmentFile(output.bam_filt, "wb", template=inbam)

# Filter criteria:
#   - mapped
#   - MAPQ >= 30
#   - NM (edit distance) <= 20
#   - total matched length (sum of 'M' CIGAR ops) >= 50
for read in inbam.fetch(until_eof=True):
    if read.is_unmapped:
        continue
    if read.mapping_quality < 30:
        continue
    nm = read.get_tag("NM") if read.has_tag("NM") else 0
    if nm > 20:
        continue
    mlen = sum(length for length, typ in read.cigartuples or [] if typ == 0)
    if mlen < 50:
        continue
    outbam.write(read)

inbam.close()
outbam.close()

# Generate flagstat on the filtered BAM
os.system(f"samtools flagstat {output.bam_filt} > {output.flagstat}")

# Sort & index the filtered BAM
os.system(f"samtools sort -o {output.sorted} {output.bam_filt}")
os.system(f"samtools index {output.sorted}")
