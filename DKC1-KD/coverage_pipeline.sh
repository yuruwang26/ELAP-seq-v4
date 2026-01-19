#!/bin/bash 

# Input arguments
BAM=$1            # e.g., IP.bam
FORWARD_BED=$2    # e.g., pos.bed
REVERSE_BED=$3    # e.g., neg.bed
OUTPUT=$4         # e.g., combined_sorted_coverage.txt

# Temporary files
FORWARD_BAM="forward_reads.bam"
REVERSE_BAM="reverse_reads.bam"
FORWARD_COV="forward_coverage.txt"
REVERSE_COV="reverse_coverage.txt"

# Step 1: Split BAM into forward and reverse strand reads using XS:A:+ (for forward) and XS:A:- (for reverse)
samtools view -h "$BAM" | awk '{if($0 ~ /^@/ || $0 ~ /XS:A:+/) print $0}' | samtools view -bS > "$FORWARD_BAM"  # Forward strand
samtools view -h "$BAM" | awk '{if($0 ~ /^@/ || $0 ~ /XS:A:-/) print $0}' | samtools view -bS > "$REVERSE_BAM"  # Reverse strand

# Step 2: Get coverage for each strand using the provided strand information in the BED files
# Use -a to include all positions, even those with zero coverage
samtools depth -a -b "$FORWARD_BED" "$FORWARD_BAM" > "$FORWARD_COV"
samtools depth -a -b "$REVERSE_BED" "$REVERSE_BAM" > "$REVERSE_COV"

# Step 3: Format coverage output to be compatible with BED format
# Forward strand: Convert samtools depth output to BED format (chrom, start, end, coverage)
awk '{print $1"\t"$2-1"\t"$2"\t"$3}' "$FORWARD_COV" > formatted_forward_coverage.bed

# Reverse strand: Convert samtools depth output to BED format (chrom, start, end, coverage)
awk '{print $1"\t"$2-1"\t"$2"\t"$3}' "$REVERSE_COV" > formatted_reverse_coverage.bed

# Step 4: Add missing positions with zero coverage using bedtools intersect
# Forward strand: Intersect with forward BED file to add missing positions with coverage 0
bedtools intersect -wa -a formatted_forward_coverage.bed -b "$FORWARD_BED" > overlap_entries.txt

bedtools intersect -a "$FORWARD_BED" -b formatted_forward_coverage.bed -v | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "0"}' > no_overlap_entries.txt

cat overlap_entries.txt no_overlap_entries.txt | sort -k1,1 -k2,2n > full_coverage_forward.txt

# Reverse strand: Intersect with reverse BED file to add missing positions with coverage 0
bedtools intersect -wa -a formatted_reverse_coverage.bed -b "$REVERSE_BED" > overlap_entries.txt

bedtools intersect -a "$REVERSE_BED" -b formatted_reverse_coverage.bed -v | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, "0"}' > no_overlap_entries.txt

cat overlap_entries.txt no_overlap_entries.txt | sort -k1,1 -k2,2n > full_coverage_reverse.txt

# Step 5: Combine the two tables and sort by chromosomal location (chrom, position)
cat full_coverage_forward.txt full_coverage_reverse.txt | sort -k1,1 -k2,2n > "$OUTPUT"

# Step 6: Cleanup
rm "$FORWARD_BAM" "$REVERSE_BAM" "$FORWARD_COV" "$REVERSE_COV" formatted_forward_coverage.bed formatted_reverse_coverage.bed full_coverage_forward.txt full_coverage_reverse.txt overlap_entries.txt no_overlap_entries.txt

# Final message
echo "Coverage pipeline complete!"
echo "Combined and sorted coverage table saved to: $OUTPUT"

