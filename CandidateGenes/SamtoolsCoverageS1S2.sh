# Get coverage info for each position of the genome
#
# copyright (c) - Manuel Delpero
#


# location of the files
#"/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/6/BFMI861-S1P_trimmed.aligned.sorted.dedup.recalibrated.bam"
#"/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/7/BFMI861-S1P_trimmed.aligned.sorted.dedup.recalibrated.bam"
#"/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/8/BFMI861-S2P_trimmed.aligned.sorted.dedup.recalibrated.bam"
#"/home/danny/NAS/Mouse/DNA/Sequencing/Alignment2020/9/BFMI861-S2P_trimmed.aligned.sorted.dedup.recalibrated.bam"
			 
#Samtools coverage in each position
cd /home/manuel/AIL_S1xS2/CoverageS1S2_Regions
nohup /home/manuel/Samtools/samtools-1.11/samtools depth -f /home/manuel/BAMS1LOC > coverageS1.coverage
nohup /home/manuel/Samtools/samtools-1.11/samtools depth -f /home/manuel/BAMS2LOC > coverageS2.coverage

# Since the file is huge we can create subsets to plot
awk '$1 == 1 {print $0}' deduped_MA605.coverage > chr1_MA605.coverage