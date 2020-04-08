# S1 subset ,merge and index
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7256-1_BFMI861-S1_S2_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 17:7725897-26054796 > chr17.bam
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7256-1_BFMI861-S1_S2_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 15:66188210-78279750 > chr15.bam
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7256-1_BFMI861-S1_S2_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 3:93516258-101097858 > chr3.bam

samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7256-2_BFMI861-S1_S4_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 17:7725897-26054796 > chr17S1_S4.bam
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7256-2_BFMI861-S1_S4_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 15:66188210-78279750 > chr15S1_S4.bam
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7256-2_BFMI861-S1_S4_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 3:93516258-101097858 > chr3S1_S4.bam

samtools merge chr17_S1.bam chr17S1_S2.bam chr17S1_S4.bam 
samtools merge chr15_S1.bam chr15S1_S2.bam chr15S1_S4.bam 
samtools merge chr3_S1.bam chr3S1_S2.bam chr3S1_S4.bam 

samtools index chr17_S1.bam
samtools index chr15_S1.bam
samtools index chr3_S1.bam

# S2 subset and merge
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7257-1_BFMI861-S2_S5_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 17:7725897-26054796 > chr17S2_S5.bam
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7257-1_BFMI861-S2_S5_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 15:66188210-78279750 > chr15S2_S5.bam
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7257-1_BFMI861-S2_S5_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 3:93516258-101097858 > chr3S2_S5.bam

samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7257-2_BFMI861-S2_S6_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 17:7725897-26054796 > chr17S2_S6.bam
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7257-2_BFMI861-S2_S6_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 15:66188210-78279750 > chr15S2_S6.bam
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7257-2_BFMI861-S2_S6_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 3:93516258-101097858 > chr3S2_S6.bam

samtools merge chr17_S2.bam chr17S2_S5.bam chr17S2_S6.bam 
samtools merge chr15_S2.bam chr15S2_S5.bam chr15S2_S6.bam 
samtools merge chr3_S2.bam chr3S2_S5.bam chr3S2_S6.bam 

samtools index chr17_S2.bam
samtools index chr15_S2.bam
samtools index chr3_S2.bam

samtools index chr17_S1

# S12 subset and merge
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7258-1_BFMI860-S12_S1_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 17:7725897-26054796 > chr17S12_S1.bam
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7258-1_BFMI860-S12_S1_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 15:66188210-78279750 > chr15S12_S1.bam 
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7258-1_BFMI860-S12_S1_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 3:93516258-101097858 > chr3S12_S1.bam

samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7258-2_BFMI860-S12_S2_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 17:7725897-26054796 > chr17S12_S2.bam 
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7258-2_BFMI860-S12_S2_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 15:66188210-78279750 > chr15S12_S2.bam
samtools view -h /home/danny/NAS/Mouse/DNA/Sequencing/MouseGenomeProject/BAM/ext_L7258-2_BFMI860-S12_S2_RP_trimmed.aligned.sorted.realigned.recalibrated.bam 3:93516258-101097858 > chr3S12_S2.bam

samtools merge chr17_S12.bam chr17S12_S1.bam chr17S12_S2.bam
samtools merge chr15_S12.bam chr15S12_S1.bam chr15S12_S2.bam  
samtools merge chr3_S12.bam chr3S12_S1.bam chr3S12_S2.bam

samtools index chr17_S12.bam && samtools index chr15_S12.bam && samtools index chr3_S12.bam

# merge the three bam files for each line subset and merge
samtools merge chr17153_S12.bam chr17_S12.bam chr15_S12.bam chr3_S12.bam
samtools merge chr17153_S1.bam chr17_S1.bam  chr15_S1.bam chr3_S1.bam 
samtools merge chr17153_S2.bam chr17_S2.bam chr15_S2.bam chr3_S2.bam

samtools index chr17153_S12.bam && samtools index chr17153_S1.bam && samtools index chr17153_S2.bam