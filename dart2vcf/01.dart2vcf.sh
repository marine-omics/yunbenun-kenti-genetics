# Note: Genome is present in an external directory. This is just the genome but with bwa indexes


./dart2vcf.py -g ../../sunny_dart/vcfconvert/aten_final_0.11.fasta ../raw_data/Report_DAc23-8597_5_moreOrders_SNP_2.csv > 8597.vcf 2>8597.err

# Reheader because sample 20 is changed to WB_36

bcftools reheader --samples renamed.txt 8597.vcf | bgzip > 8597.vcf.gz
tabix 8597.vcf.gz