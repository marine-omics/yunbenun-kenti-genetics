
wget https://marinegenomics.oist.jp/aten_2_2_ref/download/AteV2.2Ref_genome.fa.gz
gunzip AteV2.2Ref_genome.fa.gz

./dart2vcf.py -g AteV2.2Ref_genome.fa ../raw_data/Report_DAc23-8597_5_moreOrders_SNP_2.csv > 8597_chr.vcf 2>8597_chr.err

# Reheader because sample 20 is changed to WB_36

bcftools reheader --samples renamed.txt 8597_chr.vcf | bgzip > 8597_chr.vcf.gz
tabix 8597_chr.vcf.gz