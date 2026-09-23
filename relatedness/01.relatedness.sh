# Link in the vcf
# ln -s ../dart2vcf/ak.filtered.vcf.gz .

vcftools --gzvcf ak.filtered.vcf.gz --relatedness2 --keep maggie_ind.txt --out ak.filtered.maggie_ancestry
vcftools --gzvcf ak.filtered.vcf.gz --relatedness2 --keep palms_ind.txt --out ak.filtered.palms_ancestry

vcftools --gzvcf ak.filtered.vcf.gz --relatedness2 --out ak.filtered.all


ln -s ../dart2vcf/ak_r1.pop.nm.ind.vcf.gz .
vcftools --gzvcf ak_r1.pop.nm.ind.vcf.gz --relatedness2 --out ak_r1.filtered.all
