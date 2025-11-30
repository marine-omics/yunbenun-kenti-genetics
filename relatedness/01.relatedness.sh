# Link in the vcf
# ln -s ../dart2vcf/ak.filtered.vcf.gz .

vcftools --vcf ak.filtered.vcf.gz --relatedness2 --keep maggie_ind.txt --out ak.filtered.maggie_ancestry
vcftools --vcf ak.filtered.vcf.gz --relatedness2 --keep palms_ind.txt --out ak.filtered.palms_ancestry