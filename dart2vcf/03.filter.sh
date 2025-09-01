

bcftools filter -i "ID=@tmp.loci" 8597.vcf.gz -Oz > 8597.1.vcf.gz

# Retains 4919 of the original 5276 loci. Some would not have mapped or passed conversion to vcf

bcftools view -S ak.filtered.ind 8597.1.vcf.gz -Oz > ak.filtered.vcf.gz

# And 465 samples
