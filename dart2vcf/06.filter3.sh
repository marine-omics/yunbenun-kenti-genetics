cat ak_r1.filtered.loci | sed 's/-/_/' | sed 's:/:>:' | sed 's/-/:/' | sed 's/_/|F|0--/'  > tmp2.loci 
bcftools filter -i "ID=@tmp2.loci" 8597_chr.vcf.gz -Oz > 8597_chr.1.vcf.gz

# Retains 4942 of 5325 of the original loci. Some would not have mapped or passed conversion to vcf

bcftools view -c 1 -S ak_r1.filtered.ind 8597_chr.1.vcf.gz -Oz > ak_r1.pop.nm.ind.vcf.gz

