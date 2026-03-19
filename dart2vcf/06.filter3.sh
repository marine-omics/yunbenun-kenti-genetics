cat ak.pop.nm.loci | sed 's/-/_/' | sed 's:/:>:' | sed 's/-/:/' | sed 's/_/|F|0--/'  > tmp2.loci 

bcftools filter -i "ID=@tmp2.loci" 8597_chr.vcf.gz -Oz > 8597_chr.1.vcf.gz

# Retains 3801 of 4107 of the original loci. Some would not have mapped or passed conversion to vcf

bcftools view -S ak.pop.nm.ind 8597_chr.1.vcf.gz -Oz > ak.pop.nm.ind.vcf.gz

# And 403 samples