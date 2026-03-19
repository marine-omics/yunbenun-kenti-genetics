# From 
# 11467627|F|0--32:A>T
# To
# 11467627-32-A/T

cat ak.pop.nm.loci | sed 's/-/_/' | sed 's:/:>:' | sed 's/-/:/' | sed 's/_/|F|0--/'  > tmp2.loci 



bcftools filter -i "ID=@tmp2.loci" 8597.1.vcf.gz -Oz > 8597.2.vcf.gz

# Retains 4389 of the original 5251 loci. Some would not have mapped or passed conversion to vcf

bcftools view -S ak.pop.nm.ind 8597.2.vcf.gz -Oz > ak.pop.nm.ind.vcf.gz

# And 403 samples
