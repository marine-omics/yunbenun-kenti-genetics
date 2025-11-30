# From 
# 11467627|F|0--32:A>T
# To
# 11467627-32-A/T

cat ak.filtered.loci| sed 's/-/_/' | sed 's:/:>:' | sed 's/-/:/' | sed 's/_/|F|0--/'  > tmp.loci 

