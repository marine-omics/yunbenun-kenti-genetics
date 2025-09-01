
structure -K 2 -L 5275 -N 465  -m ak.mainparams.txt -i ak.filtered.str -o ak.structure.k2.out

#Extract ancestries

grep -A 466 'Inferred ancestry of individuals' ak.structure.k2.out_f | sed 's/://'> ak.structure.k2.out.anc.txt

#Write samples in order for later reading
cat ak.filtered.str| awk '{print $1}' | tail -n+2 | uniq > ak.filtered.samples.txt
