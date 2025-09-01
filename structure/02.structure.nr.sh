# Full structure analysis on filtered dataset

for k in 1 2 3 4;do
	structure -K ${k} -L 5269 -N 436  -m ak.mainparams.nr.txt -i ak.filtered.nr.str -o ak.structure.nr.k${k}.out
	#Extract ancestries
	grep -A 437 'Inferred ancestry of individuals' ak.structure.nr.k${k}.out_f | sed 's/://'> ak.structure.nr.k${k}.out.anc.txt
done

#Write samples in order for later reading
cat ak.filtered.nr.str| awk '{print $1}' | tail -n+2 | uniq > ak.filtered.samples.txt
