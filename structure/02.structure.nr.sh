# Full structure analysis on filtered dataset

for rep in $(seq 1 20);do
	mkdir -p nr.${rep}
	cd nr.${rep}
	for k in 1 2 3 4 5 6 7;do
		structure -K ${k} -L 5265 -N 427  -m ../ak.mainparams.nr.txt -i ../ak.filtered.nr.str -e ../extraparams -o ak.structure.nr.k${k}.out
		#Extract ancestries
		grep -A 428 'Inferred ancestry of individuals' ak.structure.nr.k${k}.out_f | sed 's/://'> ak.structure.nr.k${k}.out.anc.txt
	done
	# Summarise likelihoods

	grep 'Prob of Data' *.nr.*out* | awk '{print $1,$NF}' | sed -E 's/.*k([1234567]).* (.*)/\1 \2/' > k_vs_lik.txt
	cd ../
done

#Write samples in order for later reading
cat ak.filtered.nr.str| awk '{print $1}' | tail -n+2 | uniq > ak.filtered.nr.samples.txt

