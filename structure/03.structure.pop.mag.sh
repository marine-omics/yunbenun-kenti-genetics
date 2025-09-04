# Full structure analysis on filtered dataset

for rep in $(seq 1 20);do
	mkdir -p mag.${rep}
	cd mag.${rep}
	for k in 1 2 3 4 5 6 7;do
		structure -K ${k} -L 5252 -N 117  -m ../ak.mainparams.pop.mag.nm.txt -i ../ak.pop.mag.nm.str -e ../extraparams -o ak.pop.mag.nm.k${k}.out
		#Extract ancestries
		grep -A 118 'Inferred ancestry of individuals' ak.pop.mag.nm.k${k}.out_f | sed 's/://'> ak.pop.mag.nm.k${k}.out.anc.txt
	done
	# Summarise likelihoods

	grep 'Prob of Data' *.mag.nm.*out* | awk '{print $1,$NF}' | sed -E 's/.*k([1234567]).* (.*)/\1 \2/' > k_vs_lik_mag.txt
	cd ../

done

#Write samples in order for later reading
cat ak.pop.mag.nm.str| awk '{print $1}' | tail -n+2 | uniq > ak.pop.mag.nm.samples.txt

