# Full structure analysis on filtered dataset after removal of clones

run_structure(){
	k=$1
	rep=$2

	mkdir -p nr.${rep}
	cd nr.${rep}

	structure -K ${k} -L 5267 -N 426  -m ../ak.mainparams.nr.txt -i ../ak.filtered.nr.str -e ../extraparams -o ak.structure.nr.k${k}.out
	#Extract ancestries
	grep -A 427 'Inferred ancestry of individuals' ak.structure.nr.k${k}.out_f | sed 's/://'> ak.structure.nr.k${k}.out.anc.txt

	cd ../
}

export -f run_structure

parallel --bar -j 6 run_structure ::: $(seq 1 7) ::: $(seq 1 20)


for rep in $(seq 1 20);do
	cd nr.${rep}

	# Summarise likelihoods
	grep 'Prob of Data' *.nr.*out* | awk '{print $1,$NF}' | sed -E 's/.*k([1234567]).* (.*)/\1 \2/' > k_vs_lik.txt

	cd ../
done

#Write samples in order for later reading
cat ak.filtered.nr.str| awk '{print $1}' | tail -n+2 | uniq > ak.filtered.nr.samples.txt

