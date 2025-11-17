
run_structure(){
	k=$1
	rep=$2

	mkdir -p mag.${rep}
	cd mag.${rep}
	structure -K ${k} -L 2939 -N 117  -m ../ak.mainparams.pop.mag.nm.txt -i ../ak.pop.mag.nm.str -e ../extraparams -o ak.pop.mag.nm.k${k}.out
	#Extract ancestries
	grep -A 118 'Inferred ancestry of individuals' ak.pop.mag.nm.k${k}.out_f | sed 's/://'> ak.pop.mag.nm.k${k}.out.anc.txt

	cd ../
}

export -f run_structure

parallel --bar -j 8 run_structure ::: $(seq 1 7) ::: $(seq 1 20)





for rep in $(seq 1 20);do

	cd mag.${rep}

	# Summarise likelihoods
	grep 'Prob of Data' *.mag.nm.*out* | awk '{print $1,$NF}' | sed -E 's/.*k([1234567]).* (.*)/\1 \2/' > k_vs_lik_mag.txt
	cd ../

done

#Write samples in order for later reading
cat ak.pop.mag.nm.str| awk '{print $1}' | tail -n+2 | uniq > ak.pop.mag.nm.samples.txt

