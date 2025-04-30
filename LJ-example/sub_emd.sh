#!/usr/bin/

function rand(){

	min=$1
	max=$(($2-$min+1))
	num=$(($RANDOM+1000000000))
	echo $(($num%$max+$min))
}

for i in {1..20};
do
	mkdir ${i}
	cp EMD_enthalpy.lmp ${i}
        cp 1.sh ${i}
	cd ${i}

	r=$(rand 100000 1000000)
        T=400

	echo $r
	sed -i "s/tobeset_2/${r}/g" EMD_enthalpy.lmp
	sed -i "s/tobeset_1/${T}/g" EMD_enthalpy.lmp
	echo "the calculate temperature is ${T}K"
        sbatch 1.sh
	cd ..
done
