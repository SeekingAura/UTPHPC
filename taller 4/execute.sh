start=8
end=9
iterator=1
while [$start -lq $@] do
    
    for ((i = 2; i <= $@; i++ )); do
        if [$i -gt $@] then
		    break
	    fi
        echo "Iteracion $i"

        for ((j = 1; j <= 1; j++ )); do
            ./dartSecuential.bin "$i"
            ./dartParallel.bin "$i"
            ./dartParallelAtomic.bin "$i"
            ./dartParallelCritical.bin "$i"
            ./dartParallelFor.bin "$i"
            ./dartParallelReduction.bin "$i"
            ./dartParallelSchedule.bin "$i"
        done
    done
    # if [$iterator -gt 1] then
    #	let iterator=iterator*2
    #)
    let start=end+iterator
    let iterator=iterator*10
    let end=iterator*10-iterator
    # let iterator=iterator/2
    if [$start -lt $@ || $start -eq $@] then
        break
    fi
done