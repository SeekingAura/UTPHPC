for ((i = 2; i <= $@; i++ )); do
    python3 "generate input data.py" $i
    for ((j = 1; j <= 1; j++ )); do
        python3 "matrices-Secuencial.py"
        python3 "matrices-ParaleloHilos.py"
        python3 "matrices-ParaleloProcesos.py"
    done
done