FOR /L %%i IN (1,1,%1) DO (
    python "generate input data.py" %%i
    FOR /L %%j IN (1,1,10) DO (
        python "matrices-Secuencial.py"
        python "matrices-Paralelo.py"
    )
)