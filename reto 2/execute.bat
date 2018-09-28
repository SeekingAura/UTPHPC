@ECHO OFF
SET /A start=1000
SET /A end=9000
SET /A iterator=1000
:while
FOR /L %%i IN (%start%,%iterator%,%end%) DO (
	IF /I %%i GTR %2 (
		GOTO end
	)
	ECHO Iteracion %%i
	FOR /L %%j IN (1,1,1) DO (
		testOpenMP.exe %1 %%i
		testSecuential.exe %1 %%i
	)
)
SET /A start=%end%+%iterator%
SET /A iterator=%iterator%*10
SET /A end=%iterator%*10-%iterator%
SET /A iterator=%iterator%
IF /I %start% LEQ %2 (
	GOTO while
)
:end
ECHO Terminado