@ECHO OFF
SET /A start=8
SET /A end=9
SET /A iterator=1
:while
FOR /L %%i IN (%start%,%iterator%,%end%) DO (
	IF /I %%i GTR %1 (
		GOTO end
	)
	ECHO Iteracion %%i
	FOR /L %%j IN (1,1,1) DO (
		dartSecuential.exe %%i
		dartParallelFor.exe %%i
		dartParallelSchedule.exe %%i
	)
)
IF /I %iterator% GTR 1 (
	SET /A iterator=%iterator%*2
)
SET /A start=%end%+%iterator%
SET /A iterator=%iterator%*10
SET /A end=%iterator%*10-%iterator%
SET /A iterator=%iterator%/2
IF /I %start% LEQ %1 (
	GOTO while
)
:end
ECHO Terminado