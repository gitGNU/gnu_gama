@rem Win NT:  runs gama for *.gkf files in the current directory
@rem xxxxxxx (version must be explicitly supplied as a parameter)
@rem ------------------------------------------------------------

@if "%1"=="" goto help

@echo @echo %%2             >  for-all-gkf_TMP.bat
@echo @gama %%2 %%2-%%1.txt >> for-all-gkf_TMP.bat

@time /T
@echo -------------------------------------------------
@for %%i in (*.gkf) do @call for-all-gkf_TMP.bat %1 %%i
@echo -------------------------------------------------
@time /T
@del for-all-gkf_TMP.bat
@goto end

:help 
@echo Usage: for-all-gkf.bat gama-version

:end