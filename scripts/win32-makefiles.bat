@rem a) compile gnu_gama_dep 
@rem =======================
@rem
@rem    bcc32 gnu_gama_dep.cpp   
@rem 
@rem    or
@rem
@rem    cl /EHsc gnu_gama_dep.cpp
@rem
@rem
@rem b) run script 'win32-makefiles' in script directory as
@rem ====================================================== 
@rem
@rem    win32-makefiles.bat  win32-borland
@rem
@rem    or
@rem
@rem    win32-makefiles.bat  win32-msvc
@rem

cd ..

scripts\gnu_gama_dep lib        %1 > gamaprog\%1\lib\Makefile
scripts\gnu_gama_dep expat      %1 > gamaprog\%1\lib\Makefile-expat
scripts\gnu_gama_dep gama-local %1 > gamaprog\%1\gama-local\Makefile

cd scripts

