# basic sed commands for Borland format to Visual C++  makefiles conversion
#
# $Id: make-borland2msc.sed,v 1.1 2001/12/07 11:45:44 cepek Exp $

s;^CPP=.*;OUTDIR=.\
INTDIR=.\
OutDir=.\
\
CPP=cl.exe\
CPP_PROJ=/nologo /W1 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /Zp1 /MT /GR\
\
LIB32=link.exe -lib\
LIB32_FLAGS=gamalib.lib\
;
s;^ALL : gamalib.lib;ALL : init gamalib.lib\
\
init:\
        @if not exist gamalib.lib cl.exe -I../../.. -c ../../../gamalib/version.cpp\
        @if not exist gamalib.lib link -lib /out:gamalib.lib version.obj\
        @if exist version.obj     del version.obj;
s/^CC=.*/CC=cl.exe/
s;^FCC=-I../../.. .*-w-[0-9]*;FCC=-I../../.. \$(CPP_PROJ);
s;^FCC=-O2 -I../../../expat.*;FCC=-I../../../expat/xmltok -I../../../expat/xmlparse /nologo /ML /W1 /GX /O2 /D "NDEBUG" /D "WIN32" /D "_CONSOLE" /YX;
s/make/nmake/
s/^LIBR=.*/LIBR=\$(LIB32) \$(LIB32_FLAGS)/
s/\$(CPP) -Egama.exe/\$(LINK)/
s;gama.exe : gama.obj.*;LINK=link kernel32.lib user32.lib gdi32.lib\\\
winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib\\\
oleaut32.lib uuid.lib odbc32.lib odbccp32.lib\\\
/nodefaultlib:libc /nologo /subsystem:console\\\
/incremental:no /machine:I386 /out:"gama.exe"\
\
gama.exe : gama.obj $(GAMALIB);







