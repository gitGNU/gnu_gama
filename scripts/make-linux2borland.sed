# basic sed commands for linux to borland makefiles conversion
#
# $Id: make-linux2borland.sed,v 1.2 2001/12/20 19:49:43 cepek Exp $

s/\.o/.obj/g
s/^CPP=.*/CPP=bcc32/
s/^CC=.*/CC=bcc32/
s/-Wall.*-pedantic/-w-8026 -w-8027 -w-8004 -tWR #-tWM/  
s/gamalib\.a/gamalib.lib/
s/.*touch.*//
s/.*ranlib.*//
s/LIBR=.*/LIBR=tlib \/P4096 gamalib.lib +/
s/^gama : /gama.exe : /
s/-o gama/-Egama.exe/
s/ : gama.cpp.*/ : ..\/..\/linux\/gama\/gama.cpp/
s/-I. -c gama.cpp/-I..\/..\/linux\/gama -c ..\/..\/linux\/gama\/gama.cpp/

# additional parameters needed for parser expat

s/-DXML_NS/-A -Od -w-8008 -w-8065 -w-8066 -w-8057 -DXML_NS/











