# $Id: Makefile,v 1.11 2003/01/06 17:44:12 cepek Exp $
#
# this Makefile and all files in ./scripts were tested on Debian GNU/Linux 2.2
#

.PHONY : archive gamalib


all:
	@if [ ! -f gamaprog/linux/gama-local/Makefile ] || \
	    [ ! -f gamaprog/linux/lib/Makefile  ];   \
	then          \
	   make help; \
	else          \
	   make make-project; \
	fi

help:
	@echo
	@echo -e make "    " \
	"\t#" runs make on makefiles created by \`make build\'
	@echo
	@echo -e make build \
	"\t#" creates all Makefiles, builds GaMaLib and program GaMa
	@echo
	@echo "     aditional" \
	options can be passed to g++ with export CXXFLAGS=\"...\"
	@echo
	@echo -e make example \
	"\t#" runs program gama on an example input data
	@echo -e make archive \
	"\t#" creates tar archive from working directory
	@echo -e make help    \
	"\t#" this screen
	@echo

build:
#	make clean
	./scripts/Build_GaMa

gamalib:
	( cd gamaprog/linux/lib;   make; make -f Makefile-expat )

make-project:
	@if [ ! -f gamaprog/linux/gama-local/Makefile ]; then make build; fi
	@if [ ! -f gamaprog/linux/lib/Makefile        ]; then make build; fi
	@if [ ! -f gamaprog/linux/lib/Makefile-expat  ]; then make build; fi
	( cd gamaprog/linux/lib;   make; make -f Makefile-expat )
	( cd gamaprog/linux/gama-local ; make )

example:
	@if [ ! -x  gamaprog/linux/gama-local/gama-local ]; \
		then make make-project; fi
	(cd gamaprog/linux/gama-local/examples; \
	../gama-local gama-xml.gkf test. ; cat test.txt )

archive:
	./scripts/GaMaLib_archive


clean:
	rm -f `find gamalib gamaprog scripts -name *\.[o]`
	rm -f `find gamalib gamaprog scripts -name *\.[a]`
	rm -f `find gamalib -name demo*`
	rm -f `find . -name *[\~\#]*`
	rm -f `find gamaprog/linux/gama-local/examples -name test*`

# -----------------------------------------------------------------------

