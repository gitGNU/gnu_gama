# $Id: Makefile,v 1.2 2002/04/02 21:44:30 cepek Exp $
#
# this Makefile and all files in ./scripts were tested on Debian GNU/Linux 2.2
#

.PHONY : archive


all:
	@if [ ! -f gamaprog/linux/gama/Makefile ] || \
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
	make clean
	./scripts/Build_GaMa $(GaMaLib_CXX)

make-project:
	@if [ ! -f gamaprog/linux/gama/Makefile      ]; then make build; fi
	@if [ ! -f gamaprog/linux/lib/Makefile       ]; then make build; fi
	@if [ ! -f gamaprog/linux/lib/Makefile-expat ]; then make build; fi
	( cd gamaprog/linux/lib;   make; make -f Makefile-expat )
	( cd gamaprog/linux/gama ; make )

example:
	@if [ ! -x  gamaprog/linux/gama/gama ]; then make make; fi
	(cd gamaprog/linux/gama/examples; \
	../gama gama-xml.gkf test. ; cat test.txt )

archive:
	./scripts/GaMaLib_archive


clean:
	rm -f `find gamalib gamaprog -name *\.[o]`
	rm -f `find gamalib gamaprog -name demo`
	rm -f `find . -name *[\~\#]*`
	rm -f `find gamaprog/linux/gama/examples -name test*`

# -----------------------------------------------------------------------

