# $Id: Makefile,v 1.14 2003/03/19 19:54:09 cepek Exp $
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
	"\t#" runs make on makefiles created by \`make dep\'
	@echo
	@echo -e make dep \
	"\t#" creates all Makefiles
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

dep:
	./scripts/build-makefiles

build:
	./scripts/build-dictionaries
	./scripts/build-ellipsoids

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

