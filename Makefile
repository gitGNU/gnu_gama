# $Id: Makefile,v 1.22 2004/06/21 16:10:13 cepek Exp $
#
# this Makefile and all files in ./scripts were tested on Debian GNU/Linux 3.0
#


.PHONY : archive


all:
	( cd gamaprog/linux/lib; make -f Makefile-expat; make )
	( cd gamaprog/linux/gama-local ; make )
	( cd gamaprog/linux/gama-g3 ; make )

help:
	@echo
	@echo -e make "    " \
	"\t#" run make on makefiles created by \`make dep\'
	@echo
	@echo -e make dep \
	"\t#" create all Makefiles
	@echo
	@echo "     aditional" \
	options can be passed to g++ with export CXXFLAGS=\"...\"
	@echo
	@echo -e make dep-expat-1.1\
	"\t#"  create alternative Makefiles 
	@echo -e "           "\
	"\t\t#    for old version 1.1 of expat XML parser"
	@echo -e make example \
	"\t\t#" runs program gama on an example input data
	@echo -e make archive \
	"\t\t#" creates tar archive from working directory
	@echo -e make help    \
	"\t\t#" this screen
	@echo

dep:
	./scripts/build-makefiles

dep-expat-1.1:
	./scripts/build-makefiles-expat-1.1

build:
	./scripts/build-dictionaries
	./scripts/build-ellipsoids

example:
	@if [ ! -x  gamaprog/linux/gama-local/gama-local ]; \
		then make make-project; fi
	(cd gamaprog/linux/gama-local/examples; \
	../gama-local gama-xml.gkf test. ; cat test.txt )

archive:
	./scripts/build-archive

clean:
	rm -f `find gnu_gama gamalib gamaprog scripts -name *\.[o]`
	rm -f `find gnu_gama gamalib gamaprog scripts -name *\.[a]`
	rm -f `find gnu_gama gamalib -name leak.out`
	rm -f `find gnu_gama gamalib -name demo*`
	rm -f `find gamaprog/linux/gama-local/examples -name test*`
	rm -f `find gamaprog/linux/gama-g3/examples -name test*`
#	rm -f `find . -name "*[\~\#]*"`

clean-scripts:
	make clean
	(cd scripts; rm -f gnu_gama_dep slovnikar ellipsoids_xml \
		utf8-ascii )

# -----------------------------------------------------------------------

