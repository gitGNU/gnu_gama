EDITION := $(shell awk '/^@set EDITION/{print $$3}' version.texi)
FILES   := $(shell ls Makefile ChangeLog *.texi *.txt *.eps *.png *.pdf)


all:
	makeinfo gama-unstable
	makeinfo --html gama-unstable
	cp *.png gama-unstable

archive:
	@tar czvf gama-unstable-$(EDITION).tar.gz $(FILES)
	@echo
	@ls -sh   gama-unstable-$(EDITION).tar.gz
	@echo

clean:
	rm -rf gama-unstable
	rm -f  *.aux *.toc *.cp *.fn *.ky *.log *.vr *.pg *.tp
	rm -f  *.dvi *.info gama-unstable.pdf
	rm -f  *~
