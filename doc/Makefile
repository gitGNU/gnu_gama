EDITION := $(shell awk '/^@set EDITION/{print $$3}' version.texi)
FILES   := $(shell ls Makefile ChangeLog *.texi *.txt *.eps *.png *.pdf)


all:
	makeinfo gama
	makeinfo --html gama
	cp *.png gama

archive:
	@tar czvf gama-doc-$(EDITION).tar.gz $(FILES)
	@echo
	@ls -sh   gama-doc-$(EDITION).tar.gz
	@echo

clean:
	rm -rf gama
	rm -f  *.aux *.toc *.cp *.fn *.ky *.log *.vr *.pg *.tp
	rm -f  *.dvi *.info gama.pdf
	rm -f  *~
