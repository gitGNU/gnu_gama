EDITION := $(shell awk '/^@set EDITION/{print $$3}' version.texi)
FILES   := $(shell ls Makefile ChangeLog *.texi *.txt *.eps *.png *.pdf)


all:
	make info

info:
	makeinfo gama

html:
	rm -rf gama
	makeinfo --html gama
	cp *.png gama

dvi:
	tex gama.texi
	texindex gama.??
	tex gama.texi
	dvips gama

pdf:
	pdftex gama.texi
	texindex gama.??
	pdftex gama.texi

archive:
	@tar czvf gama-doc-$(EDITION).tar.gz $(FILES)
	@echo
	@ls -sh   gama-doc-$(EDITION).tar.gz
	@echo

clean:
	rm -f  *.aux *.toc *.cp *.fn *.ky *.log *.vr *.pg *.tp
	rm -f  *~
