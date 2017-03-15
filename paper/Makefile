TITLE="manuscript"

.PHONY: all
all: md tex mdpdf texpdf clean

md:
	pandoc -s $(TITLE).tex -o $(TITLE).md

tex:
	pandoc -s $(TITLE).md -o $(TITLE).tex

pdfmd:
	pandoc $(TITLE).md --latex-engine=xelatex -o $(TITLE).pdf

pdftex:

clean:
