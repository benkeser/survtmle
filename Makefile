TITLE = README

md:
	Rscript -e "rmarkdown::render('$(TITLE).Rmd')"
