TITLE = README

md:
	Rscript -e "knitr::knit('$(TITLE).Rmd')"
