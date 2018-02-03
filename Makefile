md:
	Rscript -e "rmarkdown::render('README.Rmd', output_file = 'README.md')"

site:
	Rscript -e "pkgdown::build_site()"

fastcheck:
	Rscript -e "devtools::check(build_args = '--no-build-vignettes')"

check:
	Rscript -e "devtools::check()"

test:
	Rscript -e "devtools::test()"

doc:
	Rscript -e "devtools::document()"

cov:
	Rscipt -e "covr::package_coverage(type = 'all', combine_types = FALSE, line_exclusions = list('R/plots.R', 'R/printing.R', 'R/utils.R'))"

build:
		Rscript -e "devtools::build()"

style:
		Rscript -e "styler::style_pkg()"

