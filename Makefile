.DEFAULT_GOAL := all
SUBDIRS := $(wildcard content/*/.knitr)
.PHONY: all hugo watch benchmark clean ${SUBDIRS}

all: hugo
	@:

hugo: ${SUBDIRS} | public
	hugo

${SUBDIRS}:
	$(MAKE) -C $@

watch:
	hugo -w

benchmark:
	hugo --templateMetrics --templateMetricsHints

public:
	$(eval REMOTE_URL := $(shell git remote get-url origin))
	git clone -b master --single-branch ${REMOTE_URL} $@

clean:
	ls -d public/* | grep -Ev 'offline|slides|hpc-|jbrowse' | xargs $(RM) -r

.PHONY: hex-stickers

HEXSRC := submodules/hex-stickers
HEXDST := static/_img/hex-stickers
TIDYVERSE := covr devtools dplyr forcats ggplot2 knitr lubridate pipe pkgdown purrr quarto readr readxl rmarkdown roxygen2 stringr testthat tibble tidyr tidyverse usethis
TIDYVERSE := $(addprefix ${HEXDST}/, ${TIDYVERSE})
TIDYVERSE := $(addsuffix .webp, ${TIDYVERSE})

hex-stickers: ${TIDYVERSE}
	@:

${HEXDST}/%.webp: ${HEXSRC}/PNG/%.png | ${HEXDST}
	cwebp -lossless $< -o $@

${HEXDST}:
	mkdir -p ${HEXDST}
