.DEFAULT_GOAL := all
SUBDIRS := $(wildcard content/*/.knitr)
.PHONY: all development public watch benchmark hugo clean ${SUBDIRS}

all: development public
	@:

development public: ${SUBDIRS}
	hugo --environment $@

public-clean: ${SUBDIRS}
	hugo --environment public --cleanDestinationDir

${SUBDIRS}:
	$(MAKE) -C $@

watch:
	hugo -w

benchmark:
	hugo --templateMetrics --templateMetricsHints

pull:
	git pull
	git -C public pull

init:
	git clone -b master --single-branch $$(git remote get-url origin) public
	git submodule update --init --single-branch

echo:
	echo $$(git remote get-url origin)

clean:
	ls -d development/* | grep -Ev 'offline|slides|hpc-|jbrowse' | xargs $(RM) -r

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
