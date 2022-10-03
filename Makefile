.DEFAULT_GOAL := all
.PHONY: all hugo watch clean

all: hugo
	@:

hugo: | public
	hugo

watch:
	hugo -w

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
