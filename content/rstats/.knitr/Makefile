DESTINATION := ..
SOURCES := $(wildcard *.Rmd)
KNITTED := $(addprefix ${DESTINATION}/,$(SOURCES:.Rmd=.md))

.DEFAULT_GOAL := all
.PHONY: all clean

all: ${KNITTED}
	@:

${DESTINATION}/%.md: %.Rmd setup.R
	Rscript --no-init-file -e "wtl::exec_knit('$<', '$@')"
# 	Rscript --no-init-file -e "litedown::fuse('$<', '$@')"

clean:
	$(RM) -r ${KNITTED} ../.cache ../*__files ${DESTINATION}/figure
