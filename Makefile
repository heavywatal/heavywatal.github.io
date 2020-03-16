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
	ls -d public/* | grep -Ev 'offline|slides|hpc-' | xargs $(RM) -r
