BUILDDIR := public
DEPLOYDIR := ~/Sites

all:
	$(MAKE) build
	$(MAKE) install

build:
	hugo -D

install:
	rsync -au --delete --exclude='.git' ${BUILDDIR}/ ${DEPLOYDIR}/

clean:
	rm -rf ${BUILDDIR}

server:
	hugo server -D -w
