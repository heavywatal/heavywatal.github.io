BUILDDIR := public
DEPLOYDIR := ~/Sites

all:
	$(MAKE) build
	sed -i -e 's/index.html//' ${BUILDDIR}/sitemap.xml
	$(MAKE) install

build:
	hugo -D

install:
	rsync -au --delete --exclude='.git' ${BUILDDIR}/ ${DEPLOYDIR}/

clean:
	rm -r ${BUILDDIR}

server:
	hugo server -D -w
