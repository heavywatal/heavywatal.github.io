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

sass:
	node-sass --output-style expanded -o themes/nonblog/static/css themes/nonblog/src/theme.scss
