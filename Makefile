BUILDDIR := public

all:
	$(MAKE) build
	sed -i -e 's/index.html//' ${BUILDDIR}/sitemap.xml

build:
	hugo -D

clean:
	ls -d ${BUILDDIR}/* | grep -v .git | xargs rm -r

server:
	hugo server -D -w
