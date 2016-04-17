BUILDDIR := public
DEPLOYDIR := ~/Sites

build:
	hugo -D && rsync -au --delete --exclude='.git' ${BUILDDIR}/ ${DEPLOYDIR}/

clean:
	rm -rf ${BUILDDIR}

server:
	hugo server -D -w
