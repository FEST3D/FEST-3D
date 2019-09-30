DIR := ${CURDIR}
all:geometry.f90
	f2py --quiet -c -m geometry geometry.f90
	ln -fs ${CURDIR}/geometry.so SmoothBump/pp/
	ln -fs ${CURDIR}/geometry.so Lfp/pp/
	ln -fs ${CURDIR}/geometry.so Tfp/pp/
