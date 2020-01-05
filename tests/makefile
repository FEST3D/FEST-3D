DIR := ${CURDIR}
all:geometry.f90
	f2py3 --quiet -c -m geometry geometry.f90
	bash -c "mv  geometry.cpython*.so geometry.so"
	ln -fs ${CURDIR}/geometry.so SmoothBump/pp/
	ln -fs ${CURDIR}/geometry.so Lfp/pp/
	ln -fs ${CURDIR}/geometry.so Tfp/pp/
