# Makefile copied from SUNDIALS FAQ.
# FIXME: use autoconf to find C compiler and sundials installation.

SHELL = /bin/sh
 
CC       = gcc
CFLAGS   = -g -O
CPPFLAGS =
LDFLAGS  =

SUN_DISTRO = /dvl/sundials/inst
 
MY_APPS = c-examples/idaRoberts_dns

all:
	@sun_cppflags=`eval "${SUN_DISTRO}/bin/sundials-config -m ida -t s -l c -s cppflags"`; \
	 sun_ldflags=`eval "${SUN_DISTRO}/bin/sundials-config -m ida -t s -l c -s libs"`;      \
	 for i in ${MY_APPS} ; do                                                                \
	   echo "--- Making $${i} ---" ;                                                         \
	   eval "make SUN_CPPFLAGS='$${sun_cppflags}' SUN_LDFLAGS='$${sun_ldflags}' $${i}";      \
	 done

c-examples/idaRoberts_dns: c-examples/idaRoberts_dns.c
	${CC} ${CPPFLAGS} ${SUN_CPPFLAGS} ${CFLAGS} -c c-examples/idaRoberts_dns.c -o c-examples/idaRoberts_dns.o
	${CC} -o c-examples/idaRoberts_dns c-examples/idaRoberts_dns.o ${CFLAGS} ${LDFLAGS} ${SUN_LDFLAGS}

clean:
	rm -f *.o
	rm -f ${MY_APPS}
