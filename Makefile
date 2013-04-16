# Makefile copied from SUNDIALS FAQ.
# FIXME: use autoconf to find C compiler and sundials installation.

SHELL = /bin/sh

CC       = gcc
CFLAGS   = -g -O
CPPFLAGS =
LDFLAGS  =

OCAMLOPT = ocamlopt -g

SUN_DISTRO = /dvl/sundials/inst

ML_EXAMPLES = ml-examples/idaRoberts_dns
C_EXAMPLES = c-examples/idaRoberts_dns
ALL_EXAMPLES = ${ML_EXAMPLES} ${C_EXAMPLES}

SUN_CPPFLAGS = `${SUN_DISTRO}/bin/sundials-config -m ida -t s -l c -s cppflags`
SUN_LDFLAGS = `${SUN_DISTRO}/bin/sundials-config -m ida -t s -l c -s libs`

OCAML_SUN_CPPFLAGS=-ccopt "-g ${SUN_CPPFLAGS}"
OCAML_SUN_LDFLAGS=-cclib "-g ${SUN_LDFLAGS}"

all: ${ALL_EXAMPLES}

ml-examples/idaRoberts_dns: IDA.cmxa ml-examples/idaRoberts_dns.ml
	${OCAMLOPT} bigarray.cmxa $^ -o $@ ${OCAML_SUN_CPPFLAGS} ${OCAML_SUN_LDFLAGS}

c-examples/idaRoberts_dns: c-examples/idaRoberts_dns.c
	${CC} ${CPPFLAGS} ${SUN_CPPFLAGS} ${CFLAGS} -c c-examples/idaRoberts_dns.c -o c-examples/idaRoberts_dns.o
	${CC} -o c-examples/idaRoberts_dns c-examples/idaRoberts_dns.o ${CFLAGS} ${LDFLAGS} ${SUN_LDFLAGS}

IDA.cmxa: IDA_Raw.cmx IDA_stubs.o
	${OCAMLOPT} -a -o $@ $^ ${OCAML_SUN_CPPFLAGS} ${OCAML_SUN_LDFLAGS}

IDA_stubs.o: IDA_stubs.c
	${OCAMLOPT} -c ${OCAML_SUN_CPPFLAGS} ${OCAML_SUN_LDFLAGS} -o $@ $<

%.cmx: %.ml
	${OCAMLOPT} -c -o $@ $^

clean:
	rm -f *.o *.cmi *.cmx *.cmxa *.a
	rm -f ${ALL_EXAMPLES}
	rm -f *-examples/*.o *-examples/*.cmo *-examples/*.cmx *-examples/*.cmi
