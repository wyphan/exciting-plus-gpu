
MAKE = make
UTILS := spacegroup pp pp_u dx2silo
include make.inc

all:
	cd src; $(MAKE)
	for UTIL in $(UTILS); do cd utilities/$$UTIL; $(MAKE); cd ../..; done

clean:
	cd src; $(MAKE) clean
	for UTIL in $(UTILS); do cd utilities/$$UTIL; $(MAKE) clean; cd ../..; done
	rm -rf bin

test:
	cd tests; ./tests.sh

install:
	rm -rf bin; mkdir bin
	ln -s -T ../src/elk bin/elk
	for UTIL in $(UTILS); do cd utilities/$$UTIL; $(MAKE) install; cd ../..; done

docs:
	mkdir -p docs
	cd src; $(MAKE) doc; cp elk.pdf ../docs/
	cd ../utilities/spacegroup; $(MAKE) doc; cp spacegroup.pdf ../../docs/
	cp src/addons/CRPA-Calculation.pdf docs/
