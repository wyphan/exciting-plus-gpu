
MAKE = make

include make.inc

all:
	cd src; $(MAKE)
	cd utilities/spacegroup; $(MAKE)
	cd utilities/pp; $(MAKE)

clean:
	cd src; $(MAKE) clean
	cd utilities/spacegroup; $(MAKE) clean
	cd utilities/pp; $(MAKE) clean
	rm -rf bin

test:
	cd tests; ./tests.sh

install:
	rm -rf bin; mkdir bin
	ln -s -T ../src/elk bin/elk
	cd utilities/spacegroup; $(MAKE) install
	cd utilities/pp; $(MAKE) install

docs:
	mkdir -p docs
	cd src; $(MAKE) doc; cp elk.pdf ../docs/
	cd ../utilities/spacegroup; $(MAKE) doc; cp spacegroup.pdf ../../docs/
	cp src/addons/CRPA-Calculation.pdf docs/

