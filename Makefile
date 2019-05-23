
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
	mkdir bin
	ln -s -T ../src/elk bin/elk
	cd utilities/spacegroup; $(MAKE) install
	cd utilities/pp; $(MAKE) install
