SUBDIRS = $(shell ls -d examples/*)

all:
	for dir in $(SUBDIRS) ; do \
		make -C  $$dir ; \
	done

clean:
	for dir in $(SUBDIRS) ; do \
		make -C  $$dir clean; \
	done

run:
	for dir in $(SUBDIRS) ; do \
		make -C  $$dir run; \
	done

mem-debug:
	for dir in $(SUBDIRS) ; do \
		make -C  $$dir mem-debug; \
	done

pmem-debug:
	for dir in $(SUBDIRS) ; do \
		make -C  $$dir pmem-debug; \
	done