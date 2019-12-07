DIST_SUBDIRS = argtable samtools

SUBDIRS = argtable samtools

#all-local: libargtable2.a argtable2.h libbam.a bam.h bgzf.h kstring.h sam.h

all-local: libargtable2.a argtable2.h libbam.a libhts.a

libargtable2.a argtable2.h:
	@cp argtable/$@ .

#libbam.a bam.h bgzf.h kstring.h sam.h:

libbam.a libhts.a:
	@cp samtools/$@ .

clean-local:
	for file in libargtable2.a argtable2.h libbam.a libhts.a; do \
		if [ -f "$$file" ]; then \
			rm $$file; \
		fi; \
	done
