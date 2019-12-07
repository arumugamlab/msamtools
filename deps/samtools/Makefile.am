## Process this file with automake to produce Makefile.in

# software information

software=samtools
library_file = libbam.a libhts.a

include ../../versions.txt
include defs.txt

#dep_files = $(unpack_dir)/bam.h $(unpack_dir)/bgzf.h $(unpack_dir)/kstring.h $(unpack_dir)/sam.h $(unpack_dir)/libbam.a

dep_files = $(unpack_dir)/libbam.a $(unpack_dir)/htslib-${samtools_version}/libhts.a

include ../Makefile.inc
