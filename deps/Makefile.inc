# file information

untar_complete = untar.complete
build_complete = build.complete
target_dir = ./

nobase_nodist_mydata_DATA = $(library_file)
mydatadir = $(datadir)/maui

# make rules

all-local: $(dep_files)
	cp $(dep_files) $(target_dir)

$(library_file): $(untar_complete)
	@echo "Building $(software)"
	@cd $(unpack_dir) && if [ -f "./configure" ]; then ./configure $(config_options); fi && $(MAKE)

$(untar_complete): $(tar_file)
	@echo "Extracting $(software) files from tar file"
	@rm -rf $(unpack_dir)
	@tar --extract --$(package_zip_format) --file $(tar_file)
	@touch $(untar_complete)

$(tar_file):
	@echo "Downloading $(software) tar file"
	$(WGET) @WGET_ARGS@ "$(download_url)" > $@

clean-local:
	for file in $(library_file); do \
		if [ -f "$$file" ]; then \
			rm $$file; \
		fi; \
	done
	for file in $(dep_files) $(library_file); do \
		copy=`basename $$file`; \
		if [ -f "$$copy" ]; then \
			rm $$copy; \
		fi; \
	done

maintainer-clean-local:
	for file in Makefile.in $(tar_file) $(untar_complete); do \
		if [ -f "$$file" ]; then \
			rm $$file; \
		fi; \
	done
	rm -rf $(unpack_dir)
