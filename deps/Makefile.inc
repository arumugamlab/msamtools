# file information

untar_complete = untar.complete
build_complete = build.complete
dep_complete = dep.complete
target_dir = ./

# make rules

all-local: $(dep_complete)
	cp $(dep_files) $(target_dir)

$(dep_complete): $(untar_complete)
	@echo "Building $(software)"
	@cd $(unpack_dir) && if [ -f "./configure" ]; then ./configure $(config_options); fi && $(MAKE)
	@touch $(dep_complete)

$(untar_complete): $(tar_file)
	@echo "Extracting $(software) files from tar file"
	@rm -rf $(unpack_dir)
	@tar --extract --$(package_zip_format) --file $(tar_file)
	@touch $(untar_complete)

$(tar_file):
	@echo "Downloading $(software) tar file"
	$(WGET) @WGET_ARGS@ "$(download_url)" > $@

clean-local:
	for file in $(dep_files) $(library_file); do \
		copy=`basename $$file`; \
		if [ -f "$$copy" ]; then \
			rm $$copy; \
		fi; \
	done

distclean-local:
	for file in $(tar_file) $(untar_complete) $(dep_complete); do \
		if [ -f "$$file" ]; then \
			rm $$file; \
		fi; \
	done
	rm -rf $(unpack_dir)

maintainer-clean-local:
	@echo 'This command is intended for maintainers to use; it'
	@echo 'deletes files that may need special tools to rebuild.'
	for file in Makefile.in; do \
		if [ -f "$$file" ]; then \
			rm $$file; \
		fi; \
	done
