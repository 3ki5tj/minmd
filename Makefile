build::
	python3 build.py

debug::
	python3 build.py debug

nol::
	# count the number of lines (including only .c and .h files)
	git ls-files | grep -E "\.c|\.h" | grep -v notes | xargs wc | tail -n 1
	# count the number of lines (including all files)
	git ls-files | xargs wc | tail -n 1

clean::
	find . -name "*~" | xargs rm -rf
	rm -rf build/* debug/*

