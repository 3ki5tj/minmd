build::
	python3 build.py

debug::
	python3 build.py debug

clean::
	find . -name "*~" | xargs rm -rf
	rm -rf build/* debug/*

