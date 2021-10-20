build::
	python3 build.py

clean::
	find . -name "*~" | xargs rm -rf
	rm -rf build/*

