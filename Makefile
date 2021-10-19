build::
	python3 build.py

clean::
	find . -name "*~" | xargs rm -rf
	rstrip.py -vR
