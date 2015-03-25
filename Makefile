
PYTHON_DIR = "pet"

build:
	python setup.py build_ext --inplace

# docs:
# 	make -C docs html

# install:
# 	python setup.py install --record pyql_install.txt

# uninstall:
# 	cat pyql_install.txt | grep quantlib | xargs rm -rf

# tests-preload:
# 	LD_PRELOAD=/opt/QuantLib-1.1/lib/libQuantLib.so nosetests -v quantlib/test
# tests:
# 	#nosetests -v quantlib/test
# 	cd quantlib/test
# 	python -m unittest discover -v

# build_ex:
# 	g++ -I/opt/local/include/ -I/opt/local/include/boost quantlib_test2.cpp \
#     -o test2 -L/opt/local/lib/ -lQuantLib

clean:
	find ${PYTHON_DIR} -name \*.so -exec rm {} \;
	find ${PYTHON_DIR} -name \*.pyc -exec rm {} \;
	find ${PYTHON_DIR} -name \*.cpp -exec rm {} \;
	find ${PYTHON_DIR} -name \*.c -exec rm {} \;
	rm -rf build
	rm -rf dist

all: build

.PHONY: build # docs clean
