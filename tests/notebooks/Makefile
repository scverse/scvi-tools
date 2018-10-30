type ?= html

nbs = $(wildcard *.ipynb)
types = $(nbs:%.ipynb=%.$(type))

all: $(types)

%.$(type): %.ipynb
	jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=1800 $<;
	jupyter nbconvert --to $(type) $<;
