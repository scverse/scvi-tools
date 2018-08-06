nbs = $(wildcard *.ipynb)
pdfs = $(nbs:%.ipynb=%.pdf)

all: $(pdfs)

%.pdf: %.ipynb
	jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=1800 $<;
	jupyter nbconvert --to pdf $<;
