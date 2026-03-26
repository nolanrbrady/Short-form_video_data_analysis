LATEXMK := latexmk
MAIN := main.tex
LEGACY_ARTIFACTS := \
	main.aux \
	main.bbl \
	main.bcf \
	main.blg \
	main.fdb_latexmk \
	main.fls \
	main.log \
	main.out \
	main.pdf \
	main.run.xml \
	main.synctex.gz

.PHONY: pdf watch clean tidy

pdf:
	mkdir -p build
	$(LATEXMK) -pdf $(MAIN)

watch:
	mkdir -p build
	$(LATEXMK) -pdf -pvc $(MAIN)

clean:
	$(LATEXMK) -C $(MAIN)
	rm -rf build

tidy:
	rm -f $(LEGACY_ARTIFACTS)
