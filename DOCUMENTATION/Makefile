default: vector_fit_manual.pdf

LATEX_FILES= \
 Tex/introduction.tex \
 Tex/vector_fit.tex \
 Tex/software_user_guide.tex \
 Tex/examples.tex \
 Tex/bibliography.tex

vector_fit_manual.pdf: vector_fitting.tex  $(LATEX_FILES)
	pdflatex vector_fitting.tex 
	pdflatex vector_fitting.tex 
	pdflatex vector_fitting.tex 

clean:
	rm -f *.aux *.lof *.log  *.out  *.toc

clean_all:
	rm -f *.aux *.lof *.log  *.out  *.toc *.pdf

