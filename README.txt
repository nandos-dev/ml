
Preprocessin::using notepad++ to remove all information besides the gene-name, utilizing the regex r" .*"
	" " matches the string literally (case sensitive)
	.   matchs any character (except for line terminators)
	*   matches the previous token between zero and unlimited times, as many times as possible, giving 
	back as needed(via greedy) 

IN MATLAB 
	regex was used in matlab to import the txt data of genes from the pdf reading 
	(http://www.ihes.fr/~zinovyev/princmanif2006/Wang_lancet_2005.pdf)
	and into cell arrays. 