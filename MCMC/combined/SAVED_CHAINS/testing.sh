#!/bin/bash

echo "\documentclass[10pt]{article}
\usepackage{algorithmic}
\usepackage{algorithm}
\usepackage{amsmath}
\usepackage[T1]{fontenc}
\usepackage{multicol}
\usepackage{graphicx}
\usepackage[%  
    colorlinks=true,
    pdfborder={0 0 0},
    linkcolor=red
]{hyperref}

\begin{document}

\listoffigures" >> temp.tex

for S in 85 76 96 81 80 36 83 84 94 32 106 123 50 39 56 126 54 70 88 67 95 25 137 1 86 43 73 92 93 79 47 109 44 48 17 8 12 20 31 57 120 28 13
    do
        ./corner_plot.sh $S

		cat ganymede/MCMC_$S/accepted_parameters_*.txt >ganymede/MCMC_$S/combined_accepted.txt
        cat stampede/MCMC_$S/accepted_parameters_*.txt >stampede/MCMC_$S/combined_accepted.txt
		cat ganymede/MCMC_$S/combined_accepted.txt stampede/MCMC_$S/combined_accepted.txt >combined_accepted.txt
        W="$(wc -l combined_accepted.txt)"
        A="$(cut -d' ' -f1<<<$W)" 
        rm ganymede/MCMC_$S/combined_accepted.txt
        rm stampede/MCMC_$S/combined_accepted.txt
        rm combined_accepted.txt
		
        cat ganymede/MCMC_$S/rejected_parameters_*.txt >ganymede/MCMC_$S/combined_rejected.txt
        cat stampede/MCMC_$S/rejected_parameters_*.txt >stampede/MCMC_$S/combined_rejected.txt
        cat ganymede/MCMC_$S/combined_rejected.txt stampede/MCMC_$S/combined_rejected.txt >combined_rejected.txt
        W="$(wc -l combined_rejected.txt)"
        R="$(cut -d' ' -f1<<<$W)"
        rm ganymede/MCMC_$S/combined_rejected.txt
        rm stampede/MCMC_$S/combined_rejected.txt
        rm combined_rejected.txt
        
        sum=$(($A+$R))
        
        r=$(bc <<<"scale=2;$A/$sum")
        
        T="System=$S|Accepted=${A}|Rejected=${R}|AcceptanceRatio=$r"
	    echo "$T" >> temp.tex
    
        echo "\begin{figure}[h]" >> temp.tex 
        echo "\includegraphics[width=15cm]{CornerPlot_$S.eps}" >> temp.tex
        echo "\caption{$T}" >> temp.tex
        echo "\label{S$S}" >> temp.tex
        echo "\centering" >> temp.tex
        echo "\end{figure}" >> temp.tex
        echo "" >> temp.tex



done
