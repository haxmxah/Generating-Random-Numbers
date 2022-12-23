set term png 
set output "P5-20-21-fig2.png"

set title "Gràfica 2: Àtom de Rubidi"
set xlabel "x (micrometres) " 
set ylabel "Histograma"

set grid xtics
set grid ytics
set key top center

ro = 4.
g(x) = exp(-x**2/(2.*ro**2))/sqrt(2*pi*ro**2)

set style fill solid 0.6 border 

plot "P5-20-21-res.dat" index 3 u 1:2:4 w boxerrorbars lt rgb "pink" t "Histograma w. barres d'errors", "P5-20-21-res.dat" i 3 u 1:2 w histeps t"Histograma", g(x) lt rgb "red"