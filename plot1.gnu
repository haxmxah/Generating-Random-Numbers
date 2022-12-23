set term png 
set output "P5-20-21-fig1.png"

set title "Gràfica 1: Electró en una estructura unidimensional"
set xlabel "x (nm) " 
set ylabel "Histograma"

set grid xtics
set grid ytics
set key top center

L = 3.
p(x) = (5./(324.*L)) * (x/L)**2 * (9.-(x/L)**2)


set style fill solid 0.4 border 

plot "P5-20-21-res.dat" index 1 u 1:2:4 w boxerrorbars lt rgb "orange" t "Histograma w. barres d'errors", "P5-20-21-res.dat" i 1 u 1:2 w histeps t"Histograma", p(x) t"P(x)" lt rgb "red"