set term png 
set output "P20_1.png"

set title "Mètode d'acceptació i rebuig"
set xlabel "x" 
set ylabel "P(x)"

set grid xtics
set grid ytics

L = 4.
p(x)= sin(x/L)**2/(L*pi)


set style fill solid 0.4 border 

plot "P20.dat" index 1 u 1:2:4 w boxerrorbars lt rgb "blue" t "barres d'errors", "P20.dat" i 1 u 1:2 w histeps t"histograma",p(x)