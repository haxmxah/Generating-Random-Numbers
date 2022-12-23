set term png 
set output "P20_2.png"

set title "Mètode d'acceptació i rebuig"
set xlabel "x" 
set ylabel "P(x)"

set grid xtics
set grid ytics

ro = 3.
g(x)= (6./(pi*ro*(2*pi**2-3)))*(x/ro)**2*sin(x/ro)**2

set style fill solid 0.4 border 

plot "P20.dat" index 3 u 1:2:4 w boxerrorbars lt rgb "red" t "barres d'errors", "P20.dat" i 3 u 1:2 w histeps t"histograma",g(x)