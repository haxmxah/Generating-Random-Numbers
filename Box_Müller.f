c       MARTA XIULAN ARIBÓ HERRERA
c 		ESTUDI DE DISTRIBUCIONS 
C-----------------------------------------------------------------------

        PROGRAM PRACTICA_5

C----------------------- COS DEL PROGRAMA PRINCIPAL --------------------

C ---------------------- DECLARACIÓ DE VARIABLES    --------------------

        implicit none
        integer iseed, ndat1,nbox1, i, IERR,k, ndat2,nbox2
        parameter (ndat1 = 40000) !apartat 1
        parameter (ndat2 = 10000)
        parameter (nbox1 = 80) !apartat 2
        parameter (nbox2 = 70)
        parameter (iseed = 20158541) !genera variables aleatories associades al NIUB
        double precision p
        external p
        double precision xnums1(ndat1),xnums2(ndat2),a , b, M, L,ro,mu
       double precision xhis1(nbox1),vhis1(nbox1),errhis1(nbox1),boxsize
       double precision xhis2(nbox2),vhis2(nbox2),errhis2(nbox2)
        parameter (ro = 4.d0) !micrometres (es la sigma)
        parameter (L = 3.d0) !nm
        double precision integral, var, valor_mitja,desv
        open (1,file = "P5-20-21-res.dat")
        write(1,*)"#DADES DE LA PREPRÀCTICA 5 DE NUMEROS ALEATORIS"
        call srand(iseed) !inicialització de nombres aleatoris

c--------------------------- APARTAT 1 ---------------------------------

        !APARTAT 1A)
        a = -3*L
        b = +3*L
        M = 0.11d0 !mirem de forma gràfica la funció a graficar
        write(1,*) "#APARTAT 1 ELECTRÓ"
        write(1,200)"#mitjana","variació","desviació"
        call acceptrebuig (ndat1,xnums1,a,b,M,p)
c        print*, xnums1 !comprobació trivial
        call histograma(NDAT1,xnums1,a,b,NBOX1,XHIS1,VHIS1,ERRHIS1,
     1 BOXSIZE,IERR)
        write(1,300)"Posicio central","Valor","Interval","Error"
        do i = 1, nbox1 !escriptura de dades
            write(1,*) xhis1(i), vhis1(i),boxsize,errhis1(i) !guardem les dades al fitxer
        enddo
        write(1,"(a1)") !espais en blanc
        write(1,"(a1)")  
        call system ("gnuplot -p plot1.gnu") !creació de la gràfica

        !APARTAT 1B)
        write (1,*) "CALCUL DE PROBABILITATS I NORMALITZACIÓ"
        k = 12 !intervals
        !comprobació de la normalització de la densitat de probabilitat
        call simpson (a,b,k,p,integral) 
        print*, "Normalització de la integral:", integral
        write(1,"(a23,2x,f20.12)") "#Densitat normalitzada",integral

        !Probabilitat calculada entre els extrems donats
        a = -3.d0*L
        b = L/2.d0
        call simpson (a,b,k,p,integral)
        print*, "Probabilitat:", integral
        write (1,*) "#Resultat de la probabilitat P1"
        write(1,"(a15,2x,f20.12)") "Probabilitat=",integral

        write(1,"(a1)") !espais en blanc
        write(1,"(a1)")

c--------------------------- APARTAT 2 ---------------------------------

c       Hem de crear una llista de números aleatoris segons el mètode de 
c       Box-Muller
 		write(1,*) "#APARTAT 2 ATOM RUBIDI"
 		mu = 0
        call BoxMuller(ndat2,xnums2,mu,ro) !creació de variables amb istrib gaussiana
        a = -3*ro !límits
        b = +3*ro
        call histograma(NDAT2,XNUMS2,A,B,NBOX2,XHIS2,VHIS2,ERRHIS2,
     1 BOXSIZE,IERR) !creació de dades per fer histograma
        write(1,300)"Posicio central","Valor","Interval","Error"
        do i = 1, nbox2 !escriptura de dades
            write(1,*) xhis2(i), vhis2(i),boxsize,errhis2(i) !guardem les dades al fitxer
        enddo
        write(1,"(a1)") !espais en blanc
        write(1,"(a1)")  
        call system ("gnuplot -p plot2.gnu") !creació de la gràfica

c 		Podriem haver calculat dins de la subroutine BoxMuller la mitjana,
c.      la variança i la desviació, però per comoditat a l'hora d'escriure els
c 		fitxers ho farem fora.
		write(1,*) "COMPARACIÓ DE RESULTATS ESTADÍSTICS"
		!calcul de la mitjana
		valor_mitja = sum(xnums2)/ndat2
        print*, "El valor mitjà és:", valor_mitja
        !calcul de la variància
        var = 0.d0
        do i = 1,ndat2
            var = var + (xnums2(i)-valor_mitja)**2
        enddo
        var = var/ndat2
        print*, "La variància és:", var
        !càlcul de la desviació estàndard
        desv = sqrt(var)
        print*, "La desviació estàndard és: ", desv

        !escriptura de dades en el fitxer
        write(1,200)"#mitjana","variació","desviació"
        write (1,*) valor_mitja, var, desv, "Calculades"

        write (1,*) 0.d0, ro**2,ro, "Teoriques"
        write (1,*) "#Comparacio de resultats"
        write (1,*) valor_mitja-0.d0, ro**2-var,ro-desv
        write (1,*) "-------> Resultats raonables"

        write(1,"(a1)") !espais en blanc
        write(1,"(a1)") 

c--------------------------- Format de les variables -------------------

200     format(a8, 20x, a8, 20x, a10)
300     format(a15, 10x, a8, 20x, a10,20x,a5)

        END PROGRAM

c--------------------------- SUBRUTINES I FUNCIONS ---------------------

c--------------------------- BOX MULLER --------------------------------
c      Creació de nombres aleatoris amb el mètode Box Muller 
c-----------------------------------------------------------------------

        subroutine BoxMuller (ndat,xnums,mu,ro)
        implicit none
        integer ndat, i
        double precision xnums(ndat), ro, r ,phi, pi, x, mu
        parameter (pi = 4.d0*atan(1.d0))
        do i = 1, ndat
            r = sqrt(-2.d0*log(rand()))
            phi = 2*pi*rand()
            x = mu + r*sin(phi) !considerem mu= 0
            xnums(i) = x*ro
        enddo
        end subroutine

c--------------------------- FUNCIÓ P(X)--------------------------------
c       Funció unidimensional per un electró en una estructura amb forma 
c       doble-pou 
c-----------------------------------------------------------------------
        double precision function p(x)
        implicit none
        double precision x, L
        parameter (L = 3.d0) !nm

            p = (5.d0/(324.d0*L)) * ((x/L)**2) * (9.d0-(x/L)**2)

        return
        end function

C---------------------------SUBRUTINA HISTOGRAMA -----------------------
c       Subrutina que genera un histograma normalitzat de nbox caixes 
c       fent servir les dades xdata(ndat) generades dins de l'interval
c       [xa,xb]. 
c       Té per output: les xhis(ncaixes) (posicio central de la caixa)
c                      les vhis(ncaixes) (barra corresponent a la ncaixa)
c                      el boxsize (tamany de cada caixa)
c-----------------------------------------------------------------------

       subroutine HISTOGRAMA(NDAT,XDATA,XA,XB,NBOX,XHIS,VHIS,ERRHIS,
     1 BOXSIZE,IERR)

       implicit none
C INPUT/OUTPUT VARIABLES
       integer NDAT,NBOX
       double precision XDATA(NDAT),XA,XB
       double precision XHIS(NBOX),VHIS(NBOX),ERRHIS(NBOX)
       integer IERR
C 
       integer I,IBOX,ICOUNT
       double precision BOXSIZE

       if (XA.GE.XB) then
          IERR=1
          return
       endif
C BOX SIZE
         BOXSIZE=(XB-XA)/NBOX

C COUNTS NUMBER OF POINTS WITHIN THE INTERVAL XA,XB
       ICOUNT=0

C SETS ALL TO ZERO
       do I=1,NBOX
          VHIS(I)=0
          ERRHIS(I)=0
       enddo

C WE RUN THROUGH THE DATASET
       do I=1,NDAT
C CHECKS IF DATA LIES WITHIN XA,XB
         if (XDATA(I).GE.XA.AND.XDATA(I).LE.XB) then
            IBOX=INT((XDATA(I)-XA)/BOXSIZE)+1
C PUTS XB INTO THE LAST BOX, IF NEEDED
            if (IBOX.EQ.NBOX+1) IBOX=NBOX 
            VHIS(IBOX)=VHIS(IBOX)+1
            ICOUNT=ICOUNT+1
         endif
        enddo

        if (ICOUNT.EQ.0) then 
           IERR=2
           RETURN
        endif

        IERR=0
        print*,"ACCEPTED:",ICOUNT,
     1  " OUT OF:",NDAT

        do I=1,NBOX
c CENTRAL VALUE OF THE BAR
           XHIS(I)=XA+BOXSIZE/2.D0+(I-1)*BOXSIZE
C  ERROBAR, STANDARD DEVIATION OF CORRESPONDING BINOMIAL
           ERRHIS(I)=SQRT(VHIS(I)/ICOUNT*(1.D0-VHIS(I)/ICOUNT))
     1      /BOXSIZE / SQRT(DBLE(ICOUNT))
C NORMALIZED VALUE OF THE BAR
           VHIS(I)=VHIS(I)/ICOUNT/BOXSIZE
        enddo
        end subroutine

c--------------------------- SUBRUTINA ACCEPT REBUIG -------------------
c       Subrutina que retorna una llista xnums amb 2000 numeros creats 
c       aleatoriament, i ens retorna la seva mitjana, la seva desviacio
c       i la seva variació
c-----------------------------------------------------------------------

        subroutine acceptrebuig(ndat,xnums,a,b,M,fun) 
            implicit none
            integer ndat, i, comptador, iseed
            double precision xnums(ndat),fun,a,b,M,p,x
            double precision valor_mitja, var, desv,x1,x2
            !generem dos numeros aleatoris
            !generador de numeros aleatoris associats al nostre NIUB
            i = 0
            do while (i.LT.ndat+1)
                x1 = rand()
                x2 = rand()
                x = ((b-a)*x1) + a
                p = M*x2
                if (fun(x).GE.p) then
                xnums(i) = x
                    i = i + 1
                endif
            enddo
            !calcul del valor mitjà
            valor_mitja = sum(xnums)/ndat
            print*, "El valor mitjà és:", valor_mitja
            !calcul de la variància
            var = 0.d0
            do i = 1,ndat
                var = var + (xnums(i)-valor_mitja)**2
            enddo
            var = var/ndat
            print*, "La variància és:", var
            !càlcul de la desviació estàndard
            desv = sqrt(var)
            print*, "La desviació estàndard és: ", desv
            write (1,*) valor_mitja, var, desv !escriptura al fitxer
            write(1,"(a1)") !espais en blanc
            write(1,"(a1)") 
        end subroutine acceptrebuig

c-------------------------SUBRUTINA SIMPSON ----------------------------
c        Subrutine que calcula la integral amb el mètode de simpson repetit
c       Té per entrada els extrems, la funció a integrar i els k de 2^k intervals
c       Té per surtida la integral calculada
c-----------------------------------------------------------------------
        subroutine simpson (x1,x2,k,funci,integral)
        implicit none
        external funci
        double precision x1,x2,integral,h,funci, xk
        integer s,k, interval
        interval = 2**k
        h = (x2-x1)/interval !valors intervals
        integral = (funci(x1)+funci(x2))*h/3.d0 !integral extrems
        do s = 1, 2**k-1
            xk = x1 + s*h
            if (mod(s,2).EQ.0) then 
                integral = integral + (h*2.d0*funci(xk))/3.d0!integral valors interns parells
            else
                integral = integral + (h*4.d0*funci(xk))/3.d0!integral valors interns senars
            endif
        enddo
        end subroutine simpson


