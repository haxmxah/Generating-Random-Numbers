c		Marta Xiulan Aribó Herrera
c 		Fem ús del mètode acceptar i rebuig i estudiem la distribució expnencial
c---------------------------COS DEL PROGRAMA PRINCIPAL------------------
        PROGRAM PRE
c--------------------------DEFINICIÓ DE VARIABLES-----------------------
        implicit none
        integer ndat, ndat1, iseed
        double precision fun, pi
        external fun
        parameter (iseed= 20158541)!per crear numeros aleatoris a partir del nostre NIUB
        parameter (pi = 4.d0*atan(1.d0))
        parameter (ndat = 20000)
        parameter (ndat1 = 14000)
        double precision xnums(ndat), a, b, M
        integer nbox, ierr,i, nbox1
        parameter (nbox=30)
        parameter (nbox1=120)
        double precision xhis(nbox),vhis(nbox),errhis(nbox),boxsize
        double precision xhis1(nbox1),vhis1(nbox1),errhis1(nbox1)
        double precision xexpo(ndat1),valor_mitja, var, desv

        open (1, file = "P5-20-21.dat")!on guardarem totes les dades
        write(1,*)"#DADES DE LA PREPRÀCTICA 5 DE NUMEROS ALEATORIS"
        call srand(iseed) !subroutine intrínsica random      
        
c----------------------- APARTAT 0 (opcional) --------------------------
C             no l'he fet jeje, és la subroutien histogram

c----------------------------APARTAT 1 ----------------------------------
C 						MÈTODE D'ACCEPTACIÓ I REBUIG
c 	Es calculen ndat numeros distribuïts uniformement entre [a,b] amb
c aquest mètode, i es retorna la mitjana, la desviació típica i la variança.
c 	A partir de l'arrrya amb ndat numeros, es fa un histograma amb barres d'errors.
c-----------------------------------------------------------------------
        a = -pi !límits
        b = pi
        M = 0.35d0 !cota màxim que més o menys es veu a ull en gràficar p(x)
 		write(1,*) "#Apartat 1 mètode acceptació rebuig"
 		write(1,200)"#mitjana","variació","desviació"
200 	format(a8, 20x, a8, 20x, a10)
        call acceptrebuig(ndat,xnums,a,b,M,fun) 
        call histograma(ndat,xnums,a,b,nbox,xhis,vhis,errhis,boxsize,
     1 ierr)
        write(1,300)"Posicio central","Valor","Interval","Error"
300 	format(a15, 10x, a8, 20x, a10,20x,a5)
        do i = 1, nbox !escriptura de dades
            write(1,*) xhis(i), vhis(i),boxsize,errhis(i)
        enddo
        write(1,"(a1)") !espais en blanc
        write(1,"(a1)") 

        call system ("gnuplot -p plot1.gnu") !graficador de l'histograma

c----------------------------APARTAT 2 ----------------------------------
C 					DISTRIBUCIÓ EXPONENCIAL
c 	Es calculen ndat1 numeros segons la distribució exponencial i es calcula
c 	La mitjana, la desviació típica i la variança, que s'escriuen a l'arxiu
c 	A partir de l'arrrya amb ndat numeros, es fa un histograma amb barres d'errors.
c-----------------------------------------------------------------------
		write(1,*) "#Apartat 2 distribució exponencial"
        call subexpo(ndat1,pi,xexpo)

        !calcul del valor mitjà
        valor_mitja = sum(xexpo)/ndat1
        print*, "El valor mitjà és:", valor_mitja
        !calcul de la variància
        var = 0.d0
        do i = 1,ndat1
              var = var + (xexpo(i)-valor_mitja)**2
        enddo
        var = var/ndat1
        print*, "La variància és:", var
        !càlcul de la desviació estàndard
        desv = sqrt(var)
        print*, "La desviació estàndard és: ", desv

        write(1,200)"#mitjana","variació","desviació"
        write(1,*) valor_mitja,var,desv
        write(1,"(a1)") !espais en blanc
        write(1,"(a1)")

        !valors exactes segons la distribució exponencial
        valor_mitja = 1.d0/pi
        desv = valor_mitja
        var = 1.d0/(pi**2)
        write(1,400)"#mitjana","variació","desviació","(exactes)"
400	    format(a15, 10x, a8, 20x, a10,20x,a5,2x,a15)
        write(1,*) valor_mitja,var,desv
        write(1,"(a1)") !espais en blanc
        write(1,"(a1)")

        a = 0.d0
        b = 3.d0

        call histograma(ndat1,xexpo,a,b,nbox1,xhis1,vhis1,errhis1,
     1  boxsize, ierr)
        write(1,300)"Posicio central","Valor","Interval","Error"
        do i = 1, nbox1
            write(1,*) xhis1(i), vhis1(i),boxsize,errhis1(i)
        enddo
        write(1,"(a1)") !espais en blanc
        write(1,"(a1)") 

        call system ("gnuplot -p plot2.gnu")
        close(1)

        END PROGRAM
c--------------------------- SUBRUTINES I FUNCIONS ---------------------
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
C------------------------------ FUNCIO P(X) ----------------------------
c       Funció de distribució donada pel guió de pràctiques
c-----------------------------------------------------------------------

        double precision function fun(x)
        implicit none
        double precision x, pi 
        parameter (pi = 4.d0*atan(1.d0))
          fun =((5.d0/4.d0)*exp(-abs(x))*sin(x)**2)/(1.d0-exp(-pi)) 
        end function

c--------------------------- SUBRUTINA ACCEPT REBUIG -------------------
c       Subrutina que retorna una llista xnums amb 2000 numeros creats 
c       aleatoriament, i ens retorna la seva mitjana, la seva desviacio
c       i la seva variació
c-----------------------------------------------------------------------

        subroutine acceptrebuig(ndat,xnums,a,b,M,fun) 
            implicit none
            integer ndat, i, comptador, iseed,k
            double precision xnums(ndat),fun,a,b,M,p,x
            double precision valor_mitja, var, desv,x1,x2
            !generem dos numeros aleatoris
            !generador de numeros aleatoris associats al nostre NIUB
            k = 0
            do i = 1,1000
                if (k.NE.ndat) then
                  x1 = rand()
                  x2 = rand()
                  x = ((b-a)*x1) + a
                  p = M*x2
                  if (fun(x).GE.p) then
                  xnums(i) = x
                      k = k + 1
                  endif
                else
                  exit
                endif
            enddo
c            do i = 1,ndat
c                print*, i, "jeje",xnums(i)
c            enddo
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

C---------------------------SUBRUTINA EXPONENCIAL -----------------------
c       Genera un array amb ndat elements amb naturalesa de distribució
c 		exponencial
c-----------------------------------------------------------------------
       subroutine subexpo(ndat,xlambda,xexpo)
       implicit none
       integer ndat, i
       double precision xexpo(ndat), x, xlambda
       do i = 1, ndat
              x= rand()
              xexpo(i) = - (1.d0/xlambda)* log(1.d0-x)
       enddo
c       do i = 1,ndat
c                print*, i, "jeje",xexpo(i)
c       enddo
       end subroutine subexpo








