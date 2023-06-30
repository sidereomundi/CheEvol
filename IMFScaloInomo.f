      subroutine IMFScaloIno
      real zita
      real Azita,Bzita,CZITA,DZITA
      
      real UM1,UM2,UM3,UM4
      real M1,M3,MS,M4
      real A,B,C,D
      COMMON UM1,UM2,UM3,UM4
      COMMON A,B,C,D
      COMMON M1,M2,M3,M4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      SCALO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    
      UM1=-1.35                 !! esponente stelle piccole  
      UM2=-1.7                  !! esponente stelle massicce 
      M4=0.1                    !! Massa minima              
      M1=2.0                    !! Massa cambio pendenza     
      MS=100.0                  !! Massa massima                                          
      M3=1.                     !! SERVE PER IL CALCOLO DELLA NORMALIZZAZIONE "ZITA"

      ZITA=0.3
      Azita=MS**(1.+UM2)
      Bzita=M1**(1.+UM2)
      Czita=M1**(1.+UM1)
      Dzita=M3**(1.+UM1)        !mi senza zita o m3 se vuoi usare zita a 0.3!!!


      B=ZITA/(M1**(UM2-UM1)*(Czita-Dzita)/(1.+UM1)+  ! integrale da 1 a 100 dev'essere 0.3
     $     (Azita-Bzita)/(1.+UM2))

      A=B*M1**(UM2-UM1)                              ! La funzione dev'essere continua in 2Msun



      Azita=MS**(UM2)
      Bzita=M1**(UM2)
      Czita=M1**(UM1)
      Dzita=M4**(UM1)  !!! Nota cambio, sto prendendo 0.1 come estremo inferiore


      
      
      C = A* (CZITA-DZITA)/(UM1)  ! trambo
      D = B* (AZITA-BZITA)/(UM2)  ! qrambo
      
      write(*,*) a,b,c,d,um1,um2

      
      return

      end

