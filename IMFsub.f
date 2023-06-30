      subroutine IMFKroupa

      real UM1,UM2,UM3,UM4
      real MI,M1,M2,M3,MS
      real A1,A2,B1,B2,B3
      real C1,C2,C3,D1,D2
      real A,B,C,D

      COMMON UM1,UM2,UM3,UM4
      COMMON A,B,C,D
      COMMON M1,M2,M3,M4


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      KROUPA 2001
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      UM1= 0.7                  !! esponente stelle M<0.08
      UM2=-0.8                  !! esponente stelle 0.08<M<0.5
      UM3=-1.7                  !! esponente stelle 0.5<M<1
      UM4=-1.3                  !! esponente stelle M<1

      MI=0.05                   !! Massa minima              
      M1=0.08 
      M2=0.5   
      M3=1.0                    !! Massa cambio pendenza     
      MS=50.0                   !! Massa massima                                          
      

      A1=(M1**(UM1+1))/(UM1+1)
      A2=(MI**(UM1+1))/(UM1+1)
      
      B1=(M2**(UM2+1))/(UM2+1)
      B2=(M1**(UM2+1))/(UM2+1)

      C1=(M3**(UM3+1))/(UM3+1)
      C2=(M2**(UM3+1))/(UM3+1)

      D1=(MS**(UM4+1))/(UM4+1)
      D2=(M3**(UM4+1))/(UM4+1)

      B3=M1**(UM1-UM2)

      C3=M2**(UM2-UM3)*B3
      

      A=1./(A1-A2+B3*(B1-B2)+C3*(C1-C2)+C3*(D1-D2))

      B=A*B3

      C=A*C3

      D=C
      
      write (*,*) A,B,C,D


      return 
      end


      subroutine IMFScalo

      real UM1,UM2
      real MI,M1,M3,MS
      real A,B
      real zita
      real Azita,Bzita,CZITA,DZITA
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
      MI=0.1                    !! Massa minima              
      M1=2.0                    !! Massa cambio pendenza     
      MS=100.0                  !! Massa massima                                          
      
      m3=1.                     !! SERVE PER IL CALCOLO DELLA NORMALIZZAZIONE "ZITA"

      ZITA=0.3
      Azita=MS**(1.+UM2)
      Bzita=M1**(1.+UM2)
      Czita=M1**(1.+UM1)
      Dzita=M3**(1.+UM1)        !mi senza zita o m3 se vuoi usare zita a 0.3!!!


      B=ZITA/(M1**(UM2-UM1)*(Czita-Dzita)/(1.+UM1)+
     $     (Azita-Bzita)/(1.+UM2))

      A=B*M1**(UM2-UM1)

      return
      end


      subroutine IMFSalpeter

      real UM1
      real MI, MS
      real A
      

      COMMON UM1,UM2,UM3,UM4
      COMMON A,B,C,D
      COMMON M1,M2,M3,M4
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     !   SALPETER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   
!!!   
      
      UM1=-1.35                 !! esponente salpeter                       
      MI=0.1                    !! Massa minima                  
      MS=80.0                   !! Massa massima                 
      
      
      Azita=MS**(1.+UM1)                                   
      Bzita=MI**(1.+UM1)                                   
      
      
      A=(1.+UM1)/(Azita-Bzita)                           
      write (*,*) A
 
      
      return 
      end



      subroutine MULTIKROUPA(Max,Min,Result)

      real Max,Min,Result
      real UM1,UM2,UM3,UM4
      real A,B,C,D
      real M1,M2,M3

      common UM1,UM2,UM3,UM4
      common A,B,C,D
      COMMON M1,M2,M3,M4


      if (MAX.le.M1) then 
         result=(MIN)**(UM1)-
     $        (MAX)**(UM1)
         result=result/(UM1)*A 
                                !! Frazione di stelle in numero in quell'intervallo di masse 
                                !  per unità di massa
                                !     
      else
         if (MAX.le.M2) then
                           
            result=(MIN)**(UM2)-
     $           (MAX)**(UM2)
            result=result/(UM2)*B 
                                !     ! Frazione di stelle in numero in quell'intervallo di masse 
                                !     per unità di massa
         else
            if (MAX.le.M3) then
               
               result=(MIN)**(UM3)-
     $              (MAX)**(UM3)
               result=result/(UM3)*C
            else
               result=(MIN)**(UM4)-
     $              (MAX)**(UM4)
               result=result/(UM4)*D
               
            endif
            
         endif
      endif
      if (MAX.gt.50) result=0.
      
      return
      end



      subroutine MULTISCALO(Max,Min,Result)

      real Max,Min,Result
      real UM1,UM2
      real A,B
      real M1

      common UM1,UM2,UM3,UM4
      common A,B,C,D
      COMMON M1,M2,M3,M4




      if (MAX.le.M1) then 
         
         result=(MIN)**(UM1)-
     $        (MAX)**(UM1)
         result=result/(UM1)*A 
!     !  Frazione di stelle in numero in quell'intervallo di masse 
!     !  per unità di massa
!     
      else
         result=(MIN)**(UM2)-
     $        (MAX)**(UM2)
         result=result/(UM2)*B 
!     ! Frazione di stelle in numero in quell'intervallo di masse 
!     ! per unità di massa


      endif

      return
      end
      
      subroutine MULTISALPETER(Max,Min,Result)

      real Max,Min,Result
      real UM1
      real A

      common UM1,UM2,UM3,UM4
      common A,B,C,D
      COMMON M1,M2,M3,M4

      result=(MIN)**(UM1)- 
     $     (MAX)**(UM1)     
      result=result/(UM1)*A        
      
      return 
      end
