       PROGRAM GDR
! Original program from Jean-Marc Simon
! To do to improve:
               ! - create first a regtest
               !-read dl_poly: box_size, position, should know the
               !config level
               !- proper name for output file
               !- allocate dynamically
               !- obtain info from command line
               !- run it from python which will create automatically the
               !command line
               !- use new fortran standard

      
         integer ntot, maxmol
         parameter (ntot=100000, maxmol=100000)
         character*2 atom
        INTEGER I,J,K,ij
        INTEGER  ndist
        integer NH,NO,NC, nmol, ncyc, ntotmax, natmol
        real*8 MH,MO,MC, mass(10), xr,yr,zr,Rsq2
        real*8 rx(maxmol), ry(maxmol), rz(maxmol)
        real*8 rxik, ryik, rzik, dist, dist2
        real*8 massmol, gdrstep
        real*8 histo(ntot),hist_id(ntot), nvol_id(ntot)
        real*8 dens
        REAL*8 RR(ntot)
       REAL*8 LBOX,VBOX,PI, LBOXS2, Lred
  
       character*200 GDR_File, traj_file
       
       PI=3.141592654
       write(*,*) "HEY"
       Open(30,FILE="info_gr", STATUS="UNKNOWN")
       read(30,*)
       read(30,*)
       read(30,*) traj_file
       read(30,*)
       read(30,*) ncyc
       read(30,*)
       read(30,*) lbox
       read(30,*) 
       read(30,*) nmol
       read(30,*)
       read(30,*) natmol
       read(30,*)
       do ij = 1,natmol
          read(30,*) mass(ij)
       enddo
       read(30,*) 
       read(30,*) gdrstep
       close(30)

       lboxs2=LBOX/2.0
       rsq2=lboxs2*sqrt(2.0)
       
       massmol=0.0
       
       ntotmax=int(rsq2/gdrstep)
       do ij=1,natmol
          write(*,*) ij, mass(ij)
          massmol=massmol+mass(ij)
       enddo
       write(*,*) "HOY"

       ! print*, massmol
       histo=0.0
       hist_id=0.0
       
       Open(20,FILE=traj_file, STATUS="UNKNOWN")
       do j=1,NCYC
          read(20,*)
         !         read(20,*)
         !write(*,*) "Start reading, cyc=", j
          if (mod(j,100).eq.0) write(*,*) "start cycle", j
          do i = 1,nmol
             !if (mod(i,1).eq.0) then 
             !   write(*,*) "start read molecule", i , "for cycle j =", j
             !end if
             !if (mod(i,100).eq.0) write(*,*) "start read molecule", i
             rx(i)=0.0
             ry(i)=0.0
             rz(i)=0.0
             do ij=1,natmol
                !write(*,*) ij
                READ(20,*) xr,yr,zr  !!! be careful to PBC for molecules
                !write(*,*) xr, yr, zr
                rx(i)= rx(i)+mass(ij)*xr
                ry(i)= ry(i)+mass(ij)*yr
                rz(i)= rz(i)+mass(ij)*zr
             enddo
             rx(i)=rx(i)/massmol
             ry(i)=ry(i)/massmol
             rz(i)=rz(i)/massmol
          enddo
          do i=1,nmol-1
             !if (mod(i,1).eq.0) then 
             !   write(*,*) "start calc molecule", i , "for cycle j =", j
             !end if
             rxi=rx(i)
             ryi=ry(i)
             rzi=rz(i)
             do k=i+1,nmol
                rxik=rxi-rx(k)
                rxik=rxik-LBOX*anint(rxik/LBOX)
                ryik=ryi-ry(k)
                ryik=ryik-LBOX*anint(ryik/LBOX)
                rzik=rzi-rz(k)
                rzik=rzik-LBOX*anint(rzik/LBOX)
                dist2=rxik**2+ryik**2+rzik**2
                dist=dsqrt(dist2)
                if (dist.LT.rsq2) then
                   Ndist=int(dist/gdrstep) +1
                   histo(ndist)=histo(ndist)+2.0 !! 2 dsitances une pour i et une pour k
                   if (dist.LT.2) then  !! test if yes problem of distances , PBC, ...
                      print*,dist, i,k,j
                      print*, rxi, ryi, rzi
                      print*, rx(k), ry(k), rz(k)
                      stop
                   endif
                endif
             enddo
          enddo
       enddo                    !! ncyc
       histo=histo/(ncyc*nmol)
      
       Vbox = LBOX*Lbox*Lbox
       dens=real(nmol)/vbox       !!! here N-1 or N !!!
       
       close(20)
       
       do i=1,ntotmax
          rr(i)=real(i)*gdrstep
          Nvol_id(i)=dens*(4.0*pi/3.0)*(rr(i)*rr(i)*rr(i))
          if (rr(i).GT.lboxs2)then
              Lred=LBOXS2/RR(I)
              Nvol_id(i)=Nvol_id(i)*(-2+9.0/2.0*Lred-3.0/2.0*Lred**3)
           endif
        enddo
        
       hist_id(1)=Nvol_id(1)
       do i=2,ntotmax          
          hist_id(i)=Nvol_id(i)-Nvol_id(i-1)
       enddo
       do i=1,ntotmax
          histo(i)=histo(i)/hist_id(i)
       enddo
       
       OPEN(10,FILE="GDR_CM", STATUS="UNKNOWN")
       do i=1,ntotmax
          write(10,*) rr(i)-gdrstep/2.0, histo(i), hist_id(i)
       enddo
       close(10)
       end

 !  100   FORMAT(1X, 5E14.6)

