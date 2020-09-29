       PROGRAM GDR
! Original program from Jean-Marc Simon
! To do to improve:
               !-read dl_poly: box_size, position, should know the
               !config level
               !- run it from python which will create automatically the
               !command line

       IMPLICIT NONE
     
       INTEGER, PARAMETER       :: dp = kind(1.0d0) 
       CHARACTER(len=2)         :: atom
       CHARACTER(len=200)       ::  GDR_File, traj_file, output_file
       INTEGER                  :: nconfig, trash, ntot
       INTEGER                  :: I,J,K,ij, natmol
       INTEGER                  :: ndist, NH,NO,NC, nmol, ncyc, ntotmax
       REAL(KIND=dp)            :: MH,MO,MC,  xr,yr,zr,Rsq2
       REAL(KIND=dp)            :: rxik, ryik, rzik, dist, dist2
       REAL(KIND=dp)            :: massmol, gdrstep
       REAL(KIND=dp)            :: dens, mass_at, rxi, ryi, rzi
       REAL(KIND=dp)            :: LBOX,VBOX,PI, LBOXS2, Lred
       REAL(KIND=dp),DIMENSION(:), ALLOCATABLE    ::  rx, ry, rz
       REAL(KIND=dp),DIMENSION(:), ALLOCATABLE    ::  histo, hist_id
       REAL(KIND=dp),DIMENSION(:), ALLOCATABLE    ::  mass, RR, nvol_id
       
       PI=3.141592654
       Open(30,FILE="info_gr", STATUS="UNKNOWN")
       read(30,*) traj_file
       read(30,*) output_file
       read(30,*) ncyc
       read(30,*) lbox
       read(30,*) nmol
       read(30,*) natmol
       ALLOCATE(mass(natmol))
       do ij = 1,natmol
          read(30,*) mass(ij)
       enddo
       read(30,*) gdrstep
       close(30)
 
       ntot=INT(LBOX/gdrstep)
       ALLOCATE(histo(ntot))
       ALLOCATE(hist_id(ntot))
       ALLOCATE(nvol_id(ntot))
       ALLOCATE(RR(ntot))
 
       ALLOCATE(rx(nmol))
       ALLOCATE(ry(nmol))
       ALLOCATE(rz(nmol))

       lboxs2=LBOX/2.0
       rsq2=lboxs2*sqrt(2.0)
       
       massmol=0.0
       
       ntotmax=int(rsq2/gdrstep)
       do ij=1,natmol
          massmol=massmol+mass(ij)
       enddo

       ! print*, massmol
       histo=0.0
       hist_id=0.0
       
       Open(20,FILE=traj_file, STATUS="UNKNOWN")
       read(20,*)
       read(20,*) nconfig, trash, trash
       do j=1,NCYC
          if (mod(j,100).eq.0) write(*,*) "start cycle", j
          read(20,*)  
          read(20,*) lbox, trash, trash
          read(20,*)  
          read(20,*)  
          do i = 1,nmol
             rx(i)=0.0
             ry(i)=0.0
             rz(i)=0.0
             do ij=1,natmol
                READ(20,*) trash, trash, mass_at, trash
                READ(20,*) xr,yr,zr  !!! be careful to PBC for molecules
                do k=1,nconfig
                   read(20,*)  
                enddo 
                rx(i)= rx(i)+mass_at*xr
                ry(i)= ry(i)+mass_at*yr
                rz(i)= rz(i)+mass_at*zr
             enddo
             rx(i)=rx(i)/massmol
             ry(i)=ry(i)/massmol
             rz(i)=rz(i)/massmol
          enddo
          do i=1,nmol-1
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
       
       OPEN(10,FILE=output_file, STATUS="UNKNOWN")
       do i=1,ntotmax
          write(10,*) rr(i)-gdrstep/2.0, histo(i), hist_id(i)
       enddo
       close(10)
       end

 !  100   FORMAT(1X, 5E14.6)

