       PROGRAM GDR
! Original program from Jean-Marc Simon
! To do to improve:
               !- run it from python which will create automatically the
               !command line

       IMPLICIT NONE
     
       INTEGER, PARAMETER       :: dp = kind(1.0d0) 
       CHARACTER(len=2)         :: atom
       CHARACTER(len=200)       :: traj_file, output_file
       CHARACTER(len=200)       :: label, method 
       INTEGER                  :: nconfig, trash, ntot, nlines
       INTEGER                  :: I,J,K,ij, natmol, nsteps, frequency
       INTEGER                  :: ndist, nmol, ncyc, ntotmax
       REAL(KIND=dp)            :: xr,yr,zr,Rsq2
       REAL(KIND=dp)            :: rxik, ryik, rzik, dist, dist2
       REAL(KIND=dp)            :: massmol, gdrstep
       REAL(KIND=dp)            :: dens, mass_at, rxi, ryi, rzi
       REAL(KIND=dp)            :: LBOX,VBOX,PI, LBOXS2, Lred
       REAL(KIND=dp)            :: lboy, lboz, sum_ 
       REAL(KIND=dp),DIMENSION(:), ALLOCATABLE    ::  rx, ry, rz
       REAL(KIND=dp),DIMENSION(:), ALLOCATABLE    ::  histo, hist_id
       REAL(KIND=dp),DIMENSION(:), ALLOCATABLE    ::  RR, nvol_id
       
       PI=3.141592654
       Open(30,FILE="info_gdr", STATUS="UNKNOWN")
       read(30,*) traj_file
       read(30,*) output_file
       read(30,*) nsteps
       read(30,*) frequency
       read(30,*) lbox
       read(30,*) nmol
       read(30,*) natmol
       read(30,*) method
       SELECT CASE(method)
       CASE('label')
               read(30,*) label
       END SELECT
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
       histo=0.0
       hist_id=0.0
       ncyc = 0
       
       Open(20,FILE=traj_file, STATUS="UNKNOWN")
       read(20,*)
       read(20,*) nconfig, trash, trash
       do j=1,nsteps
          sum_ = 0.0
          if (mod(j,100).eq.1) write(*,*) "start cycle", j
          if (mod(j,frequency).eq.1) then
          ncyc=ncyc+1
          read(20,*)  
          read(20,*) lbox, lboy, lboz
          read(20,*)  
          read(20,*)  
          do i = 1,nmol
             rx(i)=0.0
             ry(i)=0.0
             rz(i)=0.0
             massmol = 0.0
             do ij=1,natmol
                READ(20,*) atom, trash, mass_at, lboy
                READ(20,*) xr,yr,zr  !!! be careful to PBC for molecules
                do k=1,nconfig
                   read(20,*)  
                enddo
                SELECT CASE(method)
                CASE('label')
                    if (atom.EQ.label) then
                        rx(i) = xr
                        ry(i) = yr
                        rz(i) = zr
                    endif
                    massmol=1
                CASE('COM')
                    rx(i)= rx(i)+mass_at*xr
                    ry(i)= ry(i)+mass_at*yr
                    rz(i)= rz(i)+mass_at*zr
                    massmol = massmol + mass_at
                END SELECT
                if (atom.EQ.'C') then
                        rxi=xr
                        ryi=yr
                        rzi=zr
                endif
             enddo
             rx(i)=rx(i)/massmol
             ry(i)=ry(i)/massmol
             rz(i)=rz(i)/massmol
             sum_ = sum_ +  dsqrt((rx(i)-rxi)**2+(rx(i)-rxi)**2+(rx(i)-rxi)**2)
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
                      print*, rxik, ryik, rzik
                      print*, rx(k), ry(k), rz(k)
                      print*, rx(i), ry(i), rz(i)
                      stop
                   endif
                endif
             enddo
          enddo
          else
                  nlines=4+(2+nconfig)*nmol*natmol
                  do i=1,nlines
                     read(20,*) 
                  enddo
          end if 
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

