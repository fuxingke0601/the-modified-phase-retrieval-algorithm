!=====================================================================================
! (c) IHEP 2023. All rights reserved.
! Author: Geng Zhi
! Institute of High Energy Physics, Chinese Academy of Science (IHEP, CAS). 
! If you have any problem with the program, please contact author with the
! following 
! email: gengz@ihep.ac.cn
!======================================================================================
!未加入sayre等式

module functions
   use,intrinsic :: iso_c_binding
   use fgsl
   implicit none
   include 'fftw3.f03'

   !global variables
   real,parameter :: pi = 3.1415926
   integer nu, nv, nw      !grid of density map
   integer nrefls          !number of reflections in asymmetric unit  (including 000)

   !structure body for reflections
   !NOTE:opha and owt are used for emap, while npha and nwt are used for comparing
   type reflection
      integer h, k, l  !miller index of this reflection
      real amp, pha    !observed amplitude and phase
      real enorm, gnorm!true normalized(E) amplitude and updated E(h)~G(h)=sum(E(k)*E(h-k)) (sayre equation)
      real opha, owt   !old phase and weight of E(h) : used for emap
      real npha, nwt   !new phase and weight of E(h) : used for comparing
      real sigf        !sigma value of observed amplitude
      logical missing  !true means missing reflection
      logical weak     !true means weak reflections shifted by pi/2
   end type reflection

   type(reflection),allocatable :: hkls(:)  !memory for saving reflections information

   !asymmetry unit atom content
   type atom
      character(4) atom_name
      real num
   end type atom

   type(atom),allocatable :: atoms(:)

   !space group information
   character(11) spgrp  !name of space group
   integer nspgrp       !number code of space group 
   character(11) nampg  !name of point group
   character(6) launam  !name of Laue group
   character ltype      !lattice type
   integer nlaue        !number of Laue group
   integer nsymp        !primitive number of space group(remove centering)
   integer nsym         !total number of space group
   real rsym(4,4,192)!matrix operators of space group

   !unitcell parameter
   real cell(6)         !unit cell
   real rcell(6)        !reciprocal unit cell
   real vol , rvol      !volume of unitcell and reciprocal
   real rrr(3,3,6)      !othorgonal matrix:normally use rrr(:,:,1) as a usual convention

   logical mtz_has_pha    !input mtz file contains initial phase information
   logical use_ecal       !use normalized structure factor during charge flipping
   logical use_wilson     !correct for temperature effect before charge flipping
   integer(8) wilson_bins !number of bins for wilson plot
   real  wilson_res       !the lowest resolution truncated for wilson plot
   real frac              !fraction of density zone
   integer ntry           !try number of different initial random phases
   integer niteration     !number of charge-flipping iterations
   integer ecal_ref       !number of reflections of each bin during Ecalc

   contains

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !    pdbin  :     read in pdb files
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   subroutine pdbin(x,y,z,atms)
      implicit none
      integer ifail, ixyzin, n
      integer nasym, natms
      integer iser,irs,iz
      character(4) atnam, resnam, resno, segid, id
      character chnnam, inscod, altcod
      real occ,biso,U(6)
      real,allocatable :: x(:), y(:), z(:)
      character(4),allocatable :: atms(:)

      ifail = 0 !stop when confront an error
      ixyzin = 0
      !open xyzfile
      call xyzopen('XYZIN','INPUT',' ',ixyzin,ifail)
      if(ifail .eq. -1) stop 'Open xyzfile Failed!'

      !count total atoms
      nasym = 0
10    call xyzadvance(ixyzin,0,0,*10,*100)
      !following breaks when confronting TER
!10    call xyzadvance(ixyzin,0,1,*100,*100)
      call xyzatom(ixyzin,iser,atnam,resnam,chnnam,irs,resno,inscod,&
         altcod,segid,iz,id)
      nasym = nasym + 1

      goto 10

100   continue

      !read in unitcell and vol
      call rbcell(cell,vol)
      !get orthogonal matrix
      call rbfro1(cell,vol,rrr)
      !read in unitcell and vol in reciprocal
      call rbrcel(rcell,rvol)
      !read in space group
      call rbspgrp(spgrp)
      !get space group matrix operators and some other information
      call msymlb(22,nspgrp,spgrp,nampg,nsymp,nsym,rsym)
      !get laue group information
      call pgnlau(nampg,nlaue,launam)

      ltype = spgrp(1:1)

      write(6,'(//,a,//)')'******** Information of space group and pdb *******'
      write(6,'(a,6f8.2,/,a,a, I4,2a)') 'Unitcell: ', cell, 'space group: ', spgrp,nspgrp, ' Lattice type: ', ltype
      write(6,'(2a,/,a,I4)') 'Point Group: ', nampg, 'Total numbers of SP: ', nsym
      write(6,'(a,I4)') 'Number of SP(exclude centering)', nsymp
      write(6,'(a,I4)') 'Total numbers of atoms in asymetric unit: ', nasym
      write(6,'(2a,i4)') 'Laue Group: ', launam,nlaue
      write(6,'(//,a,//)')'***************************************************'

      natms = nasym * nsym
      allocate(x(natms))
      allocate(y(natms))
      allocate(z(natms))
      allocate(atms(natms))

      !read in atoms
      call xyzrewd(ixyzin)
      n = 0
20    call xyzadvance(ixyzin,0,0,*20,*200)
      n = n + 1
      call xyzatom(ixyzin,iser,atnam,resnam,chnnam,irs,resno,inscod,&
         altcod,segid,iz,atms(n))
      !read in fractional coordinates('F') or cartesian('O')
      call xyzcoord(ixyzin,'F','U',x(n),y(n),z(n),occ,biso,u)
      goto 20

200   continue

      call xyzclose(ixyzin)

      return
   end subroutine pdbin


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! expand_atoms: read in atomic scattering factor and expand asymmetric atoms
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   subroutine expand_atoms(x,y,z,atms,asf)
      implicit none
      real a(4), b(4), c !C-M coefficient
      real iwt           !atomic weight
      real ielec         !electron number
      real cu(2), mo(2)  !f'and f'' at Cu or Mo target wavelength
      integer ifail, idum, lhist, eof, i, j, nasym
      character(20) line
      real,allocatable :: x(:), y(:), z(:), asf(:,:)
      character(4),allocatable :: atms(:)

!      call ccpdpn(lhist,'atomsf.lib','READONLY','F',idum,ifail)

      allocate(asf(9,size(x)))

      nasym = size(x)/nsym
      outer:do i = 1, nasym
         call sfread2(adjustl(atms(i)),5,a,b,c,iwt,ielec,cu,mo,ifail)
         asf(:,i) = (/a, b, c/)
      end do outer

      !Now expand atoms into full space using space group and fractional
      !coordinates
      do i = 1, nasym
         do j = 2, nsym
            x((j-1)*nasym+i) = rsym(1,1,j)*x(i)+rsym(1,2,j)*y(i)+rsym(1,3,j)*z(i)+rsym(1,4,j)
            y((j-1)*nasym+i) = rsym(2,1,j)*x(i)+rsym(2,2,j)*y(i)+rsym(2,3,j)*z(i)+rsym(2,4,j)
            z((j-1)*nasym+i) = rsym(3,1,j)*x(i)+rsym(3,2,j)*y(i)+rsym(3,3,j)*z(i)+rsym(3,4,j)
            atms((j-1)*nasym+i) = atms(i)
            asf(:,(j-1)*nasym+i) = asf(:,i)
         end do
      end do

!      close(lhist)

      return
   end subroutine expand_atoms


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !        sfall    :    calculate structure factor in asymmetric unit
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   subroutine sfall(x,y,z,asf)
      implicit none
      real,allocatable :: x(:), y(:), z(:), asf(:,:)
      real :: res = 1.0     !resolution of to be calculated reflections
      integer h, k, l, hmax, kmax, lmax, hmin, kmin,lmin
      integer ihkl(3), jhkl(3), isym
      integer nhkl, n
      logical,allocatable :: asymhkl(:,:,:)
      real q, f             !atmoic scattering factor
      real tmp, sum1, sum2
      complex sf            !structure factor
      integer isysab

      !BASIC VARIABLES FOR PARSER
      integer,parameter :: maxpar = 200 !permited max splited segments of line
      character(4) key            !header label information in run.exe 
      character(400) line         !line for each parser to read each command
      integer ibeg(maxpar),iend(maxpar) !split each line divided by space and store 
                                  !each segment(namely,first locality and last
                                  !locality in the line)
      integer ityp(maxpar), idec(maxpar)!ityp denotes real(2) or character(1)
      integer ntok                !number of individual segments
      character(4) cvalue(maxpar) !characters of individual segment
      real fvalue(maxpar)         !real number of individual segment(begin from second)
      logical lend


      real adata(5)         !store h k l F and Phi for writing mtz 
      integer mtzout        !mtz file identifier
      character(30) title   !title of mtz file
      integer nlprg0        !total colums of mtz
      character(5) lsprg0(200) !column labels representing each column of mtz
      character ctprg0(200) !single character column type
      real :: bfactor =  25  !temperator factor
      real :: egauss = -10   !Gaussian noise egauss*sqrt(SF)

      real(8), external :: ZBQLNOR

      !CALL PARSER TO READ IN LABIN INFORMATION
10    continue
      line = ' '
      key = ' '
      ntok = maxpar
      call parser(key,line,ibeg,iend,ityp,fvalue,cvalue,idec,ntok,lend,.FALSE.)
      if(lend) goto 50   !end of reading,that is,no more lines
      if(key .eq. 'TITL')then !title of mtz
         title = line
         goto 10
      else if(key .eq. 'RESO')then
         res = fvalue(2)
         goto 10
      else if(key .eq. 'BFAC')then
         bfactor = fvalue(2)
         goto 10
      else if(key .eq. 'EGAU')then
         egauss = fvalue(2)
         goto 10
      else if(key .eq. 'END')then
         goto 50
      else
         write(6,'(a,a)') 'ERROR:not recognized label:',key
         stop
      end if

50    continue

      nhkl = 0
      hmax=nint(cell(1)/res); kmax=nint(cell(2)/res); lmax=nint(cell(3)/res)
      !dependent on specific spacegroup
      hmin =-hmax ; kmin =-kmax; lmin =-lmax
      write(6,'(3(a,I6))') 'Hmax: ', hmax, ' Kmax: ', kmax, ' Lmax: ', lmax
      write(6,'(3(a,I6))') 'Hmin: ', hmin, ' Kmin: ', kmin, ' Lmin: ', lmin

      allocate(asymhkl(hmin:hmax,kmin:kmax,lmin:lmax))
      asymhkl = .false.

      call epsln(nsym,nsymp,rsym,0)
      !choose asymmetric unit
      call asuset(spgrp,nspgrp,nampg,nsym,rsym,nsymp,nlaue,.true.)
      do l = lmin, lmax
      do k = kmin, kmax
      do h = hmin, hmax
         call rbrecip(h,k,l,q)
         if(sqrt(q) .gt. 1/res) cycle
         if(h.eq.0 .and. k.eq.0 .and. l.eq.0) cycle
         ihkl = (/h,k,l/)
         call asuput(ihkl,jhkl,isym)
         !exclude system absence
         call sysab(jhkl,isysab)
         if(isysab.eq.1) cycle
         asymhkl(jhkl(1),jhkl(2),jhkl(3)) = .true.
      end do
      end do
      end do

      nhkl = 0
      do l = lmin, lmax
      do k = kmin, kmax
      do h = hmin, hmax
         if(asymhkl(h,k,l)) nhkl = nhkl + 1
      end do
      end do
      end do

      write(6,'(a,i6)') 'Asymmetric unit contains indexs number: ', nhkl

      !write header into mtz
      mtzout = 1
      title = 'structures factors from pdb'
      nlprg0 = 5
      data lsprg0 /'H','K','L','F','PHI',195*' '/
      data ctprg0 /'H','H','H','F','P',195*' '/
      !open mtz file for writing
      call lwopen(mtzout,'HKLOUT')
      !write title
      call lwtitl(mtzout,title,0)
      !write cell information
      call lwcell(mtzout,cell)
      !store column labels in the mtz header
      call lwclab(mtzout,lsprg0,nlprg0,ctprg0,0)
      !write symmetry information
      call lwsymm(mtzout,nsym,nsymp,rsym,ltype,nspgrp,spgrp,nampg)
      !store sort order to be used
      call lwsort(mtzout,(/1,2,3,0,0/))

      sum1=0; sum2 = 0
      !calculate sf in asymmetric unit
      do l = lmin, lmax
      do k = kmin, kmax
      do h = hmin, hmax
         if( .not. asymhkl(h,k,l)) cycle
         adata(1:3)=real((/h,k,l/))
         call rbrecip(h,k,l,q)
         q = q /4

         sf = cmplx(0,0)

         do n = 1, size(x)
            f = sum(asf(1:4,n)*exp(-asf(5:8,n)*q)) + asf(9,n)
            sf = sf + f*cmplx(cos(2*pi*(h*x(n)+k*y(n)+l*z(n))),sin(2*pi*(h*x(n)+k*y(n)+l*z(n))))
         end do

         !add debye temperator factor
         sf = sf * exp(-bfactor*q)

         tmp = cabs(sf)
         adata(5) = atan2(aimag(sf),real(sf))/pi*180
         if(adata(5) .lt. 0) adata(5) = adata(5) + 360.

         if(egauss .gt. 0)then
            adata(4) = ZBQLNOR(dble(tmp),dble(egauss*sqrt(tmp)))
            sum1 = sum1+abs(tmp-adata(4))
            sum2 = sum2+tmp
         else
            adata(4) = tmp
         end if

         !write each reflection into mtz
         call lwrefl(mtzout,adata)

      end do
      end do
      end do

      call lwclos(mtzout,1)

      if(egauss .gt. 0)then
         write(6,'(/,a,/)')'************************************************************************'
         write(6,'(a,f6.2)') 'Added gaussian error can be quantified by R-factor : ', sum1/sum2
         write(6,'(/,a,/)')'************************************************************************'
      end if

      return
   end subroutine sfall


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !read in reflections and symmetry information
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   subroutine mtzin
      implicit none
      !BASIC VARIABLES FOR READING MTZ
      integer,parameter :: mcols = 200  !permited max columns of mtz
      integer,parameter :: maxpar = 200 !permited max splited segments of line
      integer mtzind        !mtz file identifier
      character(30) title   !title of mtz file
      integer nlprg0        !total colums of mtz
      character(5) lsprg0(mcols) !column labels representing each column of mtz
      character ctprg0(mcols) !single character column type
      integer lookup(mcols)
      real maxres, minres   !max and min resolution of mtz(virtually 1/d^2)
      real weak, cutoff     !fraction of weak reflections
      real(8),allocatable :: intensity(:)

      !BASIC VARIABLES FOR PARSER
      character(4) key            !header label information in run.exe 
      character(400) line         !line for each parser to read each command
      integer ibeg(maxpar),iend(maxpar) !split each line divided by space and store 
                                  !each segment(namely,first locality and last
                                  !locality in the line)
      integer ityp(maxpar), idec(maxpar)!ityp denotes real(2) or character(1)
      integer ntok                !number of individual segments
      character(4) cvalue(maxpar) !characters of individual segment
      real fvalue(maxpar)         !real number of individual segment(begin from second)
      real recin(mcols), s
      logical lend, mtzeof, logmss(mcols)
      integer ifail, i

      character(10) versnx
      integer numcol, nreflx, hmax, kmax, lmax
      real rngmtz(2,mcols)

      mtzind = 1
      use_ecal = .false.
      use_wilson = .false.

      !OPEN MTZ AND READING SPACE GROUP INFORMATION
      call lropen(mtzind,'HKLIN',0,ifail)
      call lrinfo(mtzind,versnx,numcol,nreflx,rngmtz)

      !get space group information from mtz
      call lrcell(mtzind,cell)
      !get maximum resolution
      call lrrsol(mtzind,minres,maxres)
      !get space group and other symmetry information
      call lrsymi(mtzind,nsymp,ltype,nspgrp,spgrp,nampg)
      !get symmetry operaters
      call lrsymm(mtzind,nsym,rsym)
      !generate laue group
      call pgnlau(nampg,nlaue,launam)
      !generate volume and orthogonal matrix
      call rbfro1(cell,vol,rrr)

      allocate(hkls(nreflx+1))
      allocate(intensity(nreflx+1))

      !CALL PARSER TO READ IN LABIN INFORMATION
10    continue
      line = ' '
      key = ' '
      ntok = maxpar
      call parser(key,line,ibeg,iend,ityp,fvalue,cvalue,idec,ntok,lend,.FALSE.)
      if(lend) goto 50   !end of reading,that is,no more lines
      if(key .eq. 'TITL')then !title of mtz
         title = line
         goto 10
      else if(key .eq. 'LABI')then

         mtz_has_pha = .false.
         do i = 2, ntok
            if (trim(adjustl(cvalue(i))).eq.'PHI')then
                mtz_has_pha = .true.
                exit
            end if
         end do

         data lsprg0 /'H','K','L','F',196*' '/
         data ctprg0 /'H','H','H','E',196*' '/
         data lookup /4*-1,196*0/

         if(mtz_has_pha)then
            nlprg0 = 6
            lsprg0(5) = 'SIGF'
            ctprg0(5) = 'Q'
            lookup(5) = -1
            lsprg0(6) = 'PHI'
            ctprg0(6) = 'P'
            lookup(6) = -1
         else
            nlprg0 = 5
            lsprg0(5) = 'SIGF'
            ctprg0(5) = 'Q'
            lookup(5) = -1
         end if

         call lkyin(mtzind,lsprg0,nlprg0,ntok,line,ibeg,iend)

         goto 10

      else if(key .eq. 'RESO')then
         maxres = 1./(fvalue(2))**2
         goto 10
      else if(key .eq. 'WILS')then
         use_wilson = .true.
         goto 10
      else if(key .eq. 'NWLS')then
         wilson_bins = nint(fvalue(2))
         goto 10
      else if(key .eq. 'RWLS')then
         wilson_res = fvalue(2)
         goto 10
      else if(key .eq. 'FRAC')then
         frac = fvalue(2)
         goto 10
      else if(key .eq. 'WEAK')then
         weak = fvalue(2)
         goto 10
      else if(key .eq. 'NTRY')then
         ntry = nint(fvalue(2))
         goto 10
      else if(key .eq. 'NITR')then
         niteration = nint(fvalue(2))
         goto 10
      else if(key .eq. 'CONS')then

         allocate(atoms((ntok-1)/2))
         do i = 1, (ntok-1)/2
            atoms(i)%atom_name = cvalue(i*2)
            atoms(i)%num = fvalue(i*2+1)
         end do
         goto 10

      else if(key .eq. 'ECAL')then
         use_ecal = .true.
         goto 10
      else if(key .eq. 'NECA')then
         ecal_ref = nint(fvalue(2))
         goto 10
      else if(key .eq. 'END')then
         goto 50
      else
         write(6,'(a,a)') 'ERROR:not recognized label:',key
         stop
      end if

50    continue

      !set up column assignments(namely,ctprg0) based on LABIN
      call lrassn(mtzind,lsprg0,nlprg0,lookup,ctprg0)

      !READ IN REFLECTIONS FROM MTZ
      nrefls = 0
      hmax = 0; kmax = 0; lmax = 0

110   continue
      call lrrefl(mtzind,s,recin,mtzeof)
      if(mtzeof) goto 230
      !check whether beyond resolution limited
      if(s .gt. maxres) goto 110

      nrefls = nrefls + 1
      !check for missing number flags
      call lrrefm(mtzind,logmss)

      if(logmss(lookup(4)))then
         !missing reflections
         hkls(nrefls)%h = nint(recin(1))
         hkls(nrefls)%k = nint(recin(2))
         hkls(nrefls)%l = nint(recin(3))
         hkls(nrefls)%amp = 0    !initial amplitude of missing spots are set as 0
         hkls(nrefls)%sigf = 0
         hkls(nrefls)%missing = .true.
      else
         !observed reflections
         hkls(nrefls)%h = nint(recin(1))
         hkls(nrefls)%k = nint(recin(2))
         hkls(nrefls)%l = nint(recin(3))
         hkls(nrefls)%amp = recin(lookup(4))
         hkls(nrefls)%sigf = recin(lookup(5))
         hkls(nrefls)%missing = .false.
         if(mtz_has_pha)then
            hkls(nrefls)%opha = recin(lookup(6))/180*pi
         end if
      end if

      hmax = max(hmax,hkls(nrefls)%h)
      kmax = max(kmax,hkls(nrefls)%k)
      lmax = max(lmax,hkls(nrefls)%l)

      goto 110

      !end of mtz
230   continue

      !consider adding systematic absent reflections(observed and must be zero)

      !including (0 0 0) index
      nrefls = nrefls + 1
      hkls(nrefls)%h = 0; hkls(nrefls)%k = 0; hkls(nrefls)%l = 0
      hkls(nrefls)%amp = 0; hkls(nrefls)%missing = .true.

      write(6,'(a,I7)') 'Reading MTZ completed,total reflections:',nrefls
      call lrclos(mtzind)

!      nu = 4*hmax; nv = 4*kmax; nw = 4*lmax
!      nu = 2*nint(cell(1)*sqrt(maxres)); nv = 2*nint(cell(2)*sqrt(maxres)); nw = 2*nint(cell(3)*sqrt(maxres))
!      nu = 2*cell(1)*sqrt(maxres); nv = 2*cell(2)*sqrt(maxres); nw = 2*cell(3)*sqrt(maxres)
      nu = 2*hmax+1; nv = 2*kmax+1; nw = 2*lmax+1
      !nu=55; nv=33; nw=153
      !sort all intensities including missing and label the weak reflections
      do i = 1, nrefls
         intensity(i) = hkls(i)%amp
      end do
      call quicksort(intensity,1,nrefls,1,1,nrefls)
      cutoff = intensity(nint(weak*nrefls)+1)
      do i = 1, nrefls
         if(hkls(i)%amp .lt. cutoff)then
            hkls(i)%weak = .True.
         else
            hkls(i)%weak = .False.
         end if
      end do

      write(6,'(/,a,/)') '*************************************************************'
      write(6,'(/,3(a,I4),/)') 'HMAX : ', hmax, ' KMAX : ', kmax, ' LMAX : ', lmax
      write(6,'(/,a,3I5,a,/)') 'Grid of density map : (',nu,nv,nw ,')'
      write(6,'(/,a,F8.2/)') 'Amplitude cutoff considered to be weak : ', cutoff
      write(6,'(/,a,/)') '*************************************************************'

      return
   end subroutine mtzin

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !  wilson      :     calculate debye-waller temperator factor
   !                    nbins : number of sampling grids sin(theta)/lambda
   !                    reslow  : the lowest resolution used for statistics
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   subroutine wilson(nbins,reslow)
      implicit none
      integer(8) nbins
      integer n
      integer i, j, ifail, status
      real q, v, qmax, qmin, qlow, qinc, ed(nrefls)
      real reslow, ress(nbins), wt(nbins)
      real a(4), b(4), c
      real iwt, ielec, cu(2), mo(2)
      integer mul, eps
      real(8) I_abs(142), t(142), tt, f, xx
      real(8) btemp, scal, cov00, cov01, cov11, chisq
      real hist(nbins), q2(nbins)      !In(<Iobs>/<Iabs>) and (sin(theta)/lambda)**2
      integer num(nbins)               !number of reflections in one bin

      type(fgsl_interp_accel) :: acc
      type(fgsl_spline) spline

      I_abs = 0
      !first calculate <Iabs> from unitcell content and make a sample table
      do i = 1, 142
	 t(i) = 0.01*i     !sin(theta)/lambda (increment is 0.01)
	 tt = t(i) * t(i)  !(sin(theta)/lambda)**2
	 do j = 1, size(atoms)
	    call sfread2(adjustl(atoms(j)%atom_name),5,a,b,c,iwt,ielec,cu,mo,ifail)
	    f = sum(a*exp(-b*tt))+c
	    I_abs(i) = I_abs(i) + f**2*atoms(j)%num*nsym
	 end do
      end do

      !next calculate average normalized structure factor : Iobs/(<Iabs>*epsilon)
      acc = fgsl_interp_accel_alloc()
      spline = fgsl_spline_alloc(fgsl_interp_cspline,142_8)
      status = fgsl_spline_init(spline,t,I_abs,142_8)
      qmax = 0; qmin = 100

      do i = 1, nrefls
         if(hkls(i)%missing) cycle !not include missing reflections
         call resolution(hkls(i)%h,hkls(i)%k,hkls(i)%l,q,v)
	 call muleps(hkls(i)%h,hkls(i)%k,hkls(i)%l,mul,eps)
	 if(i .eq. nrefls)then !for (0 0 0), not interpolation
	    ed(i) = (hkls(i)%amp)**2/I_abs(1)/eps
	    exit
	 end if
	 !interpolate <Iabs> from sample table
	 q = q / 2        !sin(theta)/lambda
	 if(i .ne. nrefls)then
	    qmax = max(qmax,q)
	    qmin = min(qmin,q)
         end if
         if(q .le. t(1))then
            xx = I_abs(1)
         else if(q .ge. t(142))then
            xx = I_abs(142)
         else
	    xx = fgsl_spline_eval(spline,dble(q),acc)
         end if
	 ed(i) = (hkls(i)%amp)**2/xx/eps
      end do

      call fgsl_spline_free(spline)
      call fgsl_interp_accel_free(acc)


      !now make histogram of log(<Iobs>/(<Iabs>*epsilon))
      qlow = 1./reslow/2
      qmin = max(qlow,qmin)
      qinc = (qmax**2-qmin**2) / nbins    !width of one bin (based on (sin(theta)/lambda)**2)
      hist = 0; q2 = 0; num = 0
      !Iobs/(<Ical>*epsilon) is binned as follows
      do i = 1, nrefls
         if(hkls(i)%missing) cycle !not include missing reflections
         call resolution(hkls(i)%h,hkls(i)%k,hkls(i)%l,q,v)
	 call muleps(hkls(i)%h,hkls(i)%k,hkls(i)%l,mul,eps)
	 if(q/2 .lt. qmin) cycle
	 n = nint(((q/2)**2-qmin**2)/qinc)
	 if(n .lt. 1) n = 1
	 if(n .gt. nbins) n = nbins
	 num(n) = num(n) + mul
	 q2(n) = q2(n) + (q/2)**2*mul
	 hist(n) = hist(n) + ed(i)*mul
      end do
      !and smooth histograms by adding two adjacent hists
      do i = 1, nbins-1
         hist(i) = hist(i) + hist(i+1)
	 num(i) = num(i) + num(i+1)
	 q2(i) = q2(i) + q2(i+1)
	 num(i) = max(num(i),1)
      end do
      !average hist
      do i = 1, nbins
	 q2(i) = q2(i)/num(i)             !average (sin(theta)/lambda)**2
	 ress(i) = 1./sqrt(q2(i))/2       !average resolution
	 wt(i) = num(i)*1.0/sum(num)      !weight of one bin
	 hist(i) = alog(hist(i)/num(i))   !average log(<Iobs>/(<Iabs>*epsilon))
      end do


      !next weighted linear fitting to get scale factor and temperator factor
      status = fgsl_fit_wlinear(dble(q2),1_8,dble(wt),1_8,dble(hist),1_8,nbins,scal,btemp,cov00,cov01,cov11,chisq)
      scal = exp(scal); btemp = -btemp/2

      write(6,'(/,a,/)') '**************       Wilson Statistics        ********************'
      write(6,'(a)') 'Least Square Weighted Fitting Of Temperator And Scale Factor'
      write(6,'(a,/)') 'log(<Iobs>/<Iabs>) = log(scale)-2*B*(sin(theta)/lambda)**2'
      write(6,'(a,/)')   'Fitting Result: '
      write(6,'(a,F6.4,a,F6.3,/)') 'scale(k): ', scal, '   temperator factor(B): ', btemp
      write(6,'(a,F6.4)') 'Residual: ', chisq
      write(6,'(a,/)') '******************************************************************'

      !now leave out temperator factor effect
      do i = 1, nrefls
         if(hkls(i)%missing) cycle !not include missing reflections
         call resolution(hkls(i)%h,hkls(i)%k,hkls(i)%l,q,v)
	 q = (q/2)**2
	 hkls(i)%amp = hkls(i)%amp * exp(btemp*q)
      end do

      return
   end subroutine wilson


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !   ecal      : normalize structure factor, namely :  E=(F/sqrt(epsilon))/(<F/sqrt(epsilon)>)
   !               wilson_bins and wilson_res are used for wilson
   !               ecal_ref is the average numbers in one bin
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   subroutine ecal
      implicit none
      integer(8) ecal_bins
      real bin_width
      real qmax, qmin, q, vol, emax, ee
      integer i, n, mul, eps
      real(8),allocatable :: ed(:), sed(:), t(:), ed1(:)
      real(8) xx, tzc, tza, oza, ozc
      integer,allocatable :: num(:)
      integer,allocatable :: cump(:),cumpa(:),cumpc(:)
      real sec(9), sea(9)           !sum of all centric and acentric normalized structures
      integer nsec, nsea, ncump, n_observed
      real :: tea(9) = (/0.886,1.0,1.329,2.0,3.323, 6.0,0.736,1.0,2.0/)  !theoretical values of acentric
      real :: tec(9) = (/0.798,1.0,1.596,3.0,6.383,15.0,0.968,2.0,8.0/)  !theoretical values of centric
      integer status

      type(fgsl_interp_accel) :: acc
      type(fgsl_spline) spline

      !first call wilson to divide out temperator factor
      call wilson(wilson_bins,wilson_res)

      !then set up elementary bins used for average <F**2/epsilon>
      qmax = 0; qmin = 1000; n_observed = 0
      do i = 1, nrefls
         if(hkls(i)%missing) cycle !not include missing reflections
         call resolution(hkls(i)%h,hkls(i)%k,hkls(i)%l,q,vol)
	 qmax = max(qmax,q)
	 qmin = min(qmin,q)
         n_observed = n_observed + 1
      end do
      !NOTE : bins are divided based on q**3=(2*sin(theta)/lambda)**3
      qmax = qmax**3; qmin = qmin**3
      ecal_bins = nint((n_observed+ecal_ref/2)*1.0/ecal_ref) !number of bins
      bin_width = (qmax-qmin)/ecal_bins                    !width of one bin

      !now make a histogram of <F**2/epsilon>
      allocate(ed(ecal_bins),sed(ecal_bins),num(ecal_bins),t(ecal_bins),ed1(ecal_bins))
      ed = 0; num = 0
      do i = 1, nrefls
         if(hkls(i)%missing) cycle !not include missing reflections
         call resolution(hkls(i)%h,hkls(i)%k,hkls(i)%l,q,vol)
	 call muleps(hkls(i)%h,hkls(i)%k,hkls(i)%l,mul,eps)
	 q = q**3
	 n = nint((q-qmin)/bin_width+1)
	 if(n .gt. ecal_bins) n = ecal_bins
	 ed(n) = ed(n) + mul*(hkls(i)%amp)**2/eps
	 num(n) = num(n) + mul
      end do

      !and smooth histograms by adding neigbering two bins
      if(num(1)+num(2).eq.0)then     !for the 1th bin
	 sed(1) = 0
      else
	 sed(1) = (0.45*ed(1)+0.55*ed(2))/(0.45*num(1)+0.55*num(2))
      end if
      do i = 2, ecal_bins-1          !for the middle bin
	 if(num(i-1)+num(i)+num(i+1) .eq. 0)then
	    sed(i) = 0
	 else
	    sed(i) = (0.4*ed(i)+0.3*(ed(i-1)+ed(i+1)))/(0.4*num(i)+0.3*(num(i-1)+num(i+1)))
	 end if
      end do
      if(num(ecal_bins-1)+num(ecal_bins) .eq. 0)then  !for the last bin
	 sed(ecal_bins) = 0
      else
	 sed(ecal_bins) = (0.55*ed(ecal_bins-1)+0.45*ed(ecal_bins))/(0.55*num(ecal_bins-1)+0.45*num(ecal_bins))
      end if

      !normalize structure factor
      do i = 1, ecal_bins             !calculate q**3 respectively corresponding to each bin
         t(i) = qmin + (i-1)*bin_width + bin_width/2
      end do
      acc = fgsl_interp_accel_alloc()
      spline = fgsl_spline_alloc(fgsl_interp_cspline,ecal_bins)
      status = fgsl_spline_init(spline,t,sed,ecal_bins)
      do i = 1, nrefls
         if(hkls(i)%missing) cycle !not include missing reflections
         call resolution(hkls(i)%h,hkls(i)%k,hkls(i)%l,q,vol)
	 call muleps(hkls(i)%h,hkls(i)%k,hkls(i)%l,mul,eps)
	 q = q**3
	 if(q .ge. t(ecal_bins))then
	    xx = sed(ecal_bins)
	 else if(q .le. t(1))then
	    xx = sed(1)
	 else
    	    xx = fgsl_spline_eval(spline,dble(q),acc)
	 end if
	 hkls(i)%enorm = sqrt((hkls(i)%amp)**2/xx/eps)
      end do
      call fgsl_spline_free(spline)
      call fgsl_interp_accel_free(acc)

      !summary information of normalized structure factor as follows :
      !   1. Z=E**2: distribution of  <|Z-1|> vs shell
      !   2. acentric and centric of total average: <E**(1-6)> and <|E**2-1|^(1-3)>
      !   3. N(Z)=erf(sqrt(Z/2)) for centric   and  N(Z)= 1-exp(-Z) for acentric
      !   4. parity group for <E**2>, <(E**2-1)**2> 

      !first make a histogram of <|Z-1|>
      ed1 =0
      do i = 1, nrefls
         if(hkls(i)%missing) cycle !not include missing reflections
         call resolution(hkls(i)%h,hkls(i)%k,hkls(i)%l,q,vol)
	 call muleps(hkls(i)%h,hkls(i)%k,hkls(i)%l,mul,eps)
	 q = q**3
	 n = min(nint((q-qmin)/bin_width+1),ecal_bins)
	 ed1(n) = ed1(n) + abs((hkls(i)%enorm**2)-1)*mul    !sum of |Z-1| of one bin
      end do
      ed = ed /num                                       !average of <F^2> of one bin
      ed1 = ed1 / num                                    !average of <|Z-1|> of one bin
      !now output histogram information
      write(6,'(/,a,/)') '********************   Histogram Of Binned data   *********************'
      write(6,'(a)') '   shell    resolution   nrefls   <F^2>   smoothed<F^2>    <|E^2-1|>   '
      do i = 1, ecal_bins
         write(6,'(3X,I3,6X,F6.3,7X,I4,3X,F8.2,3X,F8.2,6X,F8.4)') i, (t(i))**(-1./3.), num(i), ed(i), sed(i), ed1(i)
      end do
      write(6,'(a,/)') '***********************************************************************'

      !second sum up all the centric and acentric reflections
      sec = 0; sea = 0; nsec = 0; nsea = 0
      do i = 1, nrefls
         if(hkls(i)%missing) cycle !not include missing reflections
	 call muleps(hkls(i)%h,hkls(i)%k,hkls(i)%l,mul,eps)
	 if(ifcenter(hkls(i)%h,hkls(i)%k,hkls(i)%l))then   !for centric reflections
	    do n = 1, 6
	       sec(n) = sec(n) + mul*(hkls(i)%enorm)**n
	    end do
	    do n = 1, 3
	       sec(n+6) = sec(n+6) + mul*(abs((hkls(i)%enorm)**2-1))**n
	    end do
	    nsec = nsec + mul
	 else                                              !for acentric reflections
	    do n = 1, 6
	       sea(n) = sea(n) + mul*(hkls(i)%enorm)**n
	    end do
	    do n = 1, 3
	       sea(n+6) = sea(n+6) + mul*(abs((hkls(i)%enorm)**2-1))**n
	    end do
	    nsea = nsea + mul
	 end if
      end do
      sec = sec / nsec;  sea = sea / nsea
      !now output centric and acentric reflections summary
      write(6,'(/,a,/)') '********   Average of all centric and acentric reflections    *******'
      write(6,'(a)')     '                 Centric                            Acentric         '
      write(6,'(a)')     '          Observed    Theoretical            Observed   Theoretical  '
      write(6,'(a,5X,F6.3,5X,F6.3,19X,F6.3,5X,F6.3)')' <E>  ', sec(1),tec(1),sea(1),tea(1)
      write(6,'(a,5X,F6.3,5X,F6.3,19X,F6.3,5X,F6.3)')' <E^2>', sec(2),tec(2),sea(2),tea(2)
      write(6,'(a,5X,F6.3,5X,F6.3,19X,F6.3,5X,F6.3)')' <E^3>', sec(3),tec(3),sea(3),tea(3)
      write(6,'(a,5X,F6.3,5X,F6.3,19X,F6.3,5X,F6.3)')' <E^4>', sec(4),tec(4),sea(4),tea(4)
      write(6,'(a,5X,F6.3,5X,F6.3,19X,F6.3,5X,F6.3)')' <E^5>', sec(5),tec(5),sea(5),tea(5)
      write(6,'(a,5X,F6.3,5X,F6.3,19X,F6.3,5X,F6.3)')' <E^6>', sec(6),tec(6),sea(6),tea(6)
      write(6,'(a,3X,F6.3,5X,F6.3,19X,F6.3,5X,F6.3)')' <|Z-1|>', sec(7),tec(7),sea(7),tea(7)
      write(6,'(a,1X,F6.3,5X,F6.3,19X,F6.3,5X,F6.3)')' <|Z-1|^2>', sec(8),tec(8),sea(8),tea(8)
      write(6,'(a,1X,F6.3,5X,F6.3,19X,F6.3,5X,F6.3)')' <|Z-1|^3>', sec(9),tec(9),sea(9),tea(9)
      write(6,'(a,I5,14X,a,I5)') ' Number of centric : ', nsec, ' Number of acentric : ', nsea
      write(6,'(a,/)')   '*********************************************************************'

      !third cumulative probabily distribution
      emax = 0
      do i = 1, nrefls
         if(hkls(i)%missing) cycle !not include missing reflections
	 emax = max(hkls(i)%enorm,emax)
      end do
      ncump = min(nint(emax*1000+0.5),3000)+1
      allocate(cump(ncump),cumpa(ncump),cumpc(ncump))
      cump = 0; cumpa = 0; cumpc = 0
      do i = 1, nrefls
         if(hkls(i)%missing) cycle !not include missing reflections
	 !sort enorms by 0.001
         n = min(nint(1000.*hkls(i)%enorm+0.5),3000)+1
	 if(n .gt. ncump) n = ncump
	 cump(n) = cump(n) + 1
	 if(ifcenter(hkls(i)%h,hkls(i)%k,hkls(i)%l))then
	    cumpc(n) = cumpc(n) + 1   !centric reflections
	 else
	    cumpa(n) = cumpa(n) + 1   !acentric reflections
	 end if
      end do
      !now make cumulative distribution(from large e to small e)
      do i = 1, ncump-1
	 n = ncump - i
         cump(n) = cump(n+1)+cump(n)
	 cumpa(n) = cumpa(n+1)+cumpa(n)
	 cumpc(n) = cumpc(n+1)+cumpc(n)
      end do
      !and visualize theoretical and observed cumulative probability
      write(6,'(/,a,/)') '**** Cumulative Probability Distribution Of Normalized Intensities *****'
      write(6,'(a)')     '                  Centric                          Acentric             '
      write(6,'(a)')     '    Z      Observed     Theoreical           Observed     Theoretical   '
      do i = 1, ncump, 50
         ee = (0.001 * (i-1))**2   !Z=E**2
	 tza = 1-exp(-ee)          !theoretical acentric cumulative probability
	 tzc = fgsl_sf_erf(dble(sqrt(ee/2))) !theory centric distribution
	 oza = 1 - real(cumpa(i))/cumpa(1)
	 ozc = 1 - real(cumpc(i))/cumpc(1)
	 write(6,'(2X,F6.4,3X,F6.4,9X,F6.4,13X,F6.4,9X,F6.4)') ee, ozc, tzc, oza, tza
      end do
      write(6,'(a,/)')   '************************************************************************'

      return
   end subroutine ecal


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Charge flipping : main program to implement phase retrieval algorithm
   !  Parameters :     frac : fractional proportion of density zone
   !                   niteration : number of iterations
   !                   ntry : number of different random seeds
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   subroutine charge_flipping
      implicit none
      integer n, tt, n_multi_real, s, x, y, z, xx, yy, zz, i, j, h, k, l, m, reci, nbest !reci表示倒空间中使用哪种约束，1表示使用Fo替代Fc，2表示使用2Fo-Fc替代Fc
      integer thre_method !thre_method为1使用吴金松的阈值设点方法，固定百分比；thre_method为2使用prasa的设定方法，deta(frac)*stddev
      complex(8) sf(0:nu-1,0:nv-1,0:nw-1), map(0:nu-1,0:nv-1,0:nw-1), map_prev(0:nu-1,0:nv-1,0:nw-1), map2(0:nu-1,0:nv-1,0:nw-1)
      complex(8) mapbest(0:nu-1,0:nv-1,0:nw-1), tmp
      real(8) denstmp(nu*nv*nw)
      real phase, tmp1, tmp2, cutoff, sum1, sum2, rfactor, rbest, cphase, skew_best, cc9, cc9_best, cc9_shelx,cc9_shelx_best
      real reci_thr_diff !表示the absolute value anologue of w parameter of the Fobs+DeltaF idea (Oszlanyi&Suto, 2008)
      type(reflection) :: hkls_best(nrefls),hkls_best_10(nrefls,10)
      logical center
      type(c_ptr) plan
      integer,allocatable :: multi_index_real(:,:)

      real(8), external :: ZBQLUAB

      real:: cc9_best_10(10)=-100, cc9_min=100
      integer:: index_min, iterate_num_best(10)
      character(100)::savename,tmp_savename
      
      integer ios
      character(len=20) filename  !计算的CC或者SKEW输出的路径
      character(len=20) str_m
      real(8) skew
      
      integer charge_flip !判断正空间约束的类型，0=off (RAAR); 1="standard"; 2=with the RAAR threshold; 3=flip-mem using beta (Oszlanyi&Suto, 2008)
      real beta !RAAR算法
      call ZBQLINI(0)

      !call wilson to divide out temperator factor effect
      !if(use_wilson .and. .not.use_ecal) call wilson(wilson_bins,wilson_res)

      !call ecal to calculate normalized structure factor
      !if(use_ecal) call ecal

      write(6,'(/,a,/)') '************************  Charge Flipping  ****************************'
      write(6,'(a)')     '     ntry  rfactor    ntry  rfactor    ntry  rfactor    ntry  rfactor'

      rbest = 1.0
      skew_best=0.0001
      cc9_best=0.0001
      cc9_shelx_best=0.0001
      nbest=0

      n_multi_real=0
      do z = 0, nw-1
      do y = 0, nv-1
      do x = 0, nu-1
        do s=2,nsymp
           xx = x*rsym(1,1,s)+y*rsym(1,2,s)+z*rsym(1,3,s)
           yy = x*rsym(2,1,s)+y*rsym(2,2,s)+z*rsym(2,3,s)
           zz = x*rsym(3,1,s)+y*rsym(3,2,s)+z*rsym(3,3,s)
           if (xx .lt. 0)  xx=xx+nu-1
           if (yy .lt. 0)  yy=yy+nv-1
           if (zz .lt. 0)  zz=zz+nw-1
           if(x.eq.xx .and. y.eq.yy .and. z.eq.zz) then
              n_multi_real = n_multi_real + 1
              exit
           end if
        end do
      end do
      end do
      end do
      !print*, n_multi_real
      
      allocate(multi_index_real( n_multi_real,3))
      tt=0
      do z = 0, nw-1
      do y = 0, nv-1
      do x = 0, nu-1
        do s=2,nsymp
           xx = x*rsym(1,1,s)+y*rsym(1,2,s)+z*rsym(1,3,s)
           yy = x*rsym(2,1,s)+y*rsym(2,2,s)+z*rsym(2,3,s)
           zz = x*rsym(3,1,s)+y*rsym(3,2,s)+z*rsym(3,3,s)
           if (xx .lt. 0)  xx=xx+nu-1
           if (yy .lt. 0)  yy=yy+nv-1
           if (zz .lt. 0)  zz=zz+nw-1
           if(x.eq.xx .and. y.eq.yy .and. z.eq.zz) then
              tt=tt+1
              multi_index_real(tt,1)=x
              multi_index_real(tt,2)=y
              multi_index_real(tt,3)=z
              exit
           end if
        end do
      end do
      end do
      end do
       print*, "the number of multiplicity point: ",size(multi_index_real)
      !print*, size(multi_index_real),multi_index_real(100,1),multi_index_real(100,2),multi_index_real(100,3)
      outer: do m = 1, ntry
 
         !allocate different random phases in asymmetric unit (centric reflections must be constrained)
         do i = 1, nrefls
	    hkls(i)%gnorm = 0
            if(.not. hkls(i)%missing)then
               hkls(i)%gnorm = hkls(i)%amp  !初始的计算振幅用实验振幅代替
            end if
            if(mtz_has_pha)then
               if( hkls(i)%missing)then
                  call ifcenter2(hkls(i)%h,hkls(i)%k,hkls(i)%l,center,cphase)
                  if(center)then   !for centric reflections
                     hkls(i)%opha = cphase
                  else
                     hkls(i)%opha = sngl(ZBQLUAB(dble(-pi),dble(pi)))
                  end if
               end if
            else
               call ifcenter2(hkls(i)%h,hkls(i)%k,hkls(i)%l,center,cphase)
               if(center)then   !for centric reflections
                  hkls(i)%opha = cphase
               else
                  hkls(i)%opha = sngl(ZBQLUAB(dble(-pi),dble(pi)))
               end if
            end if
         end do
         hkls(nrefls-1)%opha = 0    !leave phase of (000) to be 0
      
       !计算初始相位得到的map图
       reci=1 !1表示使用Fo替代Fc，2表示使用2Fo-Fc替代Fc
       reci_thr_diff=0 !表示the absolute value anologue of w parameter of the Fobs+DeltaF idea (Oszlanyi&Suto, 2008)
       sf = dcmplx(0,0)
       call ApplyReciRestr(sf,reci,reci_thr_diff) !sf(0,0,0)设为0

       !next calculate density map
       plan = fftw_plan_dft_3d(nw,nv,nu,sf,map,FFTW_FORWARD,FFTW_ESTIMATE)
       call fftw_execute_dft(plan,sf,map)
       map = map/vol
      !t=0
      !do k = 0, nw-1
      !do j = 0, nv-1
      !do i = 0, nu-1
      !  if (aimag(map(i,j,k)) .gt. 0.1) then
      !     t=t+1
      !  end if
      !end do
      !end do
      !end do
      !print*,"虚部大于0.1的个数：",t
      !把迭代过程中的每一个CC或者skew存储到filename文件中  
      write(str_m,'(I10)') m
      filename="file_"//trim(adjustl(str_m))//".txt"
      open(unit=13, file=filename, iostat=ios, status="new")
      if ( ios /= 0 ) stop "Error opening file name"

      !iterations begins here
      inner: do n = 1, niteration
         
         !把上一轮的map存储起来
         map_prev=map      

         !倒空间约束
         !first expand asymmetric unit to whole space 
	 !set outer resolution and systematic absent reflections always to be 0, else comment next line to retain last cycled sf
         reci=1 !1表示使用Fo替代Fc，2表示使用2Fo-Fc替代Fc
         reci_thr_diff=0 !表示the absolute value anologue of w parameter of the Fobs+DeltaF idea (Oszlanyi&Suto, 2008)
	 sf = dcmplx(0,0)
         call ApplyReciRestr(sf,reci,reci_thr_diff) !sf(0,0,0)设为0

	 !next calculate density map
	 plan = fftw_plan_dft_3d(nw,nv,nu,sf,map,FFTW_FORWARD,FFTW_ESTIMATE)
	 call fftw_execute_dft(plan,sf,map)
	 map = map/vol
         
         !计算阈值
         thre_method=2 !thre_method为1使用吴金松的阈值设点方法，固定百分比；thre_method为2使用prasa的设定方法，deta(frac)*stddev
         if ( (n .lt. 10 ) .or. ((n .lt. 50) .and. (mod(n,3) .eq. 0)) .or. (mod(n,8) .eq. 0) ) then
            call  threshold_set(map,denstmp,cutoff,thre_method)
            !print*,"cutoff： ", cutoff
         end if
         !print*,"threshold: ",cutoff
         
         !把位于对称元素上的点的密度设为0
         do s=1,n_multi_real
            x=multi_index_real(s,1)
            y=multi_index_real(s,2)
            z=multi_index_real(s,3)
            map(x,y,z)=dcmplx(0,0)
         end do
         
         !生成一个新的map图map2用于计算相关系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         map2=dcmplx(0,0)
         do k = 0, nw-1
         do j = 0, nv-1
         do i = 0, nu-1
            if (cdabs(map(i,j,k)) .lt. 4.5*cutoff) then
               map2(i,j,k)=dcmplx(0.,0.)
            else
               map2(i,j,k)=dcmplx( 2*real(map(i,j,k))-6.5*cutoff , 0. )
            end if
         end do
         end do
         end do
         plan = fftw_plan_dft_3d(nw,nv,nu,map2,sf,FFTW_BACKWARD,FFTW_ESTIMATE)
         call fftw_execute_dft(plan,map2,sf)
         sf = sf/nu/nv/nw*vol

         do i = 1, nrefls
             tmp=sf(lmod(hkls(i)%h,nu),lmod(hkls(i)%k,nv),lmod(hkls(i)%l,nw))
             hkls(i)%gnorm = cdabs(tmp)
             hkls(i)%opha = atan2(aimag(tmp),real(tmp))
         end do
         
         call pearson(cc9)
         call cc(cc9_shelx)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	 !now enforce real space constraint
         beta=0.82
         charge_flip=0
         call ApplyDireRestr(map,map_prev,cutoff,charge_flip,beta)

	 !and ifft to get reciprocal space
	 plan = fftw_plan_dft_3d(nw,nv,nu,map,sf,FFTW_BACKWARD,FFTW_ESTIMATE)
	 call fftw_execute_dft(plan,map,sf)
	 sf = sf/nu/nv/nw*vol
         


	 !finally extract asymmetrical unit and save in opha, gnorm
         !do i = 1, nrefls
         !   tmp = dcmplx(0,0)
         !do j = 1, nsymp
         !   h = hkls(i)%h*rsym(1,1,j)+hkls(i)%k*rsym(2,1,j)+hkls(i)%l*rsym(3,1,j)
         !   k = hkls(i)%h*rsym(1,2,j)+hkls(i)%k*rsym(2,2,j)+hkls(i)%l*rsym(3,2,j)
         !   l = hkls(i)%h*rsym(1,3,j)+hkls(i)%k*rsym(2,3,j)+hkls(i)%l*rsym(3,3,j)
	 !   phase=2*pi*(hkls(i)%h*rsym(1,4,j)+hkls(i)%k*rsym(2,4,j)+hkls(i)%l*rsym(3,4,j))
	 !   tmp = tmp + sf(lmod(h,nu),lmod(k,nv),lmod(l,nw))*dcmplx(cos(phase),sin(phase))
	 !end do
         !   tmp = tmp / nsymp
	 !   hkls(i)%gnorm = cdabs(tmp)
	 !   hkls(i)%opha = atan2(aimag(tmp),real(tmp))
         do i = 1, nrefls
             tmp=sf(lmod(hkls(i)%h,nu),lmod(hkls(i)%k,nv),lmod(hkls(i)%l,nw))
             hkls(i)%gnorm = cdabs(tmp)
             hkls(i)%opha = atan2(aimag(tmp),real(tmp))
         end do
            !enforce centric constraint
!            call ifcenter2(hkls(i)%h,hkls(i)%k,hkls(i)%l,center,cphase)
!            if(center)then
!               if(acos(cos(hkls(i)%opha-cphase)) .lt. acos(cos(hkls(i)%opha-cphase-pi)))then
!                  hkls(i)%opha = cphase
!               else
!                  hkls(i)%opha = cphase+pi
!                  if(hkls(i)%opha .gt. 2*pi) hkls(i)%opha = hkls(i)%opha-2*pi
!               end if
!            end if
         !end do
         
         call skewness(denstmp,skew)
         !call pearson(cc9)
         !call cc(cc9)
         write(13,'(F7.4,2X,F7.4)') cc9_shelx,cc9
      end do inner
         close(13)
	 !calculate error metric
	 sum1 = 0; sum2 = 0
	 do i = 1, nrefls
            if(hkls(i)%missing) cycle !not include missing reflections
            if(use_ecal)then          !use normalized sf
	       sum1 = sum1 + abs(hkls(i)%enorm-hkls(i)%gnorm)
	       sum2 = sum2 + hkls(i)%enorm
            else
	       sum1 = sum1 + abs(hkls(i)%amp-hkls(i)%gnorm)
	       sum2 = sum2 + hkls(i)%amp
            end if
	 end do
	 rfactor = sum1/sum2

	 if(mod(m,4) .eq. 0)then
	    !write(6,'(4X,I4,2X,F7.4)') m, rfactor
            write(6,'(4X,I4,2X,F7.4,2X,F7.4,2X,F7.4)') m,cc9_shelx,rfactor,cc9
	 else
	    !write(6,'(4X,I4,2X,F7.4)',advance='no') m, rfactor
            write(6,'(4X,I4,2X,F7.4,2X,F7.4,2X,F7.4)') m,cc9_shelx,rfactor,cc9
	 end if
                
         if(cc9 .gt. cc9_best)then
            nbest=m
            skew_best = skew
            cc9_shelx_best=cc9_shelx
            cc9_best = cc9
            rbest = rfactor
            mapbest = map
            hkls_best = hkls(1:nrefls)
         end if

         !输出最好CC——pearson中的前10个
         cc9_min=100
         do i=1,10
            if (cc9_best_10(i)<cc9_min) then
               index_min=i
               cc9_min=cc9_best_10(i)
            end if
         end do
         !print*,index_min

         if (cc9>cc9_min) then
            cc9_best_10(index_min)=cc9
            hkls_best_10(:,index_min)=hkls(1:nrefls)
            iterate_num_best(index_min)=m
         end if
         !print*,cc9_best_10

         if(cc9 .gt. cc9_best)then
            skew_best = skew
            cc9_shelx_best=cc9_shelx
            cc9_best = cc9
            rbest = rfactor
            mapbest = map
            hkls_best = hkls(1:nrefls)
         end if
      end do outer
      !把最好的结果写入mtz
      do i=1,10
         write(6,'(a,F7.4)') 'top best 10 CC is : ', cc9_best_10(i)
         write(str_m,'(I10)') iterate_num_best(i)
         savename="recover_"//trim(adjustl(str_m))//".mtz"
         call mtzout(hkls_best_10(:,i),savename)
      end do


      write(6,'(a,I4,F8.4,F8.4,F8.4)') 'Best CC is : ', nbest,cc9_shelx_best,rbest,cc9_best
      print*, "charge_flip：",charge_flip,"阈值设定方法: ",thre_method
      write(6,'(a,/)') '**********************************************************************'

      !call ccp4mapout('MAPOUT','recovered density map',sngl(real(mapbest)),nspgrp,cell,nu,nv,nw,&
      !   0,0,0,nu-1,nv-1,nw-1)

      !call mtzout(hkls_best)

      return
   end subroutine charge_flipping
   

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !  resolution  : compute q=1/d=2*sin(theta)/lambda of hkl and volume
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   subroutine resolution(h,k,l,res,vol)
      implicit none
      integer h, k, l
      real res, vol, a, b, c
      real rad1,rad2,rad3,tmp,numerator
      real s11, s22, s33, s12, s23, s13

      a = cell(1); b = cell(2); c = cell(3)

      rad1 = cell(4)/180.*pi
      rad2 = cell(5)/180.*pi
      rad3 = cell(6)/180.*pi

      tmp = cos(rad1)**2+cos(rad2)**2+cos(rad3)**2-2*cos(rad1)*cos(rad2)*cos(rad3)

      vol = a*b*c*sqrt(1-tmp)

      s11 = ( b*c*sin(rad1) )**2
      s22 = ( a*c*sin(rad2) )**2
      s33 = ( a*b*sin(rad3) )**2
      s12 = a*b*c*c*(cos(rad1)*cos(rad2) - cos(rad3))
      s23 = a*a*b*c*(cos(rad2)*cos(rad3) - cos(rad1))
      s13 = a*b*b*c*(cos(rad3)*cos(rad1) - cos(rad2))

      numerator = s11*h*h + s22*k*k + s33*l*l + 2*s12*h*k + 2*s23*k*l + 2*s13*h*l
      
      res = sqrt(numerator)/vol

      return
   end subroutine resolution

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !muleps  :  compute unique equivalent reflections excluding multiplicity and epsilon 
   !           factor for reflection(h k l)
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   subroutine muleps(h,k,l,mul,eps)
      implicit none
      integer h, k, l, mul, eps, i, j
      integer hh, kk, ll
      mul = 1; eps = 1
      do j = 2, nsymp
         hh = h*rsym(1,1,j)+k*rsym(2,1,j)+l*rsym(3,1,j)
         kk = h*rsym(1,2,j)+k*rsym(2,2,j)+l*rsym(3,2,j)
         ll = h*rsym(1,3,j)+k*rsym(2,3,j)+l*rsym(3,3,j)
	 if(h.eq.hh .and. k.eq.kk .and. l.eq.ll) eps = eps + 1
	 if((h.eq.hh .and. k.eq.kk .and. l.eq.ll).or.(h.eq.-hh .and. k.eq.-kk .and. l.eq.-ll)) mul = mul + 1
      end do
      mul = nsymp/mul !leave out multiplicity, and mul means number of equivalent reflections with different indexes 
      return
   end subroutine muleps


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !ifcenter : test wheter a reflection is centric(true) or acentric(false)
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   function ifcenter(h,k,l)
      implicit none
      logical ifcenter
      integer h, k, l, hh, kk, ll, j
      ifcenter = .false.
      do j = 2, nsymp
         hh = h*rsym(1,1,j)+k*rsym(2,1,j)+l*rsym(3,1,j)
         kk = h*rsym(1,2,j)+k*rsym(2,2,j)+l*rsym(3,2,j)
         ll = h*rsym(1,3,j)+k*rsym(2,3,j)+l*rsym(3,3,j)
	 if(h.eq.-hh .and. k.eq.-kk .and. l.eq.-ll)then
	    ifcenter = .true.
	    return
	 end if
      end do
      return
   end function ifcenter


   subroutine ifcenter2(h,k,l,center,cphase)
      implicit none
      logical center
      integer h, k, l, hh, kk, ll, j, isym
      real cphase
      real(8), external :: ZBQLUAB

      center = .false.
      cphase = 0

      do j = 2, nsymp
         hh = h*rsym(1,1,j)+k*rsym(2,1,j)+l*rsym(3,1,j)
         kk = h*rsym(1,2,j)+k*rsym(2,2,j)+l*rsym(3,2,j)
         ll = h*rsym(1,3,j)+k*rsym(2,3,j)+l*rsym(3,3,j)
	 if(h.eq.-hh .and. k.eq.-kk .and. l.eq.-ll)then
	    center = .true.
            isym = j
	    exit
	 end if
      end do

      !phase restrictions for centric reflection: pi*h*T or pi*h*T+pi
      if(center)then
         if(ZBQLUAB(dble(0),dble(1)) .lt. 0.5)then
             cphase = pi*(h*rsym(1,4,isym)+k*rsym(2,4,isym)+l*rsym(3,4,isym))
         else
             cphase = pi*(h*rsym(1,4,isym)+k*rsym(2,4,isym)+l*rsym(3,4,isym))+pi
         end if
         !restrict cphase within (0-360)
         if(cphase .gt. 2*pi)then
             cphase = cphase - int(cphase/pi/2)*pi*2
         else if(cphase .lt. 0)then
             cphase = cphase - nint(cphase/pi/2)*pi*2
         end if
      end if
      return
   end subroutine ifcenter2


   integer function lmod(i,j)
      implicit none
      integer i, j
      if(i .ge. 0)then
	 lmod = i
      else
	 lmod = i + j
      end if
      return
   end function lmod

   recursive subroutine quicksort(a,m,n,f,s,e)
      implicit none
      integer::m,n,f,s,e!!!!!m,n are the array row and col,f is the col to
      !! sort ,s and e are the initial and end of sorting.
      real(kind=8)::a(m,n)
      integer::l,r
      real::key
      real ::temp(m)
      l=s
      r=e+1

      if(r<=l) return
      !!!!!start 
      key=a(f,s)
      do while(.true.)
       do while(.true.)
         l=l+1
         if((a(f,l)>key.or.(l>=e))) exit
        end do
       do while(.true.)
         r=r-1
         if((a(f,r)<key.or.(r<=s))) exit
         end do
       if(r<=l) exit
         temp=a(1:m,l)
         a(1:m,l)=a(1:m,r)
         a(1:m,r)=temp
        end do
         temp=a(1:m,s)
         a(1:m,s)=a(1:m,r)
         a(1:m,r)=temp
         call quicksort(a,m,n,f,s,r-1)
         call quicksort(a,m,n,f,r+1,e)
      return
   end subroutine quicksort

!========================================================================
!write out the recovered amplitude and phase together with the observed
!========================================================================
   subroutine mtzout(hkls_best,filename)
      implicit none
      integer mtzind
      character(5) lsprg0(200)!column labels denoting each column of mtz
      character ctprg0(200)   !single character column type
      real adata(6)
      type(reflection) :: hkls_best(nrefls)
      integer i
      character(100)::filename

      mtzind = 1
      data lsprg0 /'H','K','L','FO','FC','PHIC',194*' '/
      data ctprg0 /'H','H','H','F','F','P',194*' '/

      call lwopen(mtzind,filename)
      call lwtitl(mtzind,'recovered mtz',0)
      call lwcell(mtzind,cell)
      call lwclab(mtzind,lsprg0,6,ctprg0,0)
      call lwsymm(mtzind,nsym,nsymp,rsym,ltype,nspgrp,spgrp,nampg)

      do i = 1, nrefls
         adata(1:3) = (/hkls_best(i)%h, hkls_best(i)%k, hkls_best(i)%l/)
         adata(4) = hkls_best(i)%amp
         adata(5) = hkls_best(i)%gnorm
         adata(6) = hkls_best(i)%opha/pi*180
         if(adata(6) .lt. 0) adata(6) = adata(6)+360
         call lwrefl(mtzind,adata)
      end do

      call lwclos(mtzind,0)

      return
   end subroutine mtzout


!==============================================
!write density into .map
!==============================================
   subroutine ccp4mapout(name,title,map,nspgrp,cell,&
         nu,nv,nw,nu1,nv1,nw1,nu2,nv2,nw2)
      implicit none
      integer nspgrp,nu,nv,nw,nu1,nv1,nw1,nu2,nv2,nw2
      real cell(6)
      real map(0:nu-1,0:nv-1,0:nw-1)
      character(len=*) name, title
      real lsec(0:100000)
      integer ifast,imedm,islow,lfast,lmedm,lslow
      integer m1(3),m2(3),ifms(3),jfms(3),juvw(3),mxyz(3)
      integer i,iu,iv,iw
      integer nsec,mode,nsym,nsymp
      real rsym(4,4,192)
      character(len=10) namspg,nampg

      integer axis(230)
      data axis/2,2,2,2,1,1,1,1,1,2,1,1,1,1,1,2,2,2,1,2,2,1,2,207*1/

      if(axis(mod(nspgrp,1000)).eq.1) then
         jfms(1)=2
         jfms(2)=1
         jfms(3)=3
      else
         jfms(1)=3
         jfms(2)=1
         jfms(3)=2
      end if

      do i=1,3
         juvw(jfms(i))=i
      end do

      mxyz(1)=nu
      mxyz(2)=nv
      mxyz(3)=nw
      m1(juvw(1))=nu1
      m2(juvw(1))=nu2
      m1(juvw(2))=nv1
      m2(juvw(2))=nv2
      m1(juvw(3))=nw1
      m2(juvw(3))=nw2
      nsec=m2(3)-m1(3)+1
      mode=2

      call msymlb(7,nspgrp,namspg,nampg,nsym,nsymp,rsym)

      call mwrhdl(7,name,title,nsec,jfms,mxyz,m1(3),m1(1),&
         m2(1),m1(2),m2(2),cell,nspgrp,mode)
      call msywrt(7,nsym,rsym)

      lfast=m2(1)-m1(1)+1
      lmedm=m2(2)-m1(2)+1
      lslow=m2(3)-m1(3)+1

      if(lfast*lmedm.gt.100000)&
         call ccperr(1,' ccp4mapout - Mask section > lsec: recompile')
      do islow=0,lslow-1
         ifms(3)=islow+m1(3)
         do imedm=0,lmedm-1
            ifms(2)=imedm+m1(2)
            do ifast=0,lfast-1
               ifms(1)=ifast+m1(1)
               iu=ifms(juvw(1))
               iv=ifms(juvw(2))
               iw=ifms(juvw(3))
               lsec(ifast+lfast*imedm)=map(iu,iv,iw)
            end do
         end do
         call mspew(7,lsec)
      end do
      call mwclose(7)

      return
   end subroutine ccp4mapout
        
   subroutine skewness(density,skew)
        !根据密度图来判断相位的准确性，使用skewness，参考文献Acta Cryst. (2013). D69, 2039–2049
        implicit none
        real(8) density(nu*nv*nw) !将三维密度图变成一个一维图
        real(8) sum_density_2,sum_density_3,stddev,mean,skew
        integer i,n,shape_density(1:rank(density))

        !do i = 1, 10*10*10
        !    density(i)=i
        !end do
        
        !计算标准差
        shape_density=shape(density)

        mean = sum(density) / shape_density(1)
        stddev = sqrt(sum((density - mean)**2) / shape_density(1))
        
        !计算skew，去除<5sigma和>5sigma的密度
        n=0
        sum_density_2=0
        do i=1,nu*nv*nw
            if (density(i)>-5*stddev .and. density(i)<5*stddev) then
                sum_density_2=sum_density_2+density(i)**2
                sum_density_3=sum_density_3+density(i)**3
                n=n+1
            end if
        end do

        skew=(sum_density_3/n)/(sqrt(sum_density_2/n))**3
        
        return
    end subroutine skewness
    

    subroutine cc(cc9)
        !使用cc9来计算相关系数，参考文献：Acta Cryst. (2002). D58, 1772±1779
        implicit none
        integer ngood,i     !表示观测到的衍射点的个数，即观测的衍射强度不是?
        real wxy,w,wx,wy,wx2,wy2,wtmp,cc9

        ngood=0 ;wx=0; wy=0; wxy=0; wx2=0; wy2=0; w=0
        do i=1,nrefls
            if(.not. hkls(i)%missing) then
                wtmp=1/(1+0.1*(hkls(i)%sigf)*(hkls(i)%sigf))  !计算w_array = 1 / (1 + 0.1 * sigFo_array * sigFo_array)
                w=w+wtmp
                wx=wx+wtmp*hkls(i)%amp
                wy=wy+wtmp*hkls(i)%gnorm
                wxy=wxy+wtmp*(hkls(i)%amp)*(hkls(i)%gnorm)
                wx2=wx2+wtmp*(hkls(i)%amp)*(hkls(i)%amp)
                wy2=wy2+wtmp*(hkls(i)%gnorm)*(hkls(i)%gnorm)
                ngood=ngood+1
            end if
        end do

        cc9=100*(wxy*w - wx*wy)/sqrt((wx2*w - wx*wx)*(wy2*w - wy*wy))
       return
    end subroutine cc
    
   subroutine pearson(pearson_correl)
        !使用pearson来计算相关系数
        implicit none
        integer ngood,i,w,wtmp     !表示观测到的衍射点的个数，即观测的衍射强度不是?
        real wxy,wx,wy,wx2,wy2,xy,x,y,x2,y2,correl,correlw,pearson_correl

        ngood=0 ;wx=0; wy=0; wxy=0; wx2=0; wy2=0; w=0; x=0; y=0; xy=0; x2=0; y2=0
        do i=1,nrefls
           if ((.not. hkls(i)%missing) .and. (hkls(i)%sigf .gt. 0.0001) .and. (hkls(i)%amp .lt. 2.1) &
                   .and. (hkls(i)%amp .gt. 0.0001)) then
                wtmp=1
                w=w+wtmp
                wx=wx+wtmp*hkls(i)%amp
                wy=wy+wtmp*hkls(i)%gnorm
                wxy=wxy+wtmp*(hkls(i)%amp)*(hkls(i)%gnorm)
                wx2=wx2+wtmp*(hkls(i)%amp)*(hkls(i)%amp)
                wy2=wy2+wtmp*(hkls(i)%gnorm)*(hkls(i)%gnorm)
                
                if ( hkls(i)%sigf .lt. 0.608 ) then
                   ngood=ngood+1
                   x=x+hkls(i)%amp
                   y=y+hkls(i)%gnorm
                   xy=xy+(hkls(i)%amp)*(hkls(i)%gnorm)
                   x2=x2+(hkls(i)%amp)*(hkls(i)%amp)
                   y2=y2+(hkls(i)%gnorm)*(hkls(i)%gnorm)
                end if
           end if
        end do

        correlw=100*(wxy*w - wx*wy)/sqrt((wx2*w - wx*wx)*(wy2*w - wy*wy))
        correl=100*(ngood*xy - x*y)/sqrt((ngood*x2 - x*x)*(ngood*y2 - y*y))
        pearson_correl=(correlw+correl)/2
       return
    end subroutine pearson

    !阈值设定方法
    subroutine threshold_set(map,denstmp,cutoff,thre_method)
         implicit none
         integer::i,j,k,h,thre_method
         !thre_method为1使用吴金松的阈值设点方法，固定百分比；thre_method为2使用prasa的设定方法，deta(frac)*stddev
         real::cutoff,sum_denstmp,sum_denstmp_2
         real(8)::denstmp(nu*nv*nw),denstmp_2(nu*nv*nw)
         complex(8)::map(0:nu-1,0:nv-1,0:nw-1)
         h = 0
         sum_denstmp=0; sum_denstmp_2=0
         !now sort map in ascending sequence
         do k = 0, nw-1
         do j = 0, nv-1
         do i = 0, nu-1
            h = h + 1
            denstmp(h) = real(map(i,j,k))
            denstmp_2(h)=real(map(i,j,k)) * real(map(i,j,k))
            sum_denstmp=sum_denstmp+real(map(i,j,k))
            sum_denstmp_2=sum_denstmp_2+real(map(i,j,k))*real(map(i,j,k))
         end do
         end do
         end do

        !使用fsgl排序方法，吴金松文献的阈值设定方法
         if (thre_method .eq. 1) then
            !call fgsl_sort(denstmp,1,nu*nv*nw)
            call quicksort(denstmp,1,nu*nv*nw,1,1,nu*nv*nw)
            cutoff = denstmp(nint((1-frac)*nu*nv*nw))
         end if

         !使用prasa的阈值设定方法
         if (thre_method .eq. 2) then
            !denstmp_2=denstmp**2
            !cutoff=sqrt( sum(denstmp_2)/h - (sum(denstmp)/h)**2 )*frac
            cutoff=sqrt( sum_denstmp_2/h - (sum_denstmp/h)*(sum_denstmp/h) )*frac
            !print*,"mean of map is: ",sum_denstmp/h
         end if
         return
    end subroutine threshold_set

    !倒空间限制
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ApplyReciRestr(sf,reci,reci_thr_diff)
        implicit none
        integer::i,j,h,k,l,reci !reci表示倒空间中使用哪种约束，1表示使用Fo替代Fc，2表示使用2Fo-Fc替代Fc
        complex(8)::sf(0:nu-1,0:nv-1,0:nw-1)
        real::newf,phase,tmp1,tmp2,reci_thr_diff !表示the absolute value anologue of w parameter of the Fobs+DeltaF idea (Oszlanyi&Suto, 2008)
        sf=dcmplx(0,0)
        !print*,reci,reci_thr_diff
        do i = 1, nrefls-1
        do j = 1, nsymp
            h = hkls(i)%h*rsym(1,1,j)+hkls(i)%k*rsym(2,1,j)+hkls(i)%l*rsym(3,1,j)
            k = hkls(i)%h*rsym(1,2,j)+hkls(i)%k*rsym(2,2,j)+hkls(i)%l*rsym(3,2,j)
            l = hkls(i)%h*rsym(1,3,j)+hkls(i)%k*rsym(2,3,j)+hkls(i)%l*rsym(3,3,j)
            phase= hkls(i)%opha-2*pi*(hkls(i)%h*rsym(1,4,j)+hkls(i)%k*rsym(2,4,j)+hkls(i)%l*rsym(3,4,j))
            tmp1=0; tmp2=0; newf=0
            if (.not. hkls(i)%missing) then
                !表示使用2Fo-Fc替代Fc
                if (reci .eq. 2) then 
                    !进行2F0-Fc的约束
                    if (reci_thr_diff >0) then
                        newf=2*hkls(i)%amp-hkls(i)%gnorm
                        if ( newf> (hkls(i)%amp+reci_thr_diff) ) then
                            newf=hkls(i)%amp+reci_thr_diff
                        else if (newf<0.000001) then
                            newf=0.000001
                        else if (newf<hkls(i)%amp-reci_thr_diff) then
                            newf=hkls(i)%amp-reci_thr_diff
                        end if
                        tmp1=newf*cos(phase)
                        tmp2=newf*sin(phase)
                    else 
                        !print*,"good"
                        newf=2*hkls(i)%amp-hkls(i)%gnorm
                        tmp1=newf*cos(phase)
                        tmp2=newf*sin(phase)
                    end if
                end if
                !表示使用Fo替代Fc
                if (reci .eq. 1) then
                    tmp1=hkls(i)%amp*cos(phase)
                    tmp2=hkls(i)%amp*sin(phase)
                end if
           else
             tmp1 = hkls(i)%gnorm*cos(phase)
             tmp2 = hkls(i)%gnorm*sin(phase)
           end if
            sf(lmod(h,nu),lmod(k,nv),lmod(l,nw)) = dcmplx(tmp1, tmp2)
            sf(lmod(-h,nu),lmod(-k,nv),lmod(-l,nw)) = dcmplx(tmp1,-tmp2)
        end do
        end do
        !sf(0,0,0) = dcmplx(hkls(nrefls)%gnorm,0) !for (000) reflection
        !sf(0,0,0) = dcmplx(hkls(nrefls)%gnorm*cos(hkls(i)%opha),hkls(nrefls)%gnorm*sin(hkls(i)%opha))
        return
    end subroutine ApplyReciRestr
    
    !应用正空间限制
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ApplyDireRestr(map,map_prev,cutoff,charge_flip,beta)
        implicit none
        complex(8)::map(0:nu-1,0:nv-1,0:nw-1), map_prev(0:nu-1,0:nv-1,0:nw-1)
        real::cutoff, beta, cmp
        integer::charge_flip, i, j, k

        !now enforce real space constraint
         do k = 0, nw-1
         do j = 0, nv-1
         do i = 0, nu-1
            !1=使用最原始的charge flipping算法，使用最原始的判据
            if (charge_flip .eq. 1) then
                if ( real(map(i,j,k)) .lt. cutoff ) then  !none density zone
                    map(i,j,k) = -map(i,j,k)
                end if
            end if

            !2=使用RAAR判据的RAAR算法
            if (charge_flip .eq. 2) then
                cmp=real(map(i,j,k)) + real(map(i,j,k)) - real(map_prev(i,j,k)) !使用RAAR判据
                if ( cmp .lt. cutoff ) then
                    map(i,j,k) = -map(i,j,k)
                end if
            end if

            !3=flip-mem using beta (Oszlanyi&Suto, 2008)
            if (charge_flip .eq. 3) then
                if ( real(map(i,j,k)) .lt. cutoff ) then  !none density zone
                    map(i,j,k) = -map(i,j,k)
                else
                    map(i,j,k)=map(i,j,k)+beta*(map(i,j,k)-map_prev(i,j,k))
                end if
            end if

            !0=RAAR
            if (charge_flip .eq. 0) then
                cmp=cdabs(map(i,j,k)) + cdabs(map(i,j,k)) - cdabs(map_prev(i,j,k)) !使用RAAR判据
                if ( cmp .lt. cutoff ) then
                    !map(i,j,k)=map(i,j,k)-beta*dcmplx(cmp,0)
                    map(i,j,k)=dcmplx(cdabs(map(i,j,k))-beta*cmp,0)
                end if
            end if
            map(i,j,k) = dcmplx(cdabs(map(i,j,k)),0)   !reality constraint, positivity，正性约束
         end do
         end do
         end do

        return
    end subroutine ApplyDireRestr

end module functions
