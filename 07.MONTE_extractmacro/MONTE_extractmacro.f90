module MONTE_extractmacro
   !-------------------------------------------------------------------------------------------
   ! Extract the macroscopic cross-sections from an input file 
   !-------------------------------------------------------------------------------------------
   use, intrinsic :: iso_fortran_env
   use globals
   use csv_nd_with_index
   
   implicit none

   contains

!----------------------------------------------------------------------------------------------
! Reads input deck and assign initial values
!----------------------------------------------------------------------------------------------
subroutine INP_read_input(core, mat, iso, signal, stack)

   integer :: i, j, k, numlin, ifound, istk, iz
   integer :: nx = 0
   integer :: nz = 0
   integer :: ios
   real(real64) :: r
   real(real64), allocatable :: tmpr(:)
   character :: c
   character*3 :: keyword
   character*100 :: w1, w2, w3
   character*100, allocatable :: tmp(:)
   character(:), allocatable :: line
   character(len=300) :: buffer
   logical :: found
   type(core_template), intent(inout) :: core
   type(mat_template), allocatable, intent(inout) :: mat(:)
   type(iso_template), allocatable, intent(inout) :: iso(:)
   type(signal_template), allocatable, intent(inout) :: signal(:)
   type(stack_template), allocatable, intent(inout) :: stack(:)
   type(mat_template), allocatable :: tmp_mat(:)
   type(iso_template), allocatable :: tmp_iso(:)
   type(signal_template), allocatable :: tmp_sgn(:)
   type(stack_template), allocatable :: tmp_stk(:)
   type(gen_template) :: gen

   write(800,*)'BELOW IS A COPY OF THE INPUT FILE'
   write(800,*)gen%separator

   !  read input line by line and write it to output
   open(700, file='input')
      do
         read(700,'(A)',iostat=ios) buffer
         if (ios .ne. 0) exit
         line = trim(buffer)
         write(800, '(A)') line
      end do
      write(800,*)gen%separator
   close(700)
   
   
   open(700, file='input')
   !     first symbol of the card
      read(700,'(a1)',err=11, end=11)c
   !     input deck line number
      numlin = 1
   !     $ - terminal card
      do while(c .ne. '$')
   !--------------------------------------------------------------------------------------------------
   !        * - comment card
         if(c .ne. '*') then
            backspace 700
   !           card keyword
            read(700,*,err=11)keyword

   !--------------------------------------------------------------------------------------------------
   !           path to XS data
            if(keyword .eq. 'XSP')then
               backspace 700
               read(700,*,err=11)keyword, w1
               core%path2data = trim(w1)
   
   !--------------------------------------------------------------------------------------------------
   !           core material:
   !              w1 : material name
   !              w2 : isotope name
   !              w3 : name of signal containing temperature (K)
            else if(keyword .eq. 'MAT')then
               backspace 700
               read(700,*,err=11)keyword, w1, w2, r, w3

   !              process isotope
               found = .false.
               do i=1,core%niso
                  if(trim(w2) .eq. core%isoname(i))found = .true.
               end do
               if(.not. found)then
   !                 add isotope name to the array
                  core%niso = core%niso + 1
                  allocate(tmp(core%niso))
                  if(core%niso .gt. 1) tmp(1:core%niso-1) = core%isoname
                  tmp(core%niso) = trim(w2)
                  call move_alloc(tmp, core%isoname)

   !                 add instance to iso
                  allocate(tmp_iso(core%niso))
                  if(core%niso .gt. 1)tmp_iso(1:core%niso-1) = iso
                  tmp_iso(core%niso)%isoname = trim(w2)
                  call move_alloc(tmp_iso, iso)
               end if

   !              process material
               found = .false.
               do i=1,core%nmat
                  if(trim(w1) .eq. core%matname(i))then
                     found = .true.
                     ifound = i
                  end if
               end do
               if(.not. found)then
   !                 add material name to the array core%matname
                  core%nmat = core%nmat + 1
                  allocate(tmp(core%nmat))
                  if(core%nmat .gt. 1)tmp(1:core%nmat-1) = core%matname
                  tmp(core%nmat) = trim(w1)
                  call move_alloc(tmp, core%matname)

   !                 add instance to mat
                  allocate(tmp_mat(core%nmat))
                  if(core%nmat .gt. 1) tmp_mat(1:core%nmat-1) = mat
                  tmp_mat(core%nmat)%matname = trim(w1)
                  call move_alloc(tmp_mat, mat)
                  ifound = core%nmat
               end if
   !              add isotope name to the array
               mat(ifound)%niso = mat(ifound)%niso + 1
               allocate(tmp(mat(ifound)%niso))
               if(mat(ifound)%niso .gt. 1) tmp(1:mat(ifound)%niso-1) = mat(ifound)%isoname
               tmp(mat(ifound)%niso) = trim(w2)
               call move_alloc(tmp, mat(ifound)%isoname)
   !              add isotope temperature to the array
               allocate(tmp(mat(ifound)%niso))
               if(mat(ifound)%niso .gt. 1)tmp(1:mat(ifound)%niso-1) = mat(ifound)%tsignal
               tmp(mat(ifound)%niso) = trim(w3)
               call move_alloc(tmp, mat(ifound)%tsignal)
   !              add number density to the array
               allocate(tmpr(mat(ifound)%niso))
               if(mat(ifound)%niso .gt. 1) tmpr(1:mat(ifound)%niso-1) = mat(ifound)%numden
               tmpr(mat(ifound)%niso) = r
               call move_alloc(tmpr, mat(ifound)%numden)


   !--------------------------------------------------------------------------------------------------
   !           signal
            else if(keyword .eq. 'SIG')then
               backspace 700
               read(700,*,err=11)keyword, w1

   !              process signal
               found = .false.
               if(allocated(signal))then
                  do i=1,size(signal)
                     if(trim(w1) .eq. signal(i)%sgnname) found = .true.
                  end do
               else
                  allocate(signal(0))
               end if
               if(.not. found)then
   !                 add instance to signal
                  allocate(tmp_sgn(size(signal)+1))
                  tmp_sgn(1:size(signal)) = signal
                  backspace 700
                  read(700,*,err=11)keyword, w1, w2
                  tmp_sgn(size(signal)+1)%sgnname = trim(w1)
                  tmp_sgn(size(signal)+1)%sgntype = trim(w2)
                  if(trim(w2) .eq. 'CONSTANT')then
                     backspace 700
                     read(700,*,err=11)keyword, w1, w2, r
                     tmp_sgn(size(signal)+1)%r = r
                  end if
                  call move_alloc(tmp_sgn, signal)
               end if
            end if

         end if

   !        first symbol of the next line
         read(700,'(a1)',err=11, end=11)c
         numlin = numlin + 1
     end do
   
   close(700)

   return

!  input error
11 write(*,*)'****ERROR: input error at line: ', numlin
   stop
end subroutine INP_read_input


!----------------------------------------------------------------------------------------------
! Write material macroscopic XS to file
   ! Create a folder containing
   ! material_list.txt
   ! Material directories: <mat0>, <mat1>, ...
   ! In each directory, a collection of csv files, containing the following 
   !  - total.csv
   !  - absorbtion.csv
   !  - fission.csv
   !  - nu_fission.csv
   !  - chi.csv
   !  - scattering_p<0>.csv 
!----------------------------------------------------------------------------------------------




!----------------------------------------------------------------------------------------------
! Given vectors x and y, returns interpolated value yq corresponding to xq
!----------------------------------------------------------------------------------------------
function MATH_interp1(x, y, xq) result(yq)

   real(real64), intent(in) :: x(:), y(:), xq
   real(real64) :: yq, dx
   integer :: i, n

   n = size(x)

      if((xq - x(1)) * (xq - x(n)) .le. 0.0)then
!     inside bounds, proceed
      do i = 1, n - 1
         dx = x(i+1) - x(i)
         if((xq - x(i)) * (xq - x(i+1)) .le. 0.0 .and. dx .ne. 0.0)then
            yq = y(i) + (y(i+1) - y(i)) * (xq - x(i)) / dx
            return
         end if
      end do
   end if

   !  out-of-bounds
   if((xq .lt. x(1) .and. x(1) .lt. x(n)) .or. (xq .gt. x(1) .and. x(1) .gt. x(n)))then
      yq = y(1)
   else
      yq = y(n)
   end if

end function MATH_interp1

!----------------------------------------------------------------------------------------------
! Evaluates signals
!----------------------------------------------------------------------------------------------
subroutine SGN_evaluate_signals(signal)

   type(signal_template), intent(inout) :: signal

   if(signal%sgntype .eq. 'CONSTANT')signal%v = signal%r

end subroutine SGN_evaluate_signals

!----------------------------------------------------------------------------------------------
! Calculates a 1D array of sigma-zeros for each isotope of the material mat and then calculate 
! all macroscopic cross sections
!----------------------------------------------------------------------------------------------
subroutine XS_calculate_macro(mat, iso, core, signal, printmicro)

   logical, intent(in) :: printmicro ! whether to print microscopic XS to file 
   integer :: ig, i, j, k, ng, niso, iter, ii, f, t
   real(real64) :: error, s, volfrac
   real(real64), allocatable :: sig0(:,:)
   character*100, allocatable :: sgnname(:)
   type(mat_template), intent(inout):: mat
   type(core_template), intent(in) :: core
   type(iso_template), intent(in) :: iso(:)
   type(signal_template), intent(in) :: signal(:)
   type(iso_template), allocatable :: iso_t(:), iso_s(:)
   
   ! Ym note: my additions
   integer :: nl, inl, gfrom, gto
   character(len=10) :: sdummy

   integer, dimension(:), allocatable :: dims


   write(*,*)'Processing material: ', mat%matname

   !  number of energy groups
   ng = iso(1)%ng
   !  number of isotopes in material mat
   niso = size(mat%isoname)

   !  if there is a void, correct number densities
   if(mat%d_void .gt. 0.0)then
   !     volume fraction of void
      volfrac = 0.785398*mat%d_void**2 / (0.866025 * core%pitch**2)
      do i = 1, niso
         mat%numden(ii) = mat%numden(ii)/(1.0 - volfrac)
      end do
   end if

   allocate(iso_t(niso))
   allocate(iso_s(niso))
   allocate(sig0(niso,ng))

   !  allocate space for sig0 for each isotope of material mat and initialize
   allocate(mat%sig0(niso, ng))
   do ig = 1, ng
      do i = 1, niso
         mat%sig0(i,ig) = 1.0
         sig0(i,ig) = 1.0
      end do
   end do

   !  list of signals
   allocate(sgnname(size(signal)))
   do i = 1, size(signal)
      sgnname(i) = signal(i)%sgnname
   end do

   !  make interpolation for temperature
   do i = 1, niso
   !     index of isotope in the core%isoname
      j = findloc(core%isoname,mat%isoname(i), dim=1)
   !     index of temperature signal
      k = findloc(sgnname,mat%tsignal(i), dim=1)
      call XS_interpolate_temp(iso(j), signal(k)%v, iso_t(i))
   end do

   !  if material consists of only one isotope then keep sig0 = 1
   if (niso .eq. 1) then

      call XS_interpolate_sig0(iso_t(1), sig0(1,:), iso_s(1))

   else
      !  error of sig0 calculation
      error = 1.0
      !  iteration loop
      iter = 0   
      do while(error .gt. 1.0e-4)
         do i = 1, niso
            call XS_interpolate_sig0(iso_t(i), sig0(i,:), iso_s(i))
         end do
         error = 0.0
         do ig = 1, ng
            sig0(:,ig) = 0.0
            do i = 1, niso
               do ii = 1, niso
                  if(ii .ne. i) sig0(i,ig) = sig0(i,ig) + mat%numden(ii)*iso_s(ii)%sigt(ig,1,1)
               end do
               sig0(i,ig) = sig0(i,ig)/mat%numden(i)
               sig0(i,ig) = min(sig0(i,ig), 1.0e10)
               sig0(i,ig) = max(sig0(i,ig), 1.0)
         !           error control
               error = max(error, abs(mat%sig0(i,ig) - sig0(i,ig)))
               mat%sig0(i,ig) = sig0(i,ig)
            end do
         end do
         iter = iter + 1
         if(iter .gt. 1000)then
            write(*,*)'****ERROR: too many sigma-zero iterations'
            stop
         end if
      end do
   end if

   !  macroscopic total cross sections
   allocate(mat%SigT(ng))
   do ig = 1, ng
      mat%SigT(ig) = 0.0
      do i = 1, niso
         mat%SigT(ig) = mat%SigT(ig) + mat%numden(i)*iso_s(i)%sigt(ig,1,1)
      end do
   end do

   !  macroscopic fission cross sections
   allocate(mat%SigF(ng))
   do ig = 1, ng
      mat%SigF(ig) = 0.0
      do i = 1, niso
         mat%SigF(ig) = mat%SigF(ig) + mat%numden(i)*iso_s(i)%sigf(ig,1,1)
      end do
   end do

   !  macroscopic production cross sections
   allocate(mat%SigP(ng))
   do ig = 1, ng
      mat%SigP(ig) = 0.0
      do i = 1, niso
         mat%SigP(ig) = mat%SigP(ig) + mat%numden(i) * iso_s(i)%nubar(ig) * iso_s(i)%sigf(ig,1,1)
      end do
   end do

   !  fission spectrum for material
   allocate(mat%chi(ng))
   do ig = 1, ng
      mat%chi(ig) = 0.0
      do i = 1, niso
         if (iso_s(i)%sigf(ig,1,1) .gt. 0) then 
            mat%chi(ig) = mat%chi(ig) + mat%numden(i) * iso_s(i)%chi(ig) !* iso_s(i)%nubar(ig) * iso_s(i)%sigf(ig,1,1)
         end if
      end do
   end do

   s = sum(mat%chi)
   if(s .gt. 0.0)mat%chi = mat%chi/s

   !  macroscopic kerma factors
   allocate(mat%Kerma(ng))
   do ig = 1, ng
      mat%Kerma(ig) = 0.0
      do i = 1, niso
         mat%Kerma(ig) = mat%Kerma(ig) + mat%numden(i)*iso_s(i)%kerma(ig,1,1)
      end do
   end do


   !  find the maximum Legendre order among all isotopes
   ii = 0
   do i = 1, niso
      ii = max(ii, iso_s(i)%nl)
   end do

   !  allocate space for scattering matrix
   allocate(mat%SigS(ii,ng,ng))
   mat%SigS = 0.0
   do i = 1, niso
      do ii = 1, iso_s(i)%nl
         do j = 1, iso_s(i)%nonze
            f = iso_s(i)%f_e(j)
            t = iso_s(i)%t_e(j)
            mat%SigS(ii,f,t) = mat%SigS(ii,f,t) + mat%numden(i) * iso_s(i)%sige(ii,1,j,1)
         end do
         do j = 1, iso_s(i)%nonzi
            f = iso_s(i)%f_i(j)
            t = iso_s(i)%t_i(j)
            mat%SigS(ii,f,t) = mat%SigS(ii,f,t) + mat%numden(i) * iso_s(i)%sigi(ii,j,1)
         end do
      end do
   end do

   allocate(mat%SigST(ng))
   do ig = 1, ng
      mat%SigST(ig) = sum(mat%SigS(1,ig,:))
   end do

   allocate(mat%SigA(ng))
   do ig = 1, ng
      mat%SigA(ig) = mat%SigT(ig) - mat%SigST(ig)
   end do

   ! YM note: This is where I write the macroscopic properties to file 
   ! I make modifications to write the temperature corrected, 
   ! background XS corrected microscopic XS to file  
   ! create all output directories here 
   call system("mkdir -p material_xs/" // mat%matname // "/macro")

   if (printmicro) then
      do i = 1, niso
         call system("mkdir -p material_xs/" // mat%matname // "/" // iso_s(i)%isoname)
      end do
   end if 

   ! 1000(1pe12.5) -> do at most 1000 times, 1p-> x10, e12.5, scientific notation 12 char with 5 after decimal ,: (colon) terminates the format if there are no more items, thus stopping a trailing comma being written.  

   ! Total cross section macro then micro
   
   open(1000, file="material_xs/" // mat%matname // "/macro/total.csv")
   write(1000,*)'GroupIndex,Value'
   call write_csv_1D(1000, mat%SigT)
   close(1000)
   ! note: this writes out the temperature and sig0 interpolated total XS. It is also possible to write out the values of all temperature and sig0s, The information would be found in iso_t%sigt
   if (printmicro) then
      do i = 1, niso
         open(1000, file="material_xs/" // mat%matname // "/"// iso_s(i)%isoname // "/total.csv")
         write(1000,*)'GroupIndex,Value'
         call write_csv_1D(1000,iso_s(i)%sigt(:,1,1))
         close(1000)
      end do
   end if

   write(*,*)'    total'
   
   ! Get group fluxes for each isotope (as function of temperature and sig0s)
   if (printmicro) then
      do i = 1, niso
         open(1000, file="material_xs/" // mat%matname // "/"// iso_s(i)%isoname // "/scalarflux.csv")  
         write(1000,*)'GroupIndex,TempIndex,Sig0Index,Value'
         call write_csv_3D(1000,iso_s(i)%flux)
         ! call write_csv_nd_with_index(1000, reshape(iso_s(i)%flux, [size(iso_s(i)%flux)]) , shape(iso_s(i)%flux)) ! unit of file, flattened array, shape of original array
         close(1000)
      end do 
   end if 
   write(*,*)'    flux'

   ! Fission cross-section
   open(1000, file="material_xs/" // mat%matname // "/macro/fission.csv")
   write(1000,*)'GroupIndex,Value'
   call write_csv_1D(1000,mat%SigF)
   ! write(1000 , '(1000(1pe14.7,:,","))') (mat%SigF(ig), ig=1,ng) 
   close(1000)
   if (printmicro) then
      do i = 1, niso
         open(1000, file="material_xs/" // mat%matname // "/"// iso_s(i)%isoname // "/fission.csv")
         write(1000,*)'GroupIndex,Value'
         call write_csv_1D(1000,iso_s(i)%sigf(:, 1, 1))
         ! write(1000, '(1000(1pe14.7,:,","))') (iso_s(i)%sigf(ig, 1, 1), ig=1,ng) 
         close(1000)
      end do
   end if
   write(*,*)'    fiss'

   ! Kerma factor
   open(1000, file="material_xs/" // mat%matname // "/macro/kerma.csv")
   write(1000,*)'GroupIndex,Value'
   call write_csv_1D(1000,mat%Kerma)
   ! write(1000 , '(1000(1pe14.7,:,","))') (mat%Kerma(ig), ig=1,ng) 
   close(1000)
   if (printmicro) then
      do i = 1, niso
         open(1000, file="material_xs/" // mat%matname // "/"// iso_s(i)%isoname // "/kerma.csv")
         write(1000,*)'GroupIndex,Value'
         call write_csv_1D(1000,iso_s(i)%kerma(:, 1, 1))
         !write(1000, '(1000(1pe14.7,:,","))') (iso_s(i)%kerma(ig, 1, 1), ig=1,ng) 
         close(1000)
      end do
   end if 
   write(*,*)'    kerma'

   ! Nu-Fission cross-section
   open(1000, file="material_xs/" // mat%matname // "/macro/nu_fission.csv")
   write(1000,*)'GroupIndex,Value'
   call write_csv_1D(1000,mat%SigP)
   ! write(1000, '(1000(1pe14.7,:,","))') (mat%SigP(ig), ig=1,ng)       
   close(1000)
   if (printmicro) then
      do i = 1, niso
         open(1000, file="material_xs/" // mat%matname // "/"// iso_s(i)%isoname // "/nu.csv")
         write(1000,*)'GroupIndex,Value'
         call write_csv_1D(1000,iso_s(i)%nubar)
         ! write(1000, '(1000(1pe14.7,:,","))') (iso_s(i)%nubar(ig), ig=1,ng) 
         close(1000)
      end do
   end if 
   write(*,*)'    nu-fiss'

   ! Chi 
   open(1000 , file="material_xs/" // mat%matname // "/macro/chi.csv")
   write(1000,*)'GroupIndex,Value'
   call write_csv_1D(1000,mat%chi)
   ! write(1000, '(1000(1pe14.7,:,","))') (mat%chi(ig), ig=1,ng)       
   close(1000)
   if (printmicro) then
      do i = 1, niso
         open(1000, file="material_xs/" // mat%matname // "/"// iso_s(i)%isoname // "/chi.csv")
         write(1000,*)'GroupIndex,Value'
         call write_csv_1D(1000,iso_s(i)%chi)
         !write(1000, '(1000(1pe14.7,:,","))') (iso_s(i)%chi(ig), ig=1,ng) 
         close(1000)
      end do
   end if
   write(*,*)'    chi'

   ! Absorption
   open(1000 , file="material_xs/" // mat%matname // "/macro/absorption.csv")
   write(1000,*)'GroupIndex,Value'
   ! write(1000, '(1000(1pe14.7,:,","))') (mat%SigA(ig), ig=1,ng)       
   call write_csv_1D(1000,mat%SigA)
   close(1000)

   ! scattering. Convert to sparse representation, csv with the fields, LEGENDRE, FROM, TO, VALUE
   open(1000, file="material_xs/" // mat%matname // "/macro/scatter.csv")
   write(1000, *) "LegendreIndex,FromIndex,ToIndex,Value"
   nl = size(mat%SigS,1)
   do inl=1, nl
      do gfrom=1, ng
         do gto=1, ng
            s = mat%SigS(inl,gfrom,gto)
            ! only write to file if more than a certain threshold
            if (abs(s) > 1.0e-15) then
               write(1000,'(I0, ",", I0, ",", I0, ",", 1pe14.7)') inl, gfrom, gto, s
            end if
         end do ! gto      
      end do ! gfrom
   end do ! end nl 
   close(1000)
   write(*,*)'    scatter'

   ! save also isotopes and temperature of the material 
   open(1000, file="material_xs/" // mat%matname // "/isolist.csv")
   write(1000, *)"Iso,Temp,Conc"
   do i = 1, niso
      k = findloc(sgnname,mat%tsignal(i), dim=1)
      write(1000,'(A, ",",  1pe14.7, ",", 1pe14.7)')trim(mat%isoname(i)), signal(k)%v, mat%numden(i)
   end do
   close(1000)
   write(*,*)'    isolist'

   ! save sig0 values (each row is different nuclide)
   open(1000, file="material_xs/" // mat%matname // "/sig0.csv")
   write(1000,*)"Iso,GroupIndex"
   call write_csv_2D(1000,mat%sig0)
   ! do i = 1, niso
   !     write(1000, '(1000(1pe14.7,:,","))') (mat%sig0(i, ig), ig=1,ng) 
   ! end do
   close(1000)
   write(*,*)'    sig0'

   ! isotope-wise scattering (split between elastic and inelastic)
   if (printmicro) then
      do i = 1, niso
         ! elastic component 
         
         open(1000, file="material_xs/" // trim(mat%matname) // "/" // trim(iso_s(i)%isoname) // "/scatter_elas")
         write(1000,*)'LegendreIndex,From,To,Value'

         open(1005, file="material_xs/" // trim(mat%matname) // "/" // trim(iso_s(i)%isoname) // "/scatter_inelas")
         write(1005,*)'LegendreIndex,FromIndex,ToIndex,Value'

         do ii = 1, iso_s(i)%nl
            do j = 1, iso_s(i)%nonze
               f = iso_s(i)%f_e(j)
               t = iso_s(i)%t_e(j)
               write(1000,'(I0, ",", I0, ",", I0, ",", 1pe14.7)')ii, f, t, iso_s(i)%sige(ii,1,j,1)
            end do ! end j

            do j = 1, iso_s(i)%nonzi
               f = iso_s(i)%f_i(j)
               t = iso_s(i)%t_i(j)
               write(1005,'(I0, ",", I0, ",", I0, ",", 1pe14.7)')ii, f, t, iso_s(i)%sigi(ii,j,1)
            end do ! end j         

         end do ! end ii
         
         close(1000)
         close(1005)
      end do
      write(*,*)'    scatter_iso'
   end if 
   

   write(*,*)'completed material: ', mat%matname
   


   deallocate(iso_t)
   deallocate(iso_s)
   deallocate(sig0)

end subroutine XS_calculate_macro

!----------------------------------------------------------------------------------------------
! The subroutine searches matrix cards for cross sections sig from file mf for
! reaction mt, temperature index itemp, legendre order nl and returns outp(ng,nsig0), where ng
! is the number of energy groups and nsig0 is the number of sigma-zeros.
!----------------------------------------------------------------------------------------------
subroutine XS_extract_mf_mt(iso, mf, mt, itemp, nl, ncards, cards, outp)
   ! nl : requested Legendre component
   integer, intent(in) :: mf, mt, itemp, nl, ncards
   integer :: ntemp, irow, nlgn, nsig0, ig, k, irownew
   real(real64), intent(in) :: cards(:,:)
   real(real64), intent(out) :: outp(:,:)
   real(real64), allocatable :: r(:)
   type(iso_template), intent(in) :: iso

   !  find index irow of the row with required mf and mt
   ntemp = 0
   irow = 1
   do while(ntemp .lt. itemp .and. irow .lt. ncards)
      if(cards(irow,8) .eq. mf .and. cards(irow,9) .eq. mt .and. cards(irow,10) .eq. 1) ntemp = ntemp + 1
      irow = irow + 1
   end do
   !  if required mf and mt NOT found
   if(ntemp .eq. 0)then
      outp = 0.0
      return
   end if

   !  number of Legendre components
   nlgn = cards(irow-1,3)
   !  number of sigma-zeros
   nsig0 = cards(irow-1,4)
   if(nl .gt. nlgn)then
      write(*,*)'****ERROR: legendre order requested in XS_extract_mf_mt is higher than available in the library.'
      stop
   end if
   irow = irow + 1
   allocate(r(nsig0*nlgn*2))
   if(mf .eq. 3)then
      do while(cards(irow,8) .eq. mf .and. cards(irow,9) .eq. mt)
         ig = int(cards(irow-1,6))
         call XS_extract_n_words(nsig0*nlgn*2, irow, cards, irownew, r)
   !        the first nlgn*nsig0 words are flux -- skip.
         r = r(nsig0*nlgn+1:)
         do k = 1, nsig0
            outp(ig,k) = r((k - 1) * nlgn + nl) !check
         end do
         irow = irownew + 2
      end do
   else if(mf .eq. 5)then
      call XS_extract_n_words(iso%ng, irow, cards, irownew, r)
      outp(:,1) = r
   end if
   deallocate(r)

end subroutine XS_extract_mf_mt


!----------------------------------------------------------------------------------------------
! YM NJOY provides the group averaged spectra. This might be useful too
!----------------------------------------------------------------------------------------------
subroutine read_scalar_flux(itemp, ncards, cards, tflux)
   ! this returns the energy group averaged flux value (see eqn 253 of njoy) by (energy group, sigma0). Note that we are only taking the scalar flux
   integer, intent(in) :: ncards, itemp
   real(real64), intent(in) :: cards(:,:)
   real(real64), intent(out) :: tflux(:,:) ! (ng,nsig0)

   integer :: ntemp, mf, mt, ig, nlgn, nl, nsig0, irow, irownew, k
   real(real64), allocatable :: r(:)

   ! Use mf=3, mt=1 to find the flux
   mf = 3
   mt = 1
   nl = 1 ! only get scalar flux

   ntemp = 0
   irow = 1
   do while(ntemp .lt. itemp .and. irow .lt. ncards)
      if(cards(irow,8) .eq. mf .and. cards(irow,9) .eq. mt .and. cards(irow,10) .eq. 1) ntemp = ntemp + 1
      irow = irow + 1
   end do


   !  number of Legendre components
   nlgn = cards(irow-1,3)
   !  number of sigma-zeros
   nsig0 = cards(irow-1,4)

   allocate(r(nsig0*nlgn*2))

   irow = irow + 1


   do while(cards(irow,8) .eq. mf .and. cards(irow,9) .eq. mt)
      ig = int(cards(irow-1,6))
      call XS_extract_n_words(nsig0*nlgn*2, irow, cards, irownew, r)
      !        the first nlgn*nsig0 words are flux 
      do k = 1, nsig0
         tflux(ig,k) = r((k - 1) * nlgn + nl) !check
      end do
           
      irow = irownew + 2
   end do

   deallocate(r)

end subroutine read_scalar_flux







!----------------------------------------------------------------------------------------------
! Extracts Legendre-expanded scattering cross sections (mf=6) from the cards matrix
! for a given reaction mt, temperature index itemp and Legendre component nl.
! Returns 1D integer arrays f and t for from- and to-group numbers and 1D real(real64) array r for corresponding 
! cross sections
!----------------------------------------------------------------------------------------------
subroutine XS_extract_mf6(mt, itemp, nl, ncards, cards, f, t, sig)

   integer, intent(in) :: mt, itemp, nl, ncards
   integer, intent(out) :: f(:), t(:)
   integer :: irow, ntemp, nlgn, nsig0, ng2, ig2lo, nw, ig, ito, ilgn, isig0, k, l, nonz, irownew
   real(real64), intent(in) :: cards(:,:)
   real(real64), intent(out) :: sig(:,:)
   real(real64), allocatable :: r(:)

   irow = 1
   ntemp = 0
   do while(irow .lt. ncards)
   !     find the row with mf=6 & mt
      if(int(cards(irow,8)) .eq. 6 .and. int(cards(irow,9)) .eq. mt)then
   !        if this is the first line of mf=6 & mt (new block header), initialize
         if(int(cards(irow,10)) .eq. 1)then
   !           number of nonzeros
            nonz = 0
   !           number of Legendre components
            nlgn = int(cards(irow,3))
   !           number of sigma-zeros
            nsig0 = int(cards(irow,4))
   !           temperature index
            ntemp = ntemp + 1
   !           row index
            irow = irow + 1
         end if
      
   !        number of secondary positions
         ng2 = int(cards(irow,3))
   !        index to lowest nonzero group
         ig2lo = int(cards(irow,4))
   !        number of words to be read
         nw = int(cards(irow,5))
   !        current group index
         ig = int(cards(irow,6))

   !        row index
         irow = irow + 1

   !        extract nw words from cards
         call XS_extract_n_words(nw, irow, cards, irownew, r)
         irow = irownew

         if(ntemp .eq. itemp)then
   !           the first nlgn*nsig0 words are flux -- skip.
            r = r(nlgn*nsig0+1:)
            k = 0
            do ito = ig2lo, ig2lo + ng2 - 2
               nonz = nonz + 1
               f(nonz) = ig
               t(nonz) = ito
               l = 0
               do isig0 = 1, nsig0
               do ilgn = 1, nlgn
                  k = k + 1
                  if(ilgn .eq. nl)then
                     l = l + 1
                     sig(nonz,l) = r(k)
                  end if
               end do
               end do
            end do
         end if
      end if
      irow = irow + 1
   end do

end subroutine XS_extract_mf6

!----------------------------------------------------------------------------------------------
! Reads n words from row irow of matrix cards and returns them in
! vector r together with the new row number irownew, i.e. the row where the last word was read.
!----------------------------------------------------------------------------------------------
subroutine XS_extract_n_words(n, irow, cards, irownew, r)

   integer, intent(in) :: n, irow
   integer, intent(out) :: irownew
   integer :: ii, jj, nremain, k
   real(real64), intent(in) :: cards(:,:)
   real(real64), intent(out), allocatable :: r(:)

   allocate(r(n))
   irownew = irow

   !  read full lines (6 values each)
   k = 1
   do ii = 1, int(n/6)
      do jj = 1, 6
         r(k) = real(cards(irownew, jj))
         k = k + 1
      end do
      irownew = irownew + 1
   end do

   !  handle last partial line if needed
   nremain = mod(n, 6)
   if (nremain == 0) then
      irownew = irownew - 1
   else
      do jj = 1, nremain
         r(k) = real(cards(irownew, jj))
         k = k + 1
      end do
   end if

end subroutine XS_extract_n_words

!----------------------------------------------------------------------------------------------
! Returns number of Legendre components for scattering cross sections (mf=6) from the cards matrix
! for a given reaction mt.
!----------------------------------------------------------------------------------------------
function XS_extract_nlgn_mf6(mt, ncards, cards) result(nlgn)

   integer, intent(in) :: mt, ncards
   integer :: nlgn, irow, ntemp
   real(real64), intent(in) :: cards(:,:)

   irow = 1
   nlgn = 0
   do while(irow .le. ncards)
   !     find the first row with mf=6 & mt (new block header)
      if(int(cards(irow,8)) .eq. 6 .and. int(cards(irow,9)) .eq. mt .and. int(cards(irow,10)) .eq. 1)then
   !        number of Legendre components
         nlgn = int(cards(irow,3))
         exit
      end if
      irow = irow + 1
   end do

end function XS_extract_nlgn_mf6

!----------------------------------------------------------------------------------------------
! Returns number of nonzeros for Legendre-expanded scattering cross sections (mf=6) from the cards matrix
! for a given reaction mt.
!----------------------------------------------------------------------------------------------
function XS_extract_nonz_mf6(mt, ncards, cards) result(nonz)

   integer, intent(in) :: mt, ncards
   integer :: irow, ntemp, nonz, ig2lo, irownew, ito, ng2, nw
   real(real64), intent(in) :: cards(:,:)
   real(real64), allocatable :: r(:)

   irow = 1
   ntemp = 0
   nonz = 0
   do while(irow .le. ncards)
   !     find the row with mf=6 & mt
      if(int(cards(irow,8)) .eq. 6 .and. int(cards(irow,9)) .eq. mt)then
   !        if this is the first line of mf=6 & mt (new block header), initialize
         if(int(cards(irow,10)) .eq. 1)then
   !           number of nonzeros
            nonz = 0
   !           temperature index
            ntemp = ntemp + 1
   !           row index
            irow = irow + 1
         end if
         
   !        number of secondary positions
         ng2 = int(cards(irow,3))
   !        index to lowest nonzero group
         ig2lo = int(cards(irow,4))
   !        number of words to be read
         nw = int(cards(irow,5))
      
   !        row index
         irow = irow + 1
      
   !        extract nw words from cards
         call XS_extract_n_words(nw, irow, cards, irownew, r)
         irow = irownew
      
         do ito = ig2lo, ig2lo + ng2 - 2
            nonz = nonz + 1
         end do
      end if
      irow = irow + 1
   end do

end function XS_extract_nonz_mf6


!----------------------------------------------------------------------------------------------
! Check if reaction is in mf=6 or represented by mf=4 and mf=5
! for a given reaction mt.
!----------------------------------------------------------------------------------------------
! function XS_in_mf6

! end function XS_in_mf6


!----------------------------------------------------------------------------------------------
! Given isotope iso_t and background cross section sig0, returns iso_s in which all cross sections are 
! interpolated for sig0 (dimension of nsig0 becomes 1)
!----------------------------------------------------------------------------------------------
subroutine XS_interpolate_sig0(iso_t, sig0, iso_s)

   integer i, j, k
   real(real64), intent(in) :: sig0(:)
   type(iso_template), intent(in):: iso_t
   type(iso_template), intent(out):: iso_s

   !  make a copy of instance iso_t
   iso_s = iso_t

   !  total cross section
   deallocate(iso_s%sigt)
   allocate(iso_s%sigt(iso_t%ng,1,1))
   do i = 1, iso_t%ng
      iso_s%sigt(i,1,1) = MATH_interp1(log10(iso_t%sig0), iso_t%sigt(i,1,:), log10(sig0(i)))
   end do

   !  fission cross section
   deallocate(iso_s%sigf)
   allocate(iso_s%sigf(iso_t%ng,1,1))
   do i = 1, iso_t%ng
      iso_s%sigf(i,1,1) = MATH_interp1(log10(iso_t%sig0), iso_t%sigf(i,1,:), log10(sig0(i)))
   end do

   !  elastic scattering cross section
   deallocate(iso_s%sige)
   allocate(iso_s%sige(iso_t%nl,1,iso_t%nonze,1))
   do k = 1, iso_t%nl
   do j = 1, iso_t%nonze
      i = iso_t%f_e(j)
      iso_s%sige(k,1,j,1) = MATH_interp1(log10(iso_t%sig0), iso_t%sige(k,1,j,:), log10(sig0(i)))
   end do
   end do

   !  kerma factor
   deallocate(iso_s%kerma)
   allocate(iso_s%kerma(iso_t%ng,1,1))
   do i = 1, iso_t%ng
      iso_s%kerma(i,1,1) = MATH_interp1(log10(iso_t%sig0), iso_t%kerma(i,1,:), log10(sig0(i)))
      !iso_s%kerma(i,1,1) = MATH_interp1(iso_t%sig0, iso_t%kerma(i,1,:), sig0(i))
   end do

end subroutine XS_interpolate_sig0

!----------------------------------------------------------------------------------------------
! Given isotope iso and temperature temp, returns iso_t in which all cross sections are 
! interpolated for temperature temp (dimension of ntemp becomes 1)
!----------------------------------------------------------------------------------------------
subroutine XS_interpolate_temp(iso, temp, iso_t)

   integer i, j, k
   real(real64), intent(in) :: temp
   type(iso_template), intent(in):: iso
   type(iso_template), intent(out):: iso_t

   !  make a copy of instance iso
   iso_t = iso

   !  total cross section
   deallocate(iso_t%sigt)
   allocate(iso_t%sigt(iso%ng,1,iso%nsig0))
   do i = 1, iso%ng
      do j = 1, iso%nsig0
         iso_t%sigt(i,1,j) = MATH_interp1(iso%temp, iso%sigt(i,:,j), temp)
      end do
   end do

   !  fission cross section
   deallocate(iso_t%sigf)
   allocate(iso_t%sigf(iso%ng,1,iso%nsig0))
   do i = 1, iso%ng
      do j = 1, iso%nsig0
         iso_t%sigf(i,1,j) = MATH_interp1(iso%temp, iso%sigf(i,:,j), temp)
      end do
   end do

   !  elastic scattering cross section
   deallocate(iso_t%sige)
   allocate(iso_t%sige(iso%nl,1,iso%nonze,iso%nsig0))
   do i = 1, iso%nonze
      do j = 1, iso%nsig0
         do k = 1, iso%nl
            iso_t%sige(k,1,i,j) = MATH_interp1(iso%temp, iso%sige(k,:,i,j), temp)
         end do
      end do
   end do

end subroutine XS_interpolate_temp

!--------------------------------------------------------------------------------------------------
! XS: Read GENDF file for isotope iso
!--------------------------------------------------------------------------------------------------
subroutine XS_read_gendf(iso, core)

   integer :: i, j, k, ncards, itemp, nlgn, nl, mt
   integer, allocatable :: f_i(:), t_i(:)
   real(real64), allocatable :: cards(:,:), outp(:,:), r(:)
   character(:), allocatable :: isoname, fname
   character(len=12) :: w
   character (len=80) :: line = ''
   type(core_template) :: core
   type(iso_template), intent(inout) :: iso

   fname = core%path2data // '/' // iso%isoname
   write(*,*)'Processing GENDF file: ', fname

   !  READ GENDF INTO CARDS ARRAY

   !  count number of cards
   open(700, file=fname)
      read(700,'(a)')line
      ncards = 1
      do while(line(69:70) .ne. '-1')
         read(700,'(a)')line
         ncards = ncards + 1
      end do
      ncards = ncards - 1
   close(700)

   !  allocate
   allocate(cards(ncards,10))

   open(700, file=fname)
      read(700,'(a)')line
      do k = 1, ncards
         read(700,'(a)')line
   !        six words of 11-character length
         do i = 1, 6
            w = line((i-1)*11+1:i*11)
   !           if the w is float
            if(index(w,'.') .ne. 0)then
   !              start from 2 to skip leading minus
               do j=2,11                              
                  if(w(j:j) .eq. '+' .or. w(j:j) .eq. '-')then
   !                    insert missing E
                     w = w(:j-1) // 'E' // w(j:)
   !                    convert string to float
                     read(w,*)cards(k,i)
                     exit
                  end if
               end do
            else if(w .ne. '')then
   !              convert string to float, although it is integer
               read(w,*)j
               cards(k,i) = real(j)
            end if
         end do
!        append three last values (mf, mt and line number)
         read(line(71:72),*)j
         cards(k,8) = real(j)
         read(line(73:75),*)j
         cards(k,9) = real(j)
         read(line(76:80),*)j
         cards(k,10) = real(j)
      end do
   close(700)

   !  PROCESS CARDS

   !  number of energy groups (check versus input ng)
   iso%ng = int(cards(2,3))

   !  find number of base temperatures and values of base temperatures using mf=1 and mt=451
   iso%ntemp = 0
   do i = 1, ncards
      if(int(cards(i,8)) .eq. 1 .and. int(cards(i,9)) .eq. 451 .and. int(cards(i,10)) .eq. 1)then
         iso%ntemp = iso%ntemp + 1
      end if
   end do
   allocate(iso%temp(iso%ntemp))
   k = 0
   do i = 1, ncards
      if(int(cards(i,8)) .eq. 1 .and. int(cards(i,9)) .eq. 451 .and. int(cards(i,10)) .eq. 1)then
         k = k + 1
         iso%temp(k) = cards(i+1,1)
      end if
   end do

   !  find number of sigma-zeros and sigma-sero values in xs parametrization
   iso%nsig0 = cards(1,4)
   allocate(r(iso%nsig0+1))
   call XS_extract_n_words(iso%nsig0+1, 3, cards, i, r)
   iso%sig0 = r(2:)
   deallocate(r)

   allocate(outp(iso%ng,1))

   !  fission spectrum (mf = 5, mt = 18)
   allocate(iso%chi(iso%ng))
   call XS_extract_mf_mt(iso, 5, 18, 1, 1, ncards, cards, outp)
   iso%chi = outp(:,1)

   !  nubar (mf = 3, mt = 452)
   allocate(iso%nubar(iso%ng))
   call XS_extract_mf_mt(iso, 3, 452, 1, 1, ncards, cards, outp)
   iso%nubar = outp(:,1)

   deallocate(outp)

   !  kerma factors (mf = 3, mt = 301)
   allocate(iso%kerma(iso%ng,1,iso%nsig0))
   call XS_extract_mf_mt(iso, 3, 301, 1, 1, ncards, cards, iso%kerma(:,1,:))

   allocate(iso%sigf(iso%ng,iso%ntemp,iso%nsig0))
   do itemp = 1, iso%ntemp
   !     fission xs (mf = 3, mt = 18)
      call XS_extract_mf_mt(iso, 3, 18, itemp, 1, ncards, cards, iso%sigf(:,itemp,:))
   end do

   allocate(iso%sigt(iso%ng,iso%ntemp,iso%nsig0))
   allocate(iso%flux(iso%ng,iso%ntemp,iso%nsig0))
   ! get also flux 

   do itemp = 1, iso%ntemp
      !     total cross section (mf=3, mt=1)
      call XS_extract_mf_mt(iso, 3, 1, itemp, 1, ncards, cards, iso%sigt(:,itemp,:))
      !     read scalar flux (legendre index = 1)
      call read_scalar_flux(itemp, ncards, cards, iso%flux(:,itemp,:))

   end do

   !  count number of non-zeros in elastic matrix and allocate space
   iso%nonze = XS_extract_nonz_mf6(2, ncards, cards)
   iso%nl = XS_extract_nlgn_mf6(2, ncards, cards)
   allocate(iso%sige(iso%nl,iso%ntemp,iso%nonze,iso%nsig0))
   iso%sige = 0.0
   allocate(iso%f_e(iso%nonze))
   allocate(iso%t_e(iso%nonze))
   do itemp = 1, iso%ntemp
      do nlgn = 1, iso%nl
         call XS_extract_mf6(2, itemp, nlgn, ncards, cards, iso%f_e, iso%t_e, iso%sige(nlgn,itemp,:,:))
      end do
   end do

   ! Count number of non-zeros in inelastic matrix and allocate space. 
   ! YM: Note: In new versions of JEFF/ENDF-B libraries, for some isotopes, inelastic rxn MT={51,..91} is not given by MF=6 but rather MF=4 and MF=5
   
   iso%nonzi = 0
   do mt = 51, 91
      ! Routine to check if mt is in mf6 or in mf4/5 
      iso%nonzi = iso%nonzi + XS_extract_nonz_mf6(mt, ncards, cards)
   end do
   !  inelastic scattering (mt = 51... 91)
   allocate(iso%sigi(iso%nl,iso%nonzi,iso%nsig0))
   iso%sigi = 0.0
   allocate(iso%f_i(iso%nonzi))
   allocate(iso%t_i(iso%nonzi))

   iso%nonzi = 0
   do mt = 51, 91
      !     count number of non-zeros in inelastic matrix and allocate space
      ! YM note: what happens if inelastic is in mf4/5? 
      k = XS_extract_nonz_mf6(mt, ncards, cards)
      nl = XS_extract_nlgn_mf6(mt, ncards, cards)
      allocate(outp(k,iso%nsig0))
      allocate(f_i(k))
      allocate(t_i(k))
      call XS_extract_mf6(mt, 1, 1, ncards, cards, f_i, t_i, outp)
      do i = 1, k
         iso%f_i(iso%nonzi+i) = f_i(i)
         iso%t_i(iso%nonzi+i) = t_i(i)
         do j = 1,iso%nsig0
            iso%sigi(1,iso%nonzi+i,j) = outp(i,j)
         end do
      end do
      do nlgn = 2, nl
         do i = 1, k
            do j = 1,iso%nsig0
               iso%sigi(nlgn,iso%nonzi+i,j) = outp(i,j)
            end do
         end do
      end do
      iso%nonzi = iso%nonzi + k
      deallocate(outp)
      deallocate(f_i)
      deallocate(t_i)
   end do
   deallocate(cards)

end subroutine XS_read_gendf

end module MONTE_extractmacro
 