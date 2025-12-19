program main
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
   use globals
   use MONTE_extractmacro
   

   implicit none
   integer :: i, j, count_start, count_end, count_rate
   integer :: dt(8)
   real :: r
   character(20) :: timestamp
   ! whether or not to print the microscopic cross-sections 
   ! (this can lead to huge file if there are many microscopic XS or number of energy groups is large )
   logical :: printmicro 
   character(len=128) :: arg

   type(core_template) :: core
   type(gen_template) :: gen
   type(iso_template), allocatable :: iso(:)
   type(mat_template), allocatable :: mat(:)
   type(signal_template), allocatable :: signal(:)
   type(stack_template), allocatable :: stack(:)

   ! Default value
   printmicro = .false.
   j = command_argument_count()
   do i = 1, j
      call get_command_argument(i, arg)
      if (index(arg, '--printmicro=') == 1) then
         call parse_logical(arg(14:), printmicro)
      end if
   end do

   write(*,*) 'Running MONTE_extractmicro'
   write(*,*) '    printmicro = ', printmicro

   !  start output file
   open(800, file='output_extractmacro')
   !  save start time
   call system_clock(count_start)
   !  generate timestamp 
   call date_and_time(values=dt)
   write(timestamp, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
         dt(1), dt(2), dt(3), dt(5), dt(6), dt(7)
   write(800,*)gen%separator
   write(800,*)timestamp

   call INP_read_input(core, mat, iso, signal, stack)

   do i = 1, size(signal)
      call SGN_evaluate_signals(signal(i))
   end do

   do i = 1, core%niso
      call XS_read_gendf(iso(i), core)
   end do
   write(*,*)

   write(800,*)'Writing material cross-sections.'
   write(800,*)'To get microscopic XS per material, run MONTE_extractmacro.x --printmicro=true, defaults to false'
   write(800,*)'Hint: scattering matrix is [in,out]'
   write(800,*)gen%separator

   call system('mkdir -p material_xs')
   
   open(1005, file="material_xs/" // "material_list")
   write(1005,*)"Material"
   do i = 1, core%nmat
      write(1005,*)trim(mat(i)%matname)
      call system("mkdir -p material_xs/" // mat(i)%matname)
      
      ! The following method will write microscopic and macroscopic cross-sections 
      call XS_calculate_macro(mat(i), iso, core, signal, printmicro)
   end do
   close(1005)

!  generate timestamp 
   call date_and_time(values=dt)
   write(timestamp, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
         dt(1), dt(2), dt(3), dt(5), dt(6), dt(7)
!  complete output file
   write(800,*)timestamp
!  save end time and output elapsed time
   call system_clock(count_rate=count_rate)
   call system_clock(count_end)
   write(800,*)'Elapsed wall time (s):',real(count_end - count_start) / real(count_rate)
   close(800)

   contains
   subroutine parse_logical(str, val)
      character(len=*), intent(in) :: str
      logical, intent(out) :: val

      select case (adjustl(str))
      case ('true', 't', '1', 'yes', 'y')
         val = .true.
      case ('false', 'f', '0', 'no', 'n')
         val = .false.
      case default
         write(*,*) 'Invalid logical value:', trim(str)
         stop 1
      end select
   end subroutine parse_logical

end program main
