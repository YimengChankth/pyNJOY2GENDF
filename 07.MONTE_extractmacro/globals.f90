module globals
   use, intrinsic :: iso_fortran_env
   implicit none

!--------------------------------------------------------------------------------------------------
!  CORE GENERAL DATA:
   type :: core_template
      integer :: niso = 0                           !total number of isotopes
      integer :: nmat = 0                           !total number of materials
      integer :: nstk = 0                           !total number of stacks
      integer :: nx = 0                             !x-dimension of core map matrix
      integer :: ny = 0                             !y-dimension of core map matrix
      integer :: nz = 0                             !number of axial nodes
      integer :: ng = 0                             !number of energy groups
      integer, allocatable :: imap(:,:,:)           !core map (nx, ny, nz) in terms of indexes of core materials as in matname
      integer :: num_neutrons_born                  !number of neutrons at the beginning of a generation (batch size)
      integer :: num_cycles_inactive                !number of inactive cycles
      integer :: num_cycles_active                  !number of active cycles
      real(real64) :: pitch                                 !hexagonal subassembly pitch (cm)
      real(real64), allocatable :: dz(:)                    !thicknesses of axial layers (nz)
      real(real64), allocatable :: flux(:,:,:,:)            !scalar neutron flux (nx, ny, nz, ng)
      real(real64), allocatable :: power(:,:,:)             !nodal power (nx, ny, nz)
      real(real64), allocatable :: powerxy(:,:)             !plane power (nx, ny)
      real(real64), allocatable :: keff_active_cycle(:)     !multiplication factor of a cycle (num_cycles_active)
      real(real64), allocatable :: keff_expected(:)         !multiplication factor of a problem (num_cycles_active)
      real(real64), allocatable :: sigma_keff(:)            !standard deviation of k-eff (num_cycles_active)
      character*100, allocatable :: xymap(:,:)      !core map (nx, ny)
      character*100, allocatable :: isoname(:)      !total list of isotopes as the GENDF file name (niso)
      character*100, allocatable :: matname(:)      !total list of core materials (nmat)
      character(:), allocatable :: path2data        !path to the nuclear data GENDF file
   end type core_template

!--------------------------------------------------------------------------------------------------
!  GENERAL DATA:
   type :: gen_template
      character*100 :: separator = '-------------------------------------------------------------------------------------'
      real(real64) :: start_time, end_time
   end type gen_template

!--------------------------------------------------------------------------------------------------
!  ISOTOPE NUCLEAR DATA:
   type :: iso_template
      integer :: ng                                 !number of energy groups
      integer :: nins                               !number of inelastic scattering channels
      integer :: nl                                 !number of Legendre components
      integer :: nonze                              !number of nonzeros in elastic scattering matrix
      integer :: nonzi                              !number of nonzeros in inelastic scattering matrix
      integer :: nsig0                              !number of background cross sectopms in XS parametrization
      integer :: ntemp                              !number of temperatures in XS parametrization
      integer, allocatable :: f_e(:)                !index of from-group for elastic scattering (nonze)
      integer, allocatable :: f_i(:)                !index of from-group for inelastic scattering (nonzi)
      integer, allocatable :: t_e(:)                !index of to-group for elastic scattering (nonze)
      integer, allocatable :: t_i(:)                !index of to-group for inelastic scattering (nonzi)
      real(real64) :: temper                                !isotope temperature
      real(real64), allocatable :: chi(:)                   !fission spectrum (ng)
      real(real64), allocatable :: kerma(:,:,:)             !kerma-factors (ng)
      real(real64), allocatable :: flux(:,:,:)            !group flux (nl, ng, ntemp, nsig0)
      real(real64), allocatable :: nubar(:)                 !number of neutrons per fission (ng)
      real(real64), allocatable :: sig0(:)                  !background cross sections in XS parametrization (nsig0)
      real(real64), allocatable :: sige(:,:,:,:)            !elastic scattering cross section (nl, ntemp, nonze, nsig0)
      real(real64), allocatable :: sigf(:,:,:)              !fission cross section (ng, ntemp, nsig0)
      real(real64), allocatable :: sigi(:,:,:)              !inelastic scattering cross section (nonzi, nsig0)
      real(real64), allocatable :: sigt(:,:,:)              !total cross section (ng, ntemp, nsig0)
      real(real64), allocatable :: temp(:)                  !values of temperatures in XS parametrization
      character(:), allocatable :: isoname          !name of the isotope as the GENDF file name
   end type iso_template

!--------------------------------------------------------------------------------------------------
!  CORE MATERIAL NUCLEAR DATA:
   type :: mat_template
      integer :: niso = 0                           !number of isotopes in the material
      real(real64), allocatable :: chi(:)                   !fission spectrum (ng)
      real(real64) :: d_void = 0.0                          !diameter of the central void
      real(real64), allocatable :: numden(:)                !number density (niso)
      real(real64), allocatable :: sig0(:,:)                !background cross sections (niso, ng)
      real(real64), allocatable :: SigA(:)                  !macroscopic absorption cross section (ng)
      real(real64), allocatable :: SigF(:)                  !macroscopic fission cross section (ng)
      real(real64), allocatable :: SigP(:)                  !macroscopic production cross section (ng)
      real(real64), allocatable :: SigT(:)                  !macroscopic total cross section (ng)
      real(real64), allocatable :: SigS(:,:,:)              !macroscopic scattering cross section (nl,ng,ng)
      real(real64), allocatable :: SigST(:)                 !macroscopic total scattering cross section (ng)
      real(real64), allocatable :: Kerma(:)                 !macroscopic kerma factors (ng)
      character(:), allocatable :: matname          !name of the material
      character*100, allocatable :: isoname(:)      !list of names of the isotope in the material (niso)
      character*100, allocatable :: tsignal(:)      !list of names of the signals with isotope temperature (niso)
   end type mat_template

!--------------------------------------------------------------------------------------------------
!  NEUTRON PARAMETERS:
   type :: n_template
      integer :: ig                                 !energy group number (-)
      real(real64) :: x                                     !x-position (cm)
      real(real64) :: y                                     !y-position (cm)
      real(real64) :: z                                     !z-position (cm)
      real(real64) :: mu_x                                  !cosine of angle between neutron direction and x-axis (-)
      real(real64) :: mu_y                                  !cosine of angle between neutron direction and y-axis (-)
      real(real64) :: mu_z                                  !cosine of angle between neutron direction and z-axis (-)
      real(real64) :: weight                                !neutron weight (-)
      real(real64) :: weight0                               !neutron weight at the beginning of cycle (-)
   end type n_template

!--------------------------------------------------------------------------------------------------
!  SIGNAL:
   type :: signal_template
      real(real64) :: r                                     !input value
      real(real64) :: v                                     !signal value                          
      real(real64), allocatable :: x(:), y(:)               !look-up table
      character(:), allocatable :: sgn1, sgn2       !name of the input signals
      character(:), allocatable :: sgnname          !name of the signal
      character(:), allocatable :: sgntype          !type of the signal (CONSTANT)
   end type signal_template

!--------------------------------------------------------------------------------------------------
!  STACK: axial stack of materials describing a core subassembly
   type :: stack_template
      integer :: nz = 0                             !number of axial nodes
      character(:), allocatable :: stackname        !name of the stack
      character*100, allocatable :: matname(:)      !list of names of the materials in the stack (nz)
   end type stack_template


end module globals
