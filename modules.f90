!
!             Author:   Morten Hjorth-Jensen
!             ADDRESS:  Dept. Physics, University Oslo, N-0316 OSLO
!             E-MAIL:   morten.hjorth-jensen@fys.uio.no
!             LANGUAGE: F90/F95
!             LAST UPGRADE : April 2005
!
! MODULE kind_params
!
! Standard Fortran kind parameters used throughout the codebase.
! Placed at the top of modules.f90 so that every subsequent MODULE
! can USE it without circular dependencies.
!
MODULE kind_params
  IMPLICIT NONE
  INTEGER, PARAMETER, PUBLIC :: dp  = KIND(1.0D0)
  INTEGER, PARAMETER, PUBLIC :: dpc = KIND((1.0D0,1.0D0))
  INTEGER, PARAMETER, PUBLIC :: i8  = SELECTED_INT_KIND(18)
END MODULE kind_params


MODULE parallel
  USE omp_lib
  USE kind_params
  INCLUDE 'mpif.h'
  INTEGER, PUBLIC :: from, to, dest, src, ierror, iam, num_procs, master
  INTEGER, PUBLIC :: threads, size0
  INTEGER, allocatable :: subcomm_t3(:), group_t3(:,:)
  INTEGER, allocatable :: subcomm_paeom3(:), group_paeom3(:,:)
  INTEGER, allocatable :: subcomm_preom3(:), group_preom3(:,:)
  INTEGER, allocatable :: members(:), group(:), rank(:), requests(:), array_of_statuses(:,:)
  INTEGER :: status(MPI_STATUS_SIZE) 
  INTEGER :: mystatus(MPI_STATUS_SIZE)
  TYPE, PUBLIC :: topology
     INTEGER(i8) :: my_size, my_start, my_stop
     INTEGER(i8), DIMENSION(:), POINTER :: all_starts, all_stops, all_sizes
  end TYPE topology
  TYPE(topology) :: eom
end MODULE parallel


MODULE constants
  USE kind_params
  ! FROM NIST CODATA
  ! proton mass = 938.27208816(29) MeV
  ! neutron mass = 939.56542052(54) MeV
  ! hbarc = 197.3269804 MeV.fm
  ! pi = 3.1415926535897932
  ! FROM PDG 2022
  ! pion(0)_mass = 134.9768 MeV
  ! pion(+/-1)_mass = 139.57039 MeV
  REAL(dp), PARAMETER, PUBLIC :: pi = 3.1415926535897932_dp
  REAL(dp), PARAMETER, PUBLIC :: hbarc = 197.3269804_dp
  REAL(dp), PARAMETER, PUBLIC :: p_mass = 938.27208816_dp
  REAL(dp), PARAMETER, PUBLIC :: n_mass = 939.56542052_dp
  REAL(dp), PARAMETER, PUBLIC :: nuc_mass = 0.5_dp*(n_mass + p_mass)
  INTEGER, PUBLIC :: Nmax0, NNmin,NNmax, n3min,n3max, t3min,t3max
  INTEGER, PUBLIC :: n1max, N2max
  INTEGER, PUBLIC :: szmin,szmax, tzmin,tzmax, cut_3b
  INTEGER, PUBLIC :: tot_orbs, tot_orbs_n, tot_orbs_p
  INTEGER, PUBLIC :: below_ef, below_ef_n, below_ef_p
  INTEGER, PUBLIC :: test0
  REAL(dp), PUBLIC :: cc_scale, cc_level
  REAL(dp), PUBLIC :: volume, lx, ly, lz
  REAL(dp), PUBLIC :: kf, kf_n,kf_p, rho_n,rho_p
  COMPLEX(dpc) :: e0, ecorr2, ecorr3, ecorr, eccdt
  LOGICAL :: t3_switch, pre_gs0
  INTEGER, DIMENSION(-1:1), PUBLIC :: NF2
  CHARACTER(len=80) :: cc_file0, interaction
  ! eom stuff
  INTEGER, PUBLIC :: eom_states, eom_iterations
  INTEGER, PUBLIC :: nx_paeom, ny_paeom, nz_paeom, tz_paeom
  INTEGER, PUBLIC :: nx_preom, ny_preom, nz_preom, tz_preom
  INTEGER(i8), PUBLIC :: eom_ndim2, eom_ndim3, eom_ndim
  ! print stuff
  CHARACTER(len=29) :: bsp29
  !
  INTEGER, PUBLIC :: Nmax=4                    ;namelist /args/ Nmax
  INTEGER, PUBLIC :: Nocc=2                    ;namelist /args/ Nocc
  INTEGER, PUBLIC :: Pocc=0                    ;namelist /args/ Pocc
  REAL(dp),  PUBLIC :: rho=0.16d0                ;namelist /args/ rho
  INTEGER, PUBLIC :: cc_approx=0               ;namelist /args/ cc_approx
  INTEGER, PUBLIC :: tnf_approx=0              ;namelist /args/ tnf_approx
  INTEGER, PUBLIC :: t3_cut=1000               ;namelist /args/ t3_cut
  INTEGER, PUBLIC :: add_n=0                   ;namelist /args/ add_n
  INTEGER, PUBLIC :: add_p=0                   ;namelist /args/ add_p
  INTEGER, PUBLIC :: eom_approx=0              ;namelist /args/ eom_approx
  INTEGER, PUBLIC :: eom3_cut=1000             ;namelist /args/ eom3_cut
  INTEGER, PUBLIC :: pre_gs=0                  ;namelist /args/ pre_gs
  INTEGER, PUBLIC :: test=0                    ;namelist /args/ test
  INTEGER, PUBLIC :: debug_contracts_flag=1    ;namelist /args/ debug_contracts_flag
  CHARACTER(LEN=80), PUBLIC :: cc_file='ccinf' ;namelist /args/ cc_file
  
CONTAINS

  SUBROUTINE initialize_memory
    IMPLICIT NONE
    INTEGER :: ii
    DO ii = 1, 29
       bsp29(ii:ii) = char(8)
    end DO
  END SUBROUTINE initialize_memory

  
end MODULE constants


MODULE one_body_operators
  USE kind_params
  REAL(dp), PUBLIC, allocatable, dimension(:) :: x_mesh, wx_mesh
  COMPLEX(dpc), PUBLIC, allocatable, dimension(:,:,:,:) :: v2body
END MODULE one_body_operators


MODULE diis_mod
  USE kind_params
  INTEGER, PUBLIC :: n_diis, nn_diis, nstep, diis_step, diis_subspace 
  COMPLEX(dpc), PUBLIC, allocatable, dimension(:,:) :: t2_diis
  COMPLEX(dpc), PUBLIC, allocatable, dimension(:,:) :: t2_diisd
  COMPLEX(dpc), PUBLIC, allocatable, dimension(:,:) :: diis_mat
  COMPLEX(dpc), PUBLIC, allocatable, dimension(:)   :: sol_vect
end module diis_mod


MODULE single_particle_orbits
  USE kind_params
  TYPE, PUBLIC :: single_particle_descript
     INTEGER :: total_orbits
     INTEGER, DIMENSION(:), POINTER :: nx, ny, nz, tz, sz
     CHARACTER (LEN=10), DIMENSION(:), POINTER :: orbit_status, orb_type, model_space
     REAL(dp), DIMENSION(:), POINTER :: e, kx, ky, kz
     
  END TYPE single_particle_descript
  TYPE (single_particle_descript), PUBLIC :: all_orbit, neutron_data, proton_data
  
CONTAINS
  SUBROUTINE allocate_sp_array(this_array,n)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    INTEGER , INTENT(IN) :: n
    INTEGER :: I
    IF (ASSOCIATED (this_array%kx) ) DEALLOCATE(this_array%kx)
    ALLOCATE(this_array%kx(n))
    IF (ASSOCIATED (this_array%ky) ) DEALLOCATE(this_array%ky)
    ALLOCATE(this_array%ky(n))
    IF (ASSOCIATED (this_array%kz) ) DEALLOCATE(this_array%kz)
    ALLOCATE(this_array%kz(n))
    IF (ASSOCIATED (this_array%nx) ) DEALLOCATE(this_array%nx)
    ALLOCATE(this_array%nx(n))
    IF (ASSOCIATED (this_array%ny) ) DEALLOCATE(this_array%ny)
    ALLOCATE(this_array%ny(n))
    IF (ASSOCIATED (this_array%nz) ) DEALLOCATE(this_array%nz)
    ALLOCATE(this_array%nz(n))
    IF (ASSOCIATED (this_array%tz) ) DEALLOCATE(this_array%tz)
    ALLOCATE(this_array%tz(n))
    IF (ASSOCIATED (this_array%sz) ) DEALLOCATE(this_array%sz)
    ALLOCATE(this_array%sz(n))
    IF (ASSOCIATED (this_array%e) ) DEALLOCATE(this_array%e)
    ALLOCATE(this_array%e(n))
    IF (ASSOCIATED (this_array%orbit_status) ) DEALLOCATE(this_array%orbit_status)
    ALLOCATE(this_array%orbit_status(n))
    
    ! blank all characters and zero all other values
    DO i= 1, n
       this_array%orbit_status(i)= ' '
       this_array%e(i)=0.
       this_array%kx(i)=0.
       this_array%ky(i)=0.
       this_array%kz(i)=0.
       this_array%nx(i)=0
       this_array%ny(i)=0
       this_array%nz(i)=0
       this_array%tz(i)=0
       this_array%sz(i)=0
    ENDDO

  END SUBROUTINE allocate_sp_array

  SUBROUTINE deallocate_sp_array(this_array)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    DEALLOCATE(this_array%kx); DEALLOCATE(this_array%ky); DEALLOCATE(this_array%kz);
    DEALLOCATE(this_array%nx); DEALLOCATE(this_array%ny); DEALLOCATE(this_array%nz);
    DEALLOCATE(this_array%sz); DEALLOCATE(this_array%tz);
    DEALLOCATE(this_array%e); DEALLOCATE(this_array%orbit_status)
  END SUBROUTINE deallocate_sp_array

END MODULE single_particle_orbits



MODULE operator_storage
  USE kind_params

  TYPE, PUBLIC :: integer_storage
     INTEGER :: number_confs
     INTEGER, DIMENSION(:), ALLOCATABLE :: ival1
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: ival2
     INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ival3
  end TYPE integer_storage

  TYPE, PUBLIC :: block_storage
     REAL(dp), DIMENSION(:,:), ALLOCATABLE :: val
     COMPLEX(dpc), DIMENSION(:,:), ALLOCATABLE :: cval
     COMPLEX(dpc), DIMENSION(:), ALLOCATABLE :: cval1
     COMPLEX(dpc), DIMENSION(:,:,:), ALLOCATABLE :: cval3
  end TYPE block_storage

  TYPE, PUBLIC :: t3_block_view
     COMPLEX(dpc), DIMENSION(:,:), POINTER :: cval => null()
  END TYPE t3_block_view

  TYPE, PUBLIC :: complex1d_storage
     COMPLEX(dpc), DIMENSION(:), POINTER :: buf => null()
  END TYPE complex1d_storage

  TYPE, PUBLIC :: t3_superblock_storage
     TYPE(t3_block_view),   DIMENSION(:,:), ALLOCATABLE :: val2
     TYPE(complex1d_storage), DIMENSION(:), ALLOCATABLE :: pack1
  END TYPE t3_superblock_storage
  
  TYPE, PUBLIC :: superblock_storage
     TYPE(block_storage), DIMENSION(:), ALLOCATABLE :: val1
     TYPE(block_storage), DIMENSION(:,:), ALLOCATABLE :: val2
     TYPE(integer_storage), DIMENSION(:), ALLOCATABLE :: ival1
  END TYPE superblock_storage

  TYPE, PUBLIC :: index_list
     INTEGER :: n = 0
     INTEGER, ALLOCATABLE :: idx(:)
  END TYPE index_list
  
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: u_hf, tkin_mtx
  COMPLEX(dpc), DIMENSION(:,:), ALLOCATABLE :: fock_mtx
  
  INTEGER, ALLOCATABLE, PUBLIC :: number_2b(:,:)
  INTEGER, ALLOCATABLE, PUBLIC :: number_2bcross(:,:)

  INTEGER :: ch3_min, ch3_max
  INTEGER, ALLOCATABLE, PUBLIC :: climit_t3(:), klimit_t3(:)
  INTEGER, ALLOCATABLE, PUBLIC :: climits_t3(:,:), klimits_t3(:,:)
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: clist_t3(:), klist_t3(:)
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: mapping_t3(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_t3_red(:,:)
  INTEGER, ALLOCATABLE, PUBLIC :: num0_t3(:,:,:,:,:)
  TYPE (index_list), ALLOCATABLE, PUBLIC :: t3_hh_inv(:,:), t3_hp_inv(:,:)

  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_2b_configs(:,:)
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_2bcross_configs(:,:)
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_t3_configs(:,:)
  
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: t2_ccm(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: t2_ccm_eqn(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: t2_ccm_cross(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: t2_ccm_eqn_cross(:)

  TYPE (t3_superblock_storage), ALLOCATABLE, PUBLIC :: t3_ccm0(:)
  TYPE (t3_superblock_storage), ALLOCATABLE, PUBLIC :: t3_ccm(:)
  
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: number_2b_t3(:)
  TYPE (superblock_storage), ALLOCATABLE, PUBLIC :: hh_config_t3(:)
  TYPE (superblock_storage), ALLOCATABLE, PUBLIC :: pp_config_t3(:)
  INTEGER, ALLOCATABLE, PUBLIC :: number_3b_t3(:,:)
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER, ALLOCATABLE, PUBLIC :: r1_paeom_ind(:), r2_paeom_ind(:), r2_paeom_ind_cross(:)  
  COMPLEX(dpc), DIMENSION(:), ALLOCATABLE :: paeom_eigs
  COMPLEX(dpc), DIMENSION(:), ALLOCATABLE :: r1_paeom
  COMPLEX(dpc), DIMENSION(:), ALLOCATABLE :: r1_paeom_eqn
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: r2_paeom(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: r2_paeom_eqn(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: r2_paeom_cross(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: r2_paeom_eqn_cross(:)
  TYPE (superblock_storage), ALLOCATABLE, PUBLIC :: r3_paeom(:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar3b_paeom_I3(:)

  INTEGER, ALLOCATABLE, PUBLIC :: ch1_paeom3(:,:)
  INTEGER, ALLOCATABLE, PUBLIC :: ch2_paeom3(:)
  INTEGER, ALLOCATABLE, PUBLIC :: climit_paeom3(:)
  INTEGER, ALLOCATABLE, PUBLIC :: climits_paeom3(:,:)
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: clist_paeom3(:)
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: mapping_paeom3(:)
  INTEGER :: ch3_paeom_min, ch3_paeom_max
  
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_paeom_pph(:,:,:), check_my_channel_paeom_pph(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_paeom_php_cross(:,:,:), check_my_channel_paeom_php_cross(:)
  
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: r1_paeom_test(:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: r1_paeom_eqn_test(:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: r2_paeom_test(:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: r2_paeom_eqn_test(:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: r3_paeom_test(:,:,:,:,:)
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER, ALLOCATABLE, PUBLIC :: r1_preom_ind(:), r2_preom_ind(:), r2_preom_ind_cross(:)  
  COMPLEX(dpc), DIMENSION(:), ALLOCATABLE :: preom_eigs
  COMPLEX(dpc), DIMENSION(:), ALLOCATABLE :: r1_preom
  COMPLEX(dpc), DIMENSION(:), ALLOCATABLE :: r1_preom_eqn
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: r2_preom(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: r2_preom_eqn(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: r2_preom_cross(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: r2_preom_eqn_cross(:)
  TYPE (superblock_storage), ALLOCATABLE, PUBLIC :: r3_preom(:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar3b_preom_I2(:)

  INTEGER, ALLOCATABLE, PUBLIC :: ch1_preom3(:)
  INTEGER, ALLOCATABLE, PUBLIC :: klimit_preom3(:)
  INTEGER, ALLOCATABLE, PUBLIC :: klimits_preom3(:,:)
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: klist_preom3(:)
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: mapping_preom3(:)
  INTEGER :: ch3_preom_min, ch3_preom_max
  
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_preom_phh(:,:,:), check_my_channel_preom_phh(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_preom_hhp_cross(:,:,:), check_my_channel_preom_hhp_cross(:)
  
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: r1_preom_test(:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: r1_preom_eqn_test(:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: r2_preom_test(:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: r2_preom_eqn_test(:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: r3_preom_test(:,:,:,:,:)

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: v2b_pppp(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: v2b_hhhh(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: v2b_hphp(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: v2b_pphp(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: v2b_hphh(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: v2b_pphh(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: v2b_hphp_cross(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: v2b_phhp_cross(:)
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar1b_I1(:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar1b_I2(:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar1b_I3(:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar1b_I3a(:,:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I1(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I1a(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I2(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I2a(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I3(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I3a(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I4(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I5_cross(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I5a(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I5b(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I5c(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I5d(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I5e_cross(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I6(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I6a(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I7(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I7a(:)
  TYPE (block_storage), ALLOCATABLE, PUBLIC :: hbar2b_I7b(:)
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  TYPE (integer_storage), PUBLIC :: hh_config_2b
  TYPE (integer_storage), PUBLIC :: hp_config_2b
  TYPE (integer_storage), PUBLIC :: pp_config_2b
  TYPE (integer_storage), PUBLIC :: hp_config_2bcross
  TYPE (integer_storage), PUBLIC :: ph_config_2bcross
  
  TYPE (integer_storage), PUBLIC :: hh_channel_2b
  TYPE (integer_storage), PUBLIC :: hp_channel_2b
  TYPE (integer_storage), PUBLIC :: pp_channel_2b
  TYPE (integer_storage), PUBLIC :: hp_channel_2bcross
  TYPE (integer_storage), PUBLIC :: ph_channel_2bcross
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_v2b_pppp(:,:,:), check_my_channel_v2b_pppp(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_v2b_pphp(:,:,:), check_my_channel_v2b_pphp(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_v2b_hphp(:,:,:), check_my_channel_v2b_hphp(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_v2b_pphh(:,:,:), check_my_channel_v2b_pphh(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_v2b_hphh(:,:,:), check_my_channel_v2b_hphh(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_v2b_hhhh(:,:,:), check_my_channel_v2b_hhhh(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_v2b_hphp_cross(:,:,:), check_my_channel_v2b_hphp_cross(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_v2b_phhp_cross(:,:,:), check_my_channel_v2b_phhp_cross(:)
  INTEGER, ALLOCATABLE, PUBLIC :: mapping_v2b_setup(:,:,:), check_my_channel_v2b_setup(:)

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: t1_ccm_test(:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: t1_ccm_eqn_test(:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: t2_ccm_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: t2_ccm_eqn_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: t3_ccm_test(:,:,:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: t3_ccm_eqn_test(:,:,:,:,:,:)  
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: fock_test(:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: v2b_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar1b_I1_test(:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar1b_I2_test(:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar1b_I3_test(:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar1b_I3a_test(:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I1_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I1a_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I2_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I2a_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I3_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I3a_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I4_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I5_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I5a_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I5b_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I5c_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I5d_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I5e_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I6_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I6a_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I7_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I7a_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar2b_I7b_test(:,:,:,:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar3b_paeom_I3_test(:)
  COMPLEX(dpc), ALLOCATABLE, PUBLIC :: hbar3b_preom_I2_test(:)

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end MODULE operator_storage


MODULE configurations
  USE parallel
  USE single_particle_orbits
  USE operator_storage
  
  TYPE configuration_descriptor
     INTEGER :: number_confs
     INTEGER, DIMENSION(:),   POINTER :: config
     INTEGER, DIMENSION(:),   POINTER :: config_NxNyNz_Tz
  end TYPE configuration_descriptor
  
  TYPE (configuration_descriptor), PUBLIC :: channels_2b
  TYPE (configuration_descriptor), PUBLIC :: channels_2bcross
  TYPE (configuration_descriptor), PUBLIC :: channels_t3
  TYPE (configuration_descriptor), PUBLIC :: channels_paeom3
  TYPE (configuration_descriptor), PUBLIC :: channels_preom3
  
CONTAINS
  
  SUBROUTINE number_2b_confs(Nx,Ny,Nz,Tz,struct,this)
    USE constants
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: Nx,Ny,Nz,Tz,struct
    INTEGER :: a,nxa,nya,nza,tza, b,nxb,nyb,nzb,tzb
    INTEGER :: orb_low1, orb_high1, orb_low2, orb_high2, nconfs

    orb_low1 = 1
    orb_high1 = tot_orbs
    orb_low2 = 1
    orb_high2 = tot_orbs
    IF ( struct == 1 ) then
       orb_low1 = 1
       orb_high1 = below_ef
       orb_low2 = 1
       orb_high2 = below_ef
    else IF ( struct == 2 ) then
       orb_low1 = 1
       orb_high1 = below_ef
       orb_low2 = below_ef+1
       orb_high2 = tot_orbs
    else IF ( struct == 3 ) then
       orb_low1 = below_ef+1
       orb_high1 = tot_orbs
       orb_low2 = below_ef+1
       orb_high2 = tot_orbs
    end IF
    
    nconfs = 0
    DO a = orb_low1, orb_high1
       nxa = all_orbit%nx(a)
       nya = all_orbit%ny(a)
       nza = all_orbit%nz(a)
       tza = all_orbit%tz(a)
       DO b = orb_low2, orb_high2
          nxb = all_orbit%nx(b)
          nyb = all_orbit%ny(b)
          nzb = all_orbit%nz(b)
          tzb = all_orbit%tz(b)
          IF ( b <= a ) cycle

          IF ( nxa + nxb /= Nx ) cycle
          IF ( nya + nyb /= Ny ) cycle
          IF ( nza + nzb /= Nz ) cycle
          IF ( tza + tzb /= 2*Tz ) cycle
       
          nconfs=nconfs+1
       end DO
    end DO
    this%number_confs=nconfs

  end SUBROUTINE number_2b_confs

  SUBROUTINE setup_2b_confs(Nx,Ny,Nz,Tz,struct,this)
    USE constants
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: Nx,Ny,Nz,Tz,struct
    INTEGER :: a,nxa,nya,nza,tza, b,nxb,nyb,nzb,tzb, k1,k2
    INTEGER :: orb_low1, orb_high1, orb_low2, orb_high2, nconfs

    orb_low1 = 1
    orb_high1 = tot_orbs
    orb_low2 = 1
    orb_high2 = tot_orbs
    IF ( struct == 1 ) then
       orb_low1 = 1
       orb_high1 = below_ef
       orb_low2 = 1
       orb_high2 = below_ef
    else IF ( struct == 2 ) then
       orb_low1 = 1
       orb_high1 = below_ef
       orb_low2 = below_ef+1
       orb_high2 = tot_orbs
    else IF ( struct == 3 ) then
       orb_low1 = below_ef+1
       orb_high1 = tot_orbs
       orb_low2 = below_ef+1
       orb_high2 = tot_orbs
    end IF

    nconfs = 0
    DO a = orb_low1, orb_high1
       nxa = all_orbit%nx(a)
       nya = all_orbit%ny(a)
       nza = all_orbit%nz(a)
       tza = all_orbit%tz(a)
       DO b = orb_low2, orb_high2
          nxb = all_orbit%nx(b)
          nyb = all_orbit%ny(b)
          nzb = all_orbit%nz(b)
          tzb = all_orbit%tz(b)
          IF ( b <= a ) cycle

          IF ( nxa + nxb /= Nx ) cycle
          IF ( nya + nyb /= Ny ) cycle
          IF ( nza + nzb /= Nz ) cycle
          IF ( tza + tzb /= 2*Tz ) cycle

          nconfs=nconfs+1
          k2=nconfs*2
          k1=k2-1
          this%config(k1)=a
          this%config(k2)=b
       end DO
    end DO
    
    IF ( nconfs /= this%number_confs ) then
       WRITE(6,*) ' Error in configuration allocation ' ; STOP
    end IF

  end SUBROUTINE setup_2b_confs


  SUBROUTINE number_2bcross_confs(Nx,Ny,Nz,Tz,struct,this)
    USE constants
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: Nx,Ny,Nz,Tz,struct
    INTEGER :: a,nxa,nya,nza,tza, b,nxb,nyb,nzb,tzb
    INTEGER :: orb_low1, orb_high1, orb_low2, orb_high2, nconfs

    orb_low1 = 1
    orb_high1 = tot_orbs
    orb_low2 = 1
    orb_high2 = tot_orbs
    IF ( struct == 2 ) then
       orb_low1 = 1
       orb_high1 = below_ef
       orb_low2 = below_ef+1
       orb_high2 = tot_orbs
    else IF ( struct == 3 ) then
       orb_low1 = below_ef+1
       orb_high1 = tot_orbs
       orb_low2 = 1
       orb_high2 = below_ef
    end IF
    
    nconfs = 0
    DO a = orb_low1, orb_high1
       nxa = all_orbit%nx(a)
       nya = all_orbit%ny(a)
       nza = all_orbit%nz(a)
       tza = all_orbit%tz(a)
       DO b = orb_low2, orb_high2
          nxb = all_orbit%nx(b)
          nyb = all_orbit%ny(b)
          nzb = all_orbit%nz(b)
          tzb = all_orbit%tz(b)

          IF ( nxa - nxb /= Nx ) CYCLE
          IF ( nya - nyb /= Ny ) CYCLE
          IF ( nza - nzb /= Nz ) CYCLE
          IF ( tza - tzb /= 2*Tz ) CYCLE
          
          nconfs=nconfs+1
       end DO
    end DO
    this%number_confs=nconfs

  end SUBROUTINE number_2bcross_confs

  SUBROUTINE setup_2bcross_confs(Nx,Ny,Nz,Tz,struct,this)
    USE constants
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: Nx,Ny,Nz,Tz,struct
    INTEGER :: a,nxa,nya,nza,tza, b,nxb,nyb,nzb,tzb, k1,k2
    INTEGER :: orb_low1, orb_high1, orb_low2, orb_high2, nconfs

    orb_low1 = 1
    orb_high1 = tot_orbs
    orb_low2 = 1
    orb_high2 = tot_orbs
    IF ( struct == 2 ) then
       orb_low1 = 1
       orb_high1 = below_ef
       orb_low2 = below_ef+1
       orb_high2 = tot_orbs
    else IF ( struct == 3 ) then
       orb_low1 = below_ef+1
       orb_high1 = tot_orbs
       orb_low2 = 1
       orb_high2 = below_ef
    end IF
    
    nconfs = 0
    DO a = orb_low1, orb_high1
       nxa = all_orbit%nx(a)
       nya = all_orbit%ny(a)
       nza = all_orbit%nz(a)
       tza = all_orbit%tz(a)
       DO b = orb_low2, orb_high2
          nxb = all_orbit%nx(b)
          nyb = all_orbit%ny(b)
          nzb = all_orbit%nz(b)
          tzb = all_orbit%tz(b)

          IF ( nxa - nxb /= Nx ) CYCLE
          IF ( nya - nyb /= Ny ) CYCLE
          IF ( nza - nzb /= Nz ) CYCLE
          IF ( tza - tzb /= 2*Tz ) CYCLE
          
          nconfs=nconfs+1
          k2=nconfs*2
          k1=k2-1
          this%config(k1)=a
          this%config(k2)=b
       end DO
    end DO
    
    IF ( nconfs /= this%number_confs ) then
       WRITE(6,*) ' Error in configuration allocation ' ; STOP
    end IF
    
  end SUBROUTINE setup_2bcross_confs


  FUNCTION cut3b(p,q,r) result(tf)
    USE constants
    USE single_particle_orbits    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: p,q,r
    INTEGER :: Np2, Nq2, Nr2
    INTEGER :: tzp, tzq, tzr
    LOGICAL :: tf

    tf = .false.
    tzp = all_orbit%tz(p)
    tzq = all_orbit%tz(q)
    tzr = all_orbit%tz(r)
    Np2 = all_orbit%nx(p)**2 + all_orbit%ny(p)**2 + all_orbit%nz(p)**2 - NF2(tzp)
    Nq2 = all_orbit%nx(q)**2 + all_orbit%ny(q)**2 + all_orbit%nz(q)**2 - NF2(tzq)
    Nr2 = all_orbit%nx(r)**2 + all_orbit%ny(r)**2 + all_orbit%nz(r)**2 - NF2(tzr)
    IF ( Np2 <= 0 ) Np2 = abs(Np2) + 1
    IF ( Nq2 <= 0 ) Nq2 = abs(Nq2) + 1
    IF ( Nr2 <= 0 ) Nr2 = abs(Nr2) + 1
    IF ( Np2 + Nq2 + Nr2 > cut_3b+3 ) tf = .true.
    
  end FUNCTION cut3b

  ! r+pq = ch3
  FUNCTION get_ch_ch3(nx3,ny3,nz3,tz3,struct,p) result(ch)
    USE constants
    USE single_particle_orbits
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nx3,ny3,nz3,tz3, struct, p
    INTEGER :: nxp,nyp,nzp,tzp, q,r
    INTEGER :: Nx1,Ny1,Nz1,Tz1
    INTEGER :: bra,bra_confs, ch_type, ch0, ch

    ch_type = 0
    IF ( struct == 1 ) then
       ch_type = 3
    else IF ( struct == 2 ) then
       ch_type = 1
    end IF
    
    ch = 0
    nxp = all_orbit%nx(p)
    nyp = all_orbit%ny(p)
    nzp = all_orbit%nz(p)
    tzp = all_orbit%tz(p)
    DO ch0 = 1, channels_2b%number_confs
       IF ( number_2b(ch_type,ch0) == 0 ) cycle
       Nx1 = channels_2b%config_NxNyNz_Tz(4*ch0-3)
       Ny1 = channels_2b%config_NxNyNz_Tz(4*ch0-2)
       Nz1 = channels_2b%config_NxNyNz_Tz(4*ch0-1)
       Tz1 = channels_2b%config_NxNyNz_Tz(4*ch0)
       IF ( Nx1 + nxp /= nx3 ) cycle
       IF ( Ny1 + nyp /= ny3 ) cycle
       IF ( Nz1 + nzp /= nz3 ) cycle
       IF ( 2*Tz1 + tzp /= tz3 ) cycle
       bra_confs = number_2b(ch_type,ch0)
       DO bra = 1, bra_confs
          q = lookup_2b_configs(ch_type,ch0)%ival2(1,bra)
          r = lookup_2b_configs(ch_type,ch0)%ival2(2,bra)
          IF ( cut3b(p,q,r) ) cycle
          ch = ch0
          return
       end DO
    end DO
    
  end FUNCTION get_ch_ch3


  SUBROUTINE number_3b_confs(nx3,ny3,nz3,tz3,struct,nconfs)
    USE constants
    USE single_particle_orbits
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nx3,ny3,nz3,tz3, struct
    INTEGER, INTENT(OUT) :: nconfs
    INTEGER :: p,nxp,nyp,nzp,tzp, q,r
    INTEGER :: ch1,Nx1,Ny1,Nz1,Tz1
    INTEGER :: orb_low, orb_high, bra,bra_confs, ch_type
    
    nconfs = 0
    
    orb_low = 1
    orb_high = tot_orbs
    ch_type = 0
    IF ( struct == 1 ) then
       orb_low = below_ef+1
       orb_high = tot_orbs
       ch_type = 3
    else IF ( struct == 2 ) then
       orb_low = 1
       orb_high = below_ef
       ch_type = 1
    end IF

    DO p = orb_low, orb_high
       nxp = all_orbit%nx(p)
       nyp = all_orbit%ny(p)
       nzp = all_orbit%nz(p)
       tzp = all_orbit%tz(p)
       DO ch1 = 1, channels_2b%number_confs
          IF ( number_2b(ch_type,ch1) == 0 ) cycle
          Nx1 = channels_2b%config_NxNyNz_Tz(4*ch1-3)
          Ny1 = channels_2b%config_NxNyNz_Tz(4*ch1-2)
          Nz1 = channels_2b%config_NxNyNz_Tz(4*ch1-1)
          Tz1 = channels_2b%config_NxNyNz_Tz(4*ch1)
          IF ( Nx1 + nxp /= nx3 ) cycle
          IF ( Ny1 + nyp /= ny3 ) cycle
          IF ( Nz1 + nzp /= nz3 ) cycle
          IF ( 2*Tz1 + tzp /= tz3 ) cycle
          bra_confs = number_2b(ch_type,ch1)
          DO bra = 1, bra_confs
             q = lookup_2b_configs(ch_type,ch1)%ival2(1,bra)
             r = lookup_2b_configs(ch_type,ch1)%ival2(2,bra)
             IF ( cut3b(p,q,r) ) cycle
             nconfs = nconfs + 1
          end DO
          exit
       end DO
    end DO

  end SUBROUTINE number_3b_confs


  SUBROUTINE number_3b_red_confs(nx3,ny3,nz3,tz3,struct,nconfs)
    USE constants
    USE single_particle_orbits
    USE operator_storage
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nx3,ny3,nz3,tz3, struct
    INTEGER, INTENT(OUT) :: nconfs
    INTEGER :: p,nxp,nyp,nzp,tzp, q,r
    INTEGER :: ch1,Nx1,Ny1,Nz1,Tz1
    INTEGER :: orb_low, orb_high, bra,bra_confs, ch_type
    
    nconfs = 0
    
    orb_low = 1
    orb_high = tot_orbs
    ch_type = 0
    IF ( struct == 1 ) then
       orb_low = below_ef+1
       orb_high = tot_orbs
       ch_type = 3
    else IF ( struct == 2 ) then
       orb_low = 1
       orb_high = below_ef
       ch_type = 1
    end IF

    DO p = orb_low, orb_high
       nxp = all_orbit%nx(p)
       nyp = all_orbit%ny(p)
       nzp = all_orbit%nz(p)
       tzp = all_orbit%tz(p)
       DO ch1 = 1, channels_2b%number_confs
          IF ( number_2b(ch_type,ch1) == 0 ) cycle
          Nx1 = channels_2b%config_NxNyNz_Tz(4*ch1-3)
          Ny1 = channels_2b%config_NxNyNz_Tz(4*ch1-2)
          Nz1 = channels_2b%config_NxNyNz_Tz(4*ch1-1)
          Tz1 = channels_2b%config_NxNyNz_Tz(4*ch1)
          IF ( Nx1 + nxp /= nx3 ) cycle
          IF ( Ny1 + nyp /= ny3 ) cycle
          IF ( Nz1 + nzp /= nz3 ) cycle
          IF ( 2*Tz1 + tzp /= tz3 ) cycle
          bra_confs = number_2b(ch_type,ch1)
          DO bra = 1, bra_confs
             q = lookup_2b_configs(ch_type,ch1)%ival2(1,bra)
             r = lookup_2b_configs(ch_type,ch1)%ival2(2,bra)
             IF ( cut3b(p,q,r) ) cycle
             IF ( p < q ) nconfs = nconfs + 1
          end DO
          exit
       end DO
    end DO
    
  end SUBROUTINE number_3b_red_confs

  
end MODULE configurations


MODULE ang_mom_functions
  USE kind_params
  INTEGER, private, parameter:: maxjj=100
  REAL(dp), PRIVATE :: f_mb(maxjj),g_mb(maxjj),w_mb(maxjj)
  INTEGER, PRIVATE :: kh(4*maxjj)
  REAL(dp), PARAMETER, PRIVATE :: pi=3.141592654
  REAL(dp), PRIVATE :: q(maxjj,maxjj), cn(0:maxjj+1,0:maxjj+1)
 
CONTAINS
  !
  !     factorials for 3j,6j and 9j symbols            
  !     for moshinsky trans brackets and for           
  !     vector brackets                                
  !
  SUBROUTINE commons_to_angmom
    IMPLICIT NONE
    INTEGER :: l, k, i, j
    REAL(dp) :: a , sq_pi, fj, tfj, fk
    !    3j, 6j and 9j symbols
    kh=1
    kh(2*maxjj) =0
    DO l=1,maxjj
       q(l,1)=1.0d0
       q(l,l)=1.0d0
       kh(l+l+2*maxjj)=0
    end DO
    DO l=2,maxjj-1
       DO k=2,l
          q(l+1,k)=q(l,k-1)+q(l,k)
       end DO
    end DO
    !    Moshinsky brackets
    f_mb(1)=0.
    g_mb(1)=LOG(0.5D0)
    w_mb(1)=0.
    DO i=2,maxjj
       a=i-1
       f_mb(i)=f_mb(i-1)+LOG(a)
       g_mb(i)=g_mb(i-1)+LOG(a+0.5D0)
       w_mb(i)=LOG(a+a+1.)
    end DO
    !    spherical harmonics
    cn=0.
    sq_pi=1./SQRT(2.*pi)
    DO j=0,maxjj+1
       cn(0,j)=SQRT(0.5*(2.*j+1.))
    end DO
    DO j=1,maxjj+1
       tfj=2.*j
       cn(j,j)=cn(j-1,j-1)*SQRT((tfj+1.)/tfj)
    end DO
    DO j=0,maxjj+1
       fj=FLOAT(j)
       DO k=1,j-1
          fk=FLOAT(k)
          cn(k,j)=cn(k-1,j)*SQRT((fj+fk)*(fj-fk+1.))*0.5/fk
       end DO
    end DO
    cn=cn*sq_pi

  end SUBROUTINE commons_to_angmom
  
  ! Calculates CG-symbols
  REAL(dp) FUNCTION cgc(j_a,j_b,j_c,m_a,m_b,m_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j_a,j_b,j_c,m_a,m_b,m_c
    INTEGER :: iph 

    cgc = tjs(j_a,j_b,j_c,m_a,m_b,-m_c) * iph( (j_a-j_b+m_c)/2) * sqrt(j_c+1.) 
  end FUNCTION cgc

  ! Calculates 3j-symbols
  REAL(dp) FUNCTION tjs(j_a,j_b,j_c,m_a,m_b,m_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j_a,j_b,j_c,m_a,m_b,m_c
    INTEGER :: ja, jb, jc, mb, ma, mc, la, lb, lc, lt, ld, ja2, jb2, &
         jc2, i, k0, k1, k, ip
    REAL(dp) :: x, fn, p
    logical :: triag    
    
    tjs=0.
    IF ( triag( j_a, j_b, j_c) ) return 
    IF ( m_a + m_b + m_c /= 0 ) return

    ja=(j_a+m_a)/2+1
    ma=(j_a-m_a)/2+1
    jb=(j_b+m_b)/2+1
    mb=(j_b-m_b)/2+1
    jc=(j_c+m_c)/2+1
    mc=(j_c-m_c)/2+1
    la=(j_b+j_c-j_a)/2+1
    lb=(j_c+j_a-j_b)/2+1
    lc=(j_a+j_b-j_c)/2+1
    lt=(j_a+j_b+j_c)/2+1
    ld=MIN(ja,jb,jc,ma,mb,mc,la,lb,lc)
    IF(((m_a+m_b+m_c) <= 0).AND.(ld > 0)) then
       ja2=j_a+m_a
       jb2=j_b+m_b
       jc2=j_c+m_c
       i=ja2+jb2+jc2-ja2/2*2-jb2/2*2-jc2/2*2
       IF(i == 0) then
          fn=q(ja+ma-1,lc)*q(jb+mb-1,lc)/(q(lt,jc+mc-1)*q(lt+1,2) &
               *q(ja+ma-1,ja)*q(jb+mb-1,jb)*q(jc+mc-1,jc))
          k0=MAX(0,lc-ja,lc-mb)+1
          k1=MIN(lc,ma,jb)
          x=0.
          DO k=k0,k1
             x=-x-q(lc,k)*q(lb,ma-k+1)*q(la,jb-k+1)
          end DO
          ip=k1+lb+jc
          p=1-2*(ip-ip/2*2)
          tjs=p*x*SQRT(fn)
       end IF
    end IF
  end FUNCTION tjs

  ! Calculates 6j-symbols
  REAL(dp) FUNCTION sjs(j_a,j_b,j_c,l_a,l_b,l_c)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j_a,j_b,j_c,l_a,l_b,l_c
    INTEGER :: ja,jb,jc,la,lb,lc,i,mt,ma,mb,mc,na,nb,nc,ka,&
         kb,kc,l,l0,l1, ihlp, idum
    REAL(dp) :: x, fs, fss
    LOGICAL :: triag
    
    ihlp=2*maxjj-1

    sjs=0.0d0
    
    IF ( triag( j_a, j_b, j_c) ) return 
    IF ( triag( j_c, l_a, l_b) ) return 
    IF ( triag( j_a, l_b, l_c) ) return 
    IF ( triag( j_b, l_a, l_c) ) return 
    
    ja=j_a + 1
    jb=j_b + 1
    jc=j_c + 1
    la=l_a + 1
    lb=l_b + 1
    lc=l_c + 1
    i=kh(ja+jb-jc+ihlp)+kh(jb+jc-ja+ihlp)+kh(jc+ja-jb+ihlp)+kh(ja+lb-lc+ihlp) &
         +kh(lb+lc-ja+ihlp)+kh(lc+ja-lb+ihlp)+kh(la+jb-lc+ihlp)+kh(jb+lc-la+ihlp) &
         +kh(lc+la-jb+ihlp)+kh(la+lb-jc+ihlp)+kh(lb+jc-la+ihlp)+kh(jc+la-lb+ihlp)
    IF(i <= 0) then
       mt=(j_a+j_b+j_c)/2 + 2
       ma=(j_a+l_b+l_c)/2+ 2
       mb=(l_a+j_b+l_c)/2+ 2
       mc=(l_a+l_b+j_c)/2+ 2
       na=mt-ja
       nb=mt-jb
       nc=mt-jc
       ka=ma-lc
       kb=mb-lc
       kc=mc-jc

       idum=max(mt,ja+1,nc,ma,ka,mb,la+1,kb,mc,kc)

       IF(idum.gt.maxjj) then
          write(6,*) 'increase maxjj in MODULE ang_mom_functions from', maxjj, 'to', idum
          stop
       end IF

       fss=q(mt,ja+1)*q(ja,nc)/(q(ma,ja+1)*q(ja,ka)*q(mb,la+1)* &
            q(la,kb)*q(mc,la+1)*q(la,kc))
       fs=SQRT(fss)/(l_a + 1.)
       l0=MAX(mt,ma,mb,mc)+1
       l1=MIN(ma+na,mb+nb,mc+nc)
       x=0.
       DO l=l0,l1
          x=-x+q(l-1,mt)*q(na,l-ma)*q(nb,l-mb)*q(nc,l-mc)
       end DO
       sjs=-(1+2*(l1/2*2-l1))*fs*x
    end IF

  end FUNCTION sjs
  
  ! Calculates ninej-symbols
  REAL(dp) FUNCTION snj (ia,ib,ie,ic,id,IF,ig,ih,it)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ia,ib,ie,ic,id,IF,ig,ih,it
    INTEGER :: ja,jb,je,jc,jd,jf,jg,jh,jt,i,la,ld,ma,mc,na,nb,le,lf,&
         lg,me,mf,mg,ne,nf,ng,lx,mx,nx,jsi,jsf, js,is,lb, lc, &
         mb, ly, my,ny,l,l0,m0,n0,l1,m1,n1,m,n,ihx, ihlp
    REAL(dp) :: x, fn, fd, ps, fs, u, y, z, ud, p

    ihlp=2*maxjj-1

    snj=0.
    ja=ia+1
    jb=ib+1
    jc=ic+1
    jd=id+1
    je=ie+1
    jf=IF+1
    jg=ig+1
    jh=ih+1
    jt=it+1
    i=kh(ja+jb-je+ihlp)+kh(jb+je-ja+ihlp)+kh(je+ja-jb+ihlp)+kh(jc+jd-jf+ihlp) &
         +kh(jd+jf-jc+ihlp)+kh(jf+jc-jd+ihlp)+kh(jg+jh-jt+ihlp)+kh(jh+jt-jg+ihlp) &
         +kh(jt+jg-jh+ihlp)+kh(ja+jc-jg+ihlp)+kh(jc+jg-ja+ihlp)+kh(jg+ja-jc+ihlp) &
         +kh(jb+jd-jh+ihlp)+kh(jd+jh-jb+ihlp)+kh(jh+jb-jd+ihlp)+kh(je+jf-jt+ihlp) &
         +kh(jf+jt-je+ihlp)+kh(jt+je-jf+ihlp)
    IF(i <= 0) then
       la=(ie+IF+it)/2+2
       ld=(ig+ih+it)/2+2
       ma=(ia+ic+ig)/2+2
       mc=(IF+ic+id)/2+2
       na=(ib+id+ih)/2+2
       nb=(ib+ie+ia)/2+2
       le=(ie+IF-it)/2+1
       lf=(IF+it-ie)/2+1
       lg=(it+ie-IF)/2+1
       me=(ia+ic-ig)/2+1
       mf=(ic+ig-ia)/2+1
       mg=(ig+ia-ic)/2+1
       ne=(ib+id-ih)/2+1
       nf=(id+ih-ib)/2+1
       ng=(ih+ib-id)/2+1
       lx=(it+ig-ih)/2+1
       mx=(ic+id-IF)/2+1
       nx=(ib+ie-ia)/2+1
       fn=q(la,jt+1)*q(jt,lg)*q(ma,jc+1)*q(jc,mf)*q(na,jb+1)*q(jb,ne)
       fd=q(ld,jt+1)*q(jt,lx)*q(mc,jc+1)*q(jc,mx)*q(nb,jb+1)*q(jb,nx)
       jsi=MAX(abs(je-jh),abs(jg-jf),abs(ja-jd))+1
       jsf=MIN(je+jh,jg+jf,ja+jd)-1
       ps=-1-2*(jsi/2*2-jsi)
       fs=ps*SQRT(fn/fd)/FLOAT((ig+1)*(ie+1))
       u=0.
       DO js=jsi,jsf,2
          is=js-1
          lb=(ie+ih+is)/2+2
          lc=(ig+IF+is)/2+2
          mb=(ia+id+is)/2+2
          ly=(ie+ih-is)/2+1
          my=(ig+IF-is)/2+1
          ny=(ia-id+is)/2+1
          ud=q(lb,je+1)*q(je,ly)*q(lc,jg+1)*q(jg,my)*q(mb,js+1)*q(js,ny)
          l0=MAX(la,lb,lc,ld)+1
          m0=MAX(ma,mb,mc,lc)+1
          n0=MAX(na,nb,mb,lb)+1
          l1=MIN(le+ld,lf+lb,lg+lc)
          m1=MIN(me+lc,mf+mb,mg+mc)
          n1=MIN(ne+lb,nf+nb,ng+mb)
          x=0.
          DO l=l0,l1
             x=-x-q(l-1,la)*q(le,l-ld)*q(lf,l-lb)*q(lg,l-lc)
          end DO
          y=0.
          DO m=m0,m1
             y=-y-q(m-1,ma)*q(me,m-lc)*q(mf,m-mb)*q(mg,m-mc)
          end DO
          z=0.
          DO n=n0,n1
             z=-z-q(n-1,na)*q(ne,n-lb)*q(nf,n-nb)*q(ng,n-mb)
          end DO
          ihx=l1+m1+n1
          p=1+2*(ihx/2*2-ihx)
          u=u+p*x*y*z/ud
       end DO
       snj=u*fs
    end IF

  end FUNCTION snj

  ! Spherical harmonics from Num. Recipes
  DOUBLE PRECISION FUNCTION spherical_harmonics(m1,l,x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m1, l
    DOUBLE PRECISION, INTENT(IN) ::  x
    DOUBLE PRECISION, DIMENSION(0:51) :: y
    DOUBLE PRECISION :: fj, z, fac, div, sum, a, b, c
    INTEGER :: iphase, m, j
    spherical_harmonics=0.
    m=iabs(m1)
    IF(m.LT.0) m=-m1
    y(0)=1.
    IF(l.EQ.0) then
       sum=y(0)
    ELSE
       a=m-l
       b=l+m+1
       c=m+1
       z=0.5-x*0.5
       DO j=1,l-m+1
          fj=j-1
          y(j)=y(j-1)*(a+fj)*(b+fj)*z
          div=(c+fj)*(fj+1.)
          y(j)=y(j)/div
       end DO
       IF(m > 0) then
          fac=(1.-x*x)**m
          fac=SQRT(fac)
       ELSE
          fac=1.
       end IF
       sum=0.
       DO j=0,l-m
          sum=sum+y(j)
       end DO
       iphase=m
       IF(m1.LT.0) then
          iphase=0
       end IF
       sum=sum*fac*((-1)**iphase)
    end IF
    spherical_harmonics=cn(m,l)*sum

  end FUNCTION spherical_harmonics

end MODULE ang_mom_functions
