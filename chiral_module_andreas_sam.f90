
! References (NN forces)
![1] = Machleidt Phys Rep 503 (2011) 1-75 AKA EM2011
![2] = Epelbaum NPA 747 (2005) 362-424    AKA Ep2005
![3] = Machleidt PRC 91 014002 (2015)     AKA EM2015
![4] = Epelbaum EPJA 51, 53 (2015)        AKA Ep2015
![5] = Krebs EPJA 32, 127 (2007)
![6] = Machleidt Phys Scripta 91, 083007 (2016)

MODULE chiral_tables 
  COMPLEX*16, allocatable :: sigma_dot_q_tab(:,:,:,:,:)
  COMPLEX*16, allocatable :: sigma_dot_qxq_tab(:,:,:,:,:)
  COMPLEX*16, allocatable :: sigma_dot_q_ope_tab(:,:,:,:,:)
  COMPLEX*16, allocatable :: sigma1_dot_q_sigma2_dot_q_tab(:,:,:,:,:,:,:)
  COMPLEX*16 :: tau1_dot_tauXtau_tab(-1:1,-1:1,-1:1,-1:1,-1:1,-1:1)
  COMPLEX*16 :: tau_dot_tau_tab(-1:1,-1:1,-1:1,-1:1)
  COMPLEX*16 :: sigma_dot_sigma_tab(-1:1,-1:1,-1:1,-1:1)
  REAL*8 :: delta_tab(-1:1,-1:1)
end MODULE chiral_tables

MODULE chiral_constants 
  USE constants

  ! basic mesh info
  TYPE, PUBLIC :: chp_mesh_info
     INTEGER  :: amount
     REAL(8) :: xmin, xmax
  END TYPE chp_mesh_info
  ! Gauss-Legendre mesh point x, corresponding integration weight w and corresponding x*x*w-value
  TYPE, PUBLIC :: chp_gauleg_mesh_point
     SEQUENCE
     REAL(8) :: x, w, xxw
  END TYPE chp_gauleg_mesh_point
  ! mesh points and weights in momentum space
  TYPE, PUBLIC :: chp_gauleg_mesh
     TYPE(chp_mesh_info)                                    :: info
     TYPE(chp_gauleg_mesh_point), DIMENSION(:), ALLOCATABLE :: pnt
  END TYPE chp_gauleg_mesh
  
  
  ! chiral definitions
  INTEGER, PARAMETER :: LO   = 0
  INTEGER, PARAMETER :: NLO  = 2
  INTEGER, PARAMETER :: NNLO = 3
  INTEGER, PARAMETER :: N3LO = 4
  INTEGER, PARAMETER :: chp_DR = 0
  INTEGER, PARAMETER :: chp_SFR = 1    
  INTEGER, PARAMETER :: pc_machleidt2011 = 0 ! [1] 
  INTEGER, PARAMETER :: pc_epelbaum2005  = 1 ! [2] (and also [3,4,6])
  INTEGER, PARAMETER :: chp_3NF_reg_NONLOCAL = 0
  INTEGER, PARAMETER :: chp_3NF_reg_LOCAL    = 1
  INTEGER, PARAMETER :: chp_3NF_reg_NOREG    = 2
  REAL*8, PARAMETER :: sfr = 700.d0 
  REAL*8, PARAMETER :: sfr2 = sfr*sfr
  CHARACTER(LEN=2), PARAMETER :: chp_itope = 'EM'

  ! chiral parameters
  INTEGER :: which_pot
  INTEGER :: chiral_order
  INTEGER :: nexp_NN                ! exponent of the NN regulator
  REAL*8 :: lambda                  ! NN cutoff
  INTEGER :: regulatization_scheme  ! 0: dimensional regularization (DR), 1: spectral FUNCTION regularization (SFR)
  INTEGER :: power_counting         ! Which power counting scheme to use (they differ in the relativistic corrections proportional to 1/MN)
  REAL*8 :: lambdaChi               ! lambdaChi is the chiral symmetry breaking scale
  INTEGER:: chp_3NF_regulator_type  ! type of regulator: 0: non-local, 1: local, 2: no regulator
  INTEGER:: nreg_3NF_nonloc         ! order n of the non-local regulator
  REAL*8 :: lambda_3NF_nonloc       ! cutoff of the non-local 3N regulator
  LOGICAL :: chiral_delta_flag
  
  ! physical constants
  REAL*8, parameter :: twopi = pi*2.d0
  REAL*8, parameter :: pi2 = pi*pi
  ! FROM NIST CODATA
  ! proton mass = 938.27208816(29) MeV
  ! neutron mass = 939.56542052(54) MeV
  REAL*8, parameter :: proton_mass  = 938.27208816d0
  REAL*8, parameter :: neutron_mass = 939.56542052d0
  REAL*8, parameter :: delta_mass = 1232.d0
  REAL*8 :: gA = 1.29d0 
  REAL*8 :: fpi = 92.4d0 
  REAL*8 :: hA = 1.4d0
  REAL*8 :: gA2, gA4
  REAL*8 :: fpi2, fpi4, fpi_inv
  REAL*8 :: hA2, hA4

  REAL*8 :: mnuc(-1:1), mpi(-1:2)
  REAL*8 :: mnuc_inv(-1:1)! nucleon mass inversed
  REAL*8 :: mnuc2(-1:1)   ! nucleon mass squared
  REAL*8 :: mpi2(-1:2)    ! pion mass squared
  REAL*8 :: mpi3(-1:2)    ! pion mass cubed
  REAL*8 :: mpi4(-1:2)    ! mpi2 squared
  REAL*8 :: mpi5(-1:2)    ! mpi^5
  REAL*8 :: twompi(-1:2)  ! two times pion mass
  REAL*8 :: fourmpi2(-1:2)! four times pion mass squared
  REAL*8 :: delta_split     ! delta_mass - mnuc(0) splitting constant
  REAL*8 :: delta_split2, four_delta_split2
  COMPLEX*16 :: sigma_x(2,2), sigma_y(2,2), sigma_z(2,2)
    
  REAL*8 :: const(100)
  REAL*8 :: cnlo(1:7)
  REAL*8 :: c1,c2,c3,c4
  REAL*8 :: C1_3NF, C2_3NF, C3_3NF, C4_3NF
  REAL*8 :: LEC_c1, lec_c2, lec_c3, lec_c4
  REAL*8 :: LEC_C1_3NF, LEC_C2_3NF, LEC_C3_3NF, LEC_C4_3NF
  REAL*8 :: CE, CD, Econst, Dconst
  REAL*8 :: CS(-1:1), CT(-1:1)
  REAL*8 :: c3_3NF_delta, c4_3NF_delta
  REAL*8 :: lec_c3_3NF_delta, lec_c4_3NF_delta
  
  REAL*8 :: sfr_heavyside(-1:2)
  REAL*8 :: z_low
  REAL*8 :: mu_low, mu_high  
  ! Integration meshes for loop integrals involving a Delta excitation
  TYPE(chp_gauleg_mesh), PRIVATE :: chp_delta_2PE_loop_int_mesh_mu
  INTEGER, PARAMETER, PRIVATE :: nof_2PE_loop_int_mu_points = 25

CONTAINS 

  !
  ! Potentials:
  ! 1: N2LO_opt(500)                Phys. Rev. Lett. 110, 192502 (2013)
  ! 2: NNLO_sat(450)                Phys. Rev. C 91,  051301(R) (2015) 
  ! 3: DeltaNNLOgo(394)             Phys. Rev. C 102, 054301 (2020)
  ! 4: DeltaNNLOgo(450)             Phys. Rev. C 102, 054301 (2020)
  ! 5: Delta NNLO(450)              Phys. Rev. C 97,  024332 (2018)
  ! 6: Delta NNLO(500)              Phys. Rev. C 97,  024332 (2018)
  ! 0: Minnesota Potential          Minnesota Potential
  !
  SUBROUTINE init_chp_constants 
    USE parallel 
    USE constants
    IMPLICIT NONE
    INTEGER :: i
    REAL*8 :: c1s0(-1:1), c3s1(-1:1), cnlo_pw(1:7)

    ! FROM PDG 2022
    ! pion(0)_mass = 134.9768 MeV
    ! pion(+/-1)_mass = 139.57039 MeV
    mpi(-1)  = 139.57039d0
    mpi(0)   = 134.9768d0
    mpi(+1)  = 139.57039d0
    mpi(+2)  = (mpi(-1) + mpi(0) + mpi(+1))/3.d0
    
    mnuc(-1) = proton_mass
    mnuc(0)  = 0.5d0*(proton_mass+neutron_mass)
    mnuc(+1) = neutron_mass
    mnuc_inv = 1.D0 / mnuc
    delta_split = delta_mass - mnuc(0)
    delta_split2 = delta_split * delta_split
    four_delta_split2 = 4.d0 * delta_split2
    
    mnuc2(:)   = mnuc*mnuc 
    mpi2(:)    = mpi*mpi   
    mpi3(:)    = mpi*mpi2
    mpi4(:)    = mpi2*mpi2
    mpi5(:)    = mpi4*mpi 
    twompi(:)  = 2.0D0*mpi 
    fourmpi2(:)= 4.0D0*mpi2

    
    IF ( interaction == "nnloopt" ) then
       ! N2LO_opt(500): PRL 110, 192502 (2013)
       which_pot = 1
       lambda = 500.d0
       chiral_order = NNLO
       nexp_NN = 3
       chiral_delta_flag = .FALSE.
       regulatization_scheme = chp_SFR
       power_counting = pc_machleidt2011
       gA = 1.29d0
       fpi = 92.4d0
       chp_3NF_regulator_type = chp_3NF_reg_LOCAL
       lambda_3NF_nonloc = lambda
       nreg_3NF_nonloc = nexp_NN
       !
       LEC_c1 = -0.9186395d0
       LEC_c2 =  0.0000000d0 
       LEC_c3 = -3.8886875d0
       LEC_c4 =  4.3103272d0
       c1s0(-1) = -0.1513660d0 ! Ct1S0pp
       c1s0(0)  = -0.1521411d0 ! Ct1S0np
       c1s0(1)  = -0.1517647d0 ! Ct1S0nn
       c3s1(-1) = -0.1584342d0 ! Ct3S1pp
       c3s1(0)  = -0.1584342d0 ! Ct3S1np
       c3s1(1)  = -0.1584342d0 ! Ct3S1nn
       cnlo_pw(1) =  2.4040219d0 ! C_1S0 
       cnlo_pw(2) =  1.2633908d0 ! C_3P0
       cnlo_pw(3) =  0.4170455d0 ! C_1P1 
       cnlo_pw(4) = -0.7826585d0 ! C_3P1
       cnlo_pw(5) =  0.9283847d0 ! C_3S1
       cnlo_pw(6) =  0.6181414d0 ! C_3S1 - 3D1
       cnlo_pw(7) = -0.6778085d0 ! C_3P2
       cD =  0.2000000d0
       cE = -0.3600000d0

    else IF ( interaction == "nnlosat" ) then 
       ! NNLO_sat(450): PRC 91, 051301(R) (2015) 
       which_pot = 2
       lambda = 450.d0
       chiral_order = NNLO
       nexp_NN = 3
       chiral_delta_flag = .FALSE.
       regulatization_scheme = chp_SFR
       power_counting = pc_machleidt2011
       gA = 1.29d0
       fpi = 92.4d0
       chp_3NF_regulator_type = chp_3NF_reg_NONLOCAL
       lambda_3NF_nonloc = lambda
       nreg_3NF_nonloc = nexp_NN
       !
       LEC_c1 = -1.1215212d0  
       LEC_c2 =  0.0000000d0
       LEC_c3 = -3.9250059d0 
       LEC_c4 =  3.7656872d0
       c1s0(-1) = -0.1581494d0 ! Ct1S0pp
       c1s0(0)  = -0.1598224d0 ! Ct1S0np
       c1s0(1)  = -0.1591504d0 ! Ct1S0nn
       c3s1(-1) = -0.1776744d0 ! Ct3S1pp
       c3s1(0)  = -0.1776744d0 ! Ct3S1np
       c3s1(1)  = -0.1776744d0 ! Ct3S1nn
       cnlo_pw(1) =  2.5393678d0 ! C_1S0 
       cnlo_pw(2) =  1.3983656d0 ! C_3P0
       cnlo_pw(3) =  0.5559588d0 ! C_1P1 
       cnlo_pw(4) = -1.1360953d0 ! C_3P1
       cnlo_pw(5) =  1.0028927d0 ! C_3S1
       cnlo_pw(6) =  0.6007160d0 ! C_3S1 - 3D1
       cnlo_pw(7) = -0.8023003d0 ! C_3P2
       cD =  0.8168059d0
       cE = -0.0395747d0

    else IF ( interaction == "dnnlo394go" ) then
       ! DeltaNNLOgo(394): PRC 102, 054301 (2020)
       which_pot = 3
       lambda = 394.d0
       chiral_order = NNLO
       nexp_NN = 4
       chiral_delta_flag = .TRUE.
       regulatization_scheme = chp_SFR
       power_counting = pc_epelbaum2005
       gA = 1.289d0
       hA = 1.40d0
       fpi = 92.2d0
       chp_3NF_regulator_type = chp_3NF_reg_NONLOCAL
       lambda_3NF_nonloc = lambda
       nreg_3NF_nonloc = nexp_NN
       !
       LEC_c1 = -0.7400000d0  
       LEC_c2 = -0.4900000d0
       LEC_c3 = -0.6500000d0 
       LEC_c4 =  0.9600000d0
       c1s0(-1) = -0.3381420d0 ! Ct1S0pp
       c1s0(0)  = -0.3392500d0 ! Ct1S0np
       c1s0(1)  = -0.3387460d0 ! Ct1S0nn
       c3s1(-1) = -0.2598390d0 ! Ct3S1pp
       c3s1(0)  = -0.2598390d0 ! Ct3S1np
       c3s1(1)  = -0.2598390d0 ! Ct3S1nn
       cnlo_pw(1) =  2.5053890d0 ! C_1S0 
       cnlo_pw(2) =  0.7004990d0 ! C_3P0
       cnlo_pw(3) = -0.3879600d0 ! C_1P1 
       cnlo_pw(4) = -0.9648560d0 ! C_3P1
       cnlo_pw(5) =  1.0021890d0 ! C_3S1
       cnlo_pw(6) =  0.4525230d0 ! C_3S1 - 3D1
       cnlo_pw(7) = -0.8831220d0 ! C_3P2
       cD =  0.0810000d0
       cE = -0.0020000d0

    else IF ( interaction == "dnnlo450go" ) then
       ! DeltaNNLOgo(450): PRC 102, 054301 (2020)
       which_pot = 4
       lambda = 450.d0
       chiral_order = NNLO
       nexp_NN = 3
       chiral_delta_flag = .TRUE.
       regulatization_scheme = chp_SFR
       power_counting = pc_epelbaum2005
       gA = 1.289d0
       hA = 1.40d0
       fpi = 92.2d0
       chp_3NF_regulator_type = chp_3NF_reg_NONLOCAL
       lambda_3NF_nonloc = lambda
       nreg_3NF_nonloc = nexp_NN
       !
       LEC_c1 = -0.7400000d0  
       LEC_c2 = -0.4900000d0
       LEC_c3 = -0.6500000d0 
       LEC_c4 =  0.9600000d0
       c1s0(-1) = -0.3391110d0 ! Ct1S0pp
       c1s0(0)  = -0.3401140d0 ! Ct1S0np
       c1s0(1)  = -0.3398870d0 ! Ct1S0nn
       c3s1(-1) = -0.2539500d0 ! Ct3S1pp
       c3s1(0)  = -0.2539500d0 ! Ct3S1np
       c3s1(1)  = -0.2539500d0 ! Ct3S1nn
       cnlo_pw(1) =  2.5266360d0 ! C_1S0 
       cnlo_pw(2) =  0.6719080d0 ! C_3P0
       cnlo_pw(3) = -0.2194980d0 ! C_1P1 
       cnlo_pw(4) = -0.9153980d0 ! C_3P1
       cnlo_pw(5) =  0.9649900d0 ! C_3S1
       cnlo_pw(6) =  0.4457430d0 ! C_3S1 - 3D1
       cnlo_pw(7) = -0.8954050d0 ! C_3P2
       cD = -0.4540000d0
       cE = -0.1860000d0

    else IF ( interaction == "dnnlo450" ) then
       ! Delta NNLO(450): PRC 97, 024332 (2018)
       which_pot = 5
       lambda = 450.d0
       chiral_order = NNLO
       nexp_NN = 3
       chiral_delta_flag = .TRUE.
       regulatization_scheme = chp_SFR
       power_counting = pc_epelbaum2005
       gA = 1.289d0
       hA = 1.40d0
       fpi = 92.2d0
       chp_3NF_regulator_type = chp_3NF_reg_NONLOCAL
       lambda_3NF_nonloc = lambda
       nreg_3NF_nonloc = nexp_NN
       !
       LEC_c1 = -0.7400000d0  
       LEC_c2 = -0.4900000d0
       LEC_c3 = -0.6500000d0 
       LEC_c4 =  0.9600000d0
       c1s0(-1) = -0.3371370d0 ! Ct1S0pp
       c1s0(0)  = -0.3381390d0 ! Ct1S0np
       c1s0(1)  = -0.3380230d0 ! Ct1S0nn
       c3s1(-1) = -0.2293100d0 ! Ct3S1pp
       c3s1(0)  = -0.2293100d0 ! Ct3S1np
       c3s1(1)  = -0.2293100d0 ! Ct3S1nn
       cnlo_pw(1) =  2.4765890d0 ! C_1S0 
       cnlo_pw(2) =  0.6455500d0 ! C_3P0
       cnlo_pw(3) = -0.0285410d0 ! C_1P1 
       cnlo_pw(4) = -1.0223590d0 ! C_3P1
       cnlo_pw(5) =  0.6959530d0 ! C_3S1
       cnlo_pw(6) =  0.3583300d0 ! C_3S1 - 3D1
       cnlo_pw(7) = -0.8702030d0 ! C_3P2
       cD =  0.7900000d0
       cE =  0.0170000d0

    else IF ( interaction == "dnnlo500" ) then
       ! Delta NNLO(450): PRC 97, 024332 (2018)
       which_pot = 6
       lambda = 500.d0
       chiral_order = NNLO
       nexp_NN = 3
       chiral_delta_flag = .TRUE.
       regulatization_scheme = chp_SFR
       power_counting = pc_epelbaum2005
       gA = 1.289d0
       hA = 1.40d0
       fpi = 92.2d0
       chp_3NF_regulator_type = chp_3NF_reg_NONLOCAL
       lambda_3NF_nonloc = lambda
       nreg_3NF_nonloc = nexp_NN
       !
       LEC_c1 = -0.7400000d0  
       LEC_c2 = -0.4900000d0
       LEC_c3 = -0.6500000d0 
       LEC_c4 =  0.9600000d0
       c1s0(-1) = -0.3373030d0 ! Ct1S0pp
       c1s0(0)  = -0.3383200d0 ! Ct1S0np
       c1s0(1)  = -0.3382230d0 ! Ct1S0nn
       c3s1(-1) = -0.2217210d0 ! Ct3S1pp
       c3s1(0)  = -0.2217210d0 ! Ct3S1np
       c3s1(1)  = -0.2217210d0 ! Ct3S1nn
       cnlo_pw(1) =  2.4880190d0 ! C_1S0 
       cnlo_pw(2) =  0.6984540d0 ! C_3P0
       cnlo_pw(3) = -0.0126510d0 ! C_1P1 
       cnlo_pw(4) = -0.9372640d0 ! C_3P1
       cnlo_pw(5) =  0.6753530d0 ! C_3S1
       cnlo_pw(6) =  0.3544790d0 ! C_3S1 - 3D1
       cnlo_pw(7) = -0.8595260d0 ! C_3P2
       cD = -0.8200000d0
       cE = -0.3500000d0

    else IF ( interaction == "pw_nnlo_emn500" ) then 
       ! NNLO_EMN(500): 
       which_pot = 7
       lambda = 500.d0
       chiral_order = NNLO
       nexp_NN = 4
       chiral_delta_flag = .FALSE.
       ! regulatization_scheme = chp_SFR
       ! power_counting = pc_machleidt2011
       gA = 1.29d0
       fpi = 92.4d0
       chp_3NF_regulator_type = chp_3NF_reg_NONLOCAL
       lambda_3NF_nonloc = lambda
       nreg_3NF_nonloc = nexp_NN
       !
       LEC_c1 = -0.74d0
       LEC_c2 =  0.00d0
       LEC_c3 = -3.61d0 
       LEC_c4 =  2.44d0
       c1s0(-1) = 0.d0 ! Ct1S0pp
       c1s0(0)  = 0.d0 ! Ct1S0np
       c1s0(1)  = 0.d0 ! Ct1S0nn
       c3s1(-1) = 0.d0 ! Ct3S1pp
       c3s1(0)  = 0.d0 ! Ct3S1np
       c3s1(1)  = 0.d0 ! Ct3S1nn
       cnlo_pw(1) = 0.d0 ! C_1S0 
       cnlo_pw(2) = 0.d0 ! C_3P0
       cnlo_pw(3) = 0.d0 ! C_1P1 
       cnlo_pw(4) = 0.d0 ! C_3P1
       cnlo_pw(5) = 0.d0 ! C_3S1
       cnlo_pw(6) = 0.d0 ! C_3S1 - 3D1
       cnlo_pw(7) = 0.d0 ! C_3P2
       cD = -1.50d0
       cE = -0.61d0
       
    ELSE
       ! Default = Minnesota
       which_pot = 0
       IF ( iam == 0 ) write(6,*) 'Interaction: Minnesota'
       return
    end IF
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
    ! Derived constants
    gA2 = gA*gA
    gA4 = gA2*gA2
    hA2 = hA*hA
    hA4 = hA2*hA2
    fpi2 = fpi*fpi
    fpi4 = fpi2*fpi2
    fpi_inv = 1.d0/fpi

    ! lec_c1, lec_c2, lec_c3, lec_c4 are in units of GeV^-1 (as provided in papers)
    ! c1, c2, c3, c4 are in units of MeV^-1
    c1 = LEC_c1*1.0D-3
    c2 = LEC_c2*1.0D-3
    c3 = LEC_c3*1.0D-3
    c4 = LEC_c4*1.0D-3
    LEC_c1_3NF = LEC_c1
    LEC_c2_3NF = 0.d0
    LEC_c3_3NF = LEC_c3
    LEC_c4_3NF = LEC_c4
    c1_3NF = LEC_c1_3NF*1.0D-3
    c2_3NF = LEC_c2_3NF*1.0D-3
    c3_3NF = LEC_c3_3NF*1.0D-3
    c4_3NF = LEC_c4_3NF*1.0D-3

    ! Redefine the constants so to include the Delta (TPE 3N force) [see PRC 97, 024332 (2018)]
    c3_3NF_delta = -4.d0*hA2/(9.d0*delta_split) ! in MeV^-1
    c4_3NF_delta = +2.d0*hA2/(9.d0*delta_split)
    lec_c3_3NF_delta = c3_3NF_delta * 1.e3 ! in GeV^-1
    lec_c4_3NF_delta = c4_3NF_delta * 1.e3
    IF ( chiral_delta_flag ) then
       lec_c3_3NF = lec_c3_3NF + lec_c3_3NF_delta ! in GeV^-1
       lec_c4_3NF = lec_c4_3NF + lec_c4_3NF_delta
       c3_3NF = c3_3NF + c3_3NF_delta ! in MeV^-1
       c4_3NF = C4_3NF + c4_3NF_delta
    end IF
    
    ! the NLO contacts are input in units of 10^4/GeV^4
    ! Multiply by 10^-8 to convert them to MeV^-4
    DO i=1,7
       cnlo_pw(i) = cnlo_pw(i) * 1.D-08
    END DO
    ! the LO CIB contacts are input in units of 10^4/GeV^2
    ! this = 10^-2/MeV^2 
    c1s0 = c1s0*0.01d0 
    c3s1 = c3s1*0.01d0 
    ! See Eq. 4.39 p. 26
    CS = (c1s0 + 3.d0*c3s1) /16.d0/pi
    CT = (c3s1 - c1s0) /16.d0/pi
    
    cnlo(1) = (-5.d0*cnlo_pw(7)+6.d0*cnlo_pw(5)-3.d0*cnlo_pw(4)-3.d0*cnlo_pw(3)-cnlo_pw(2)+2.d0*cnlo_pw(1))/(64.d0*pi)
    cnlo(2) = ( 5.d0*cnlo_pw(7)+6.d0*cnlo_pw(5)+3.d0*cnlo_pw(4)+3.d0*cnlo_pw(3)+cnlo_pw(2)+2.d0*cnlo_pw(1))/(16.d0*pi)
    cnlo(3) = -( 2.d0*cnlo_pw(7)-2.d0*dsqrt(2.d0)*cnlo_pw(6)-2.d0*cnlo_pw(5)-3.d0*cnlo_pw(3)+cnlo_pw(2)+2.d0*cnlo_pw(1))/(64.d0*pi)
    cnlo(4) = -(-2.d0*cnlo_pw(7)-2.d0*dsqrt(2.d0)*cnlo_pw(6)-2.d0*cnlo_pw(5)+3.d0*cnlo_pw(3)-cnlo_pw(2)+2.d0*cnlo_pw(1))/(16.d0*pi)
    cnlo(5) = -(-5.d0*cnlo_pw(7)+3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(16.d0*pi)
    cnlo(6) = ( cnlo_pw(7)-6.d0*dsqrt(2.d0)*cnlo_pw(6)-3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(64.d0*pi)
    cnlo(7) = -(cnlo_pw(7)+6.d0*dsqrt(2.d0)*cnlo_pw(6)-3.d0*cnlo_pw(4)+2.d0*cnlo_pw(2))/(16.d0*pi)

    ! NNLO 3NF constants
    ! chiral symmetry breaking scale and dimensionful cE and cD couplings
    lambdaChi = 700 ! MeV
    Econst = cE/(fpi4*lambdaChi)
    Dconst = cD/(fpi2*lambdaChi)
    
    ! This is for SFR
    DO i=-1,2
       IF ( (sfr-twompi(i)) <  0.0D0 ) sfr_heavyside(i) = 0.0D0
       IF ( (sfr-twompi(i)) >= 0.0D0 ) sfr_heavyside(i) = 1.0D0
    END DO
    z_low = 2 * mpi(2) / sfr
    IF (z_low > 1) z_low = 1

    IF ( chiral_delta_flag ) then
       mu_low  =  2.0D0*mpi(2)
       mu_high = sfr
       CALL chp_setup_delta_loop_int_data(nof_2PE_loop_int_mu_points, mu_low, mu_high)
    end IF
    
    const = 0.0D0
    const(1) = 1.0D0/(twopi**3) ! 1: 1/(2pi)^3
    const(2) = gA2/(4.0D0*fpi2) ! 2: gA^2/(4*fpi^2)
    const(3) = 384.0D0*pi2*fpi4 ! 3: 384pi^2*fpi^4
    const(4) = 5.0D0*gA4-4.0D0*gA2-1.0D0 ! 4: 5gA^4-4gA^2-1
    const(5) = 23.0D0*gA4-10.0D0*gA2 - 1.0D0 ! 5: 23gA^4-10gA^2-1
    const(6) = 48.0D0*gA4 ! 6: 48gA^4
    const(7) = 3.0D0*gA4 ! 7: 3gA^4
    const(8) = 64.0D0*pi2*fpi4 ! 8: 64pi^2fpi^4
    const(9) = 3.0D0*gA2/(16.0D0*pi*fpi4) ! 9: 3gA^2/16pifpi^4
    const(10) = ga2/16.0D0 ! 10: ga^2/16
    const(11) = 2.0D0*(2.0D0*c1-c3) ! 11: 2.0D0*(2c1-c3)
    const(12) = const(7)/(256.0D0*pi*fpi4) ! 12 : const(7)/256pifpi^4
    const(13) = gA2/(128.0D0*pi*fpi4) ! 13: gA^2/(128pifpi^4)
    const(14) = gA4/(128.0D0*pi*fpi4) ! 14: gA^4/(128pifpi^4)
    const(15) = 3.0D0*gA4/(512.0D0*pi*fpi4) ! 15: 3gA^4/(512pifpi^4)
    const(16) = gA2/(32.0D0*pi*fpi4) ! 16: gA2/(32pifpi^4)
    const(17) = gA2/8.0D0 ! 17: gA2/8
    const(18) = gA4/(256.0D0*pi*fpi4) ! 18: gA4/(256pifpi^4)
    const(19) = 3.0D0*gA4/(32.0D0*pi*fpi4) ! 19: 3gA4/(32pifpi^4)
    const(20) = const(16)*(1.0D0-gA2) ! 20: const(16)*(1-gA2)
    const(21) = 1d0/(pi2*fpi4) ! 21: 1/(pi^2*fpi^4)
    const(22) = 1D0 / (pi*fpi4) ! 22: 1 / (pi * fpi^4)
    const(23) = gA/(8.0D0*fpi2) ! 23: gA/(8*fpi^2)

    sigma_x = 0.d0
    sigma_y = 0.d0
    sigma_z = 0.d0
    sigma_x(1,2) = 1.d0
    sigma_x(2,1) = 1.d0
    sigma_y(1,2) = dcmplx(0.d0,-1.d0)
    sigma_y(2,1) = dcmplx(0.d0, 1.d0)
    sigma_z(1,1) = 1.d0
    sigma_z(2,2) = -1.d0    
    DO i=-1,2
       IF ( (sfr-twompi(i)) <  0.0D0 ) sfr_heavyside(i) = 0.0D0
       IF ( (sfr-twompi(i)) >= 0.0D0 ) sfr_heavyside(i) = 1.0D0
    END DO

    IF ( iam == 0 ) then
       write(6,*) 'Interaction:  ', interaction
       write(6,"(A12,F30.16)") 'C1', cnlo(1)* 1.D8
       write(6,"(A12,F30.16)") 'C2', cnlo(2)* 1.D8
       write(6,"(A12,F30.16)") 'C3', cnlo(3)* 1.D8
       write(6,"(A12,F30.16)") 'C4', cnlo(4)* 1.D8
       write(6,"(A12,F30.16)") 'C5', cnlo(5)* 1.D8
       write(6,"(A12,F30.16)") 'C6', cnlo(6)* 1.D8
       write(6,"(A12,F30.16)") 'C7', cnlo(7)* 1.D8
       
       write(6,*) "3N force parameters"
       write(6,"(A12,F30.8,A12)") "C1", LEC_c1_3NF, "GeV^-1"
       write(6,"(A12,F30.8,A12)") "C2", LEC_c2_3NF, "GeV^-1"
       write(6,"(A12,F30.8,A12)") "C3", LEC_c3_3NF, "GeV^-1"
       write(6,"(A12,F30.8,A12)") "C4", LEC_c4_3NF, "GeV^-1"
       write(6,"(A12,F30.5)") "cD", cD
       write(6,"(A12,F30.5)") "cE", cE
       
       IF ( chp_3NF_regulator_type.eq.chp_3NF_reg_NONLOCAL ) then
          write(6,"(A30)") "Using non-local 3N regulator"
          write(6,"(A30,F30.1)") "Lambda 3N (non-local) ", lambda_3NF_nonloc 
       elseIF ( chp_3NF_regulator_type.eq.chp_3NF_reg_LOCAL ) then 
          write(6,"(A30)") "Using local 3N regulator"
          write(6,"(A30,F30.1)") "Lambda 3N (local) ", lambda
       else
          write(6,"(A30)") "Using NO 3N regulator"
       end IF
    end IF
  
  end subroutine init_chp_constants



  SUBROUTINE chp_setup_gauleg_mesh(mesh)
    IMPLICIT none
    TYPE(chp_gauleg_mesh), INTENT(INOUT) :: mesh
    INTEGER :: i, j, m, n
    REAL*8 :: x1, x2
    REAL*8 :: p1,p2,p3,pp,xl,xm,z,z1
    REAL*8, DIMENSION(:), ALLOCATABLE :: x, w
    REAL*8, PARAMETER :: EPS = 3.D-14
    
    ! allocate points and weights storages
    ALLOCATE(x(mesh%info%amount))
    ALLOCATE(w(mesh%info%amount))
    IF (ALLOCATED(mesh%pnt) ) DEALLOCATE (mesh%pnt)
    IF (mesh%info%amount <= 0 .or. mesh%info%xmin > mesh%info%xmax) then
       WRITE(*,*) ': incorrect mesh info', mesh%info; stop
    end IF
    ALLOCATE( mesh%pnt( 1:mesh%info%amount ) )
    mesh%pnt(:) = chp_gauleg_mesh_point(0.0D0, 0.0D0, 0.0D0)
    
    ! set values of local variables
    x1 = mesh%info%xmin; x2 = mesh%info%xmax; n = mesh%info%amount
    m=(n+1)/2
    xm=0.5D0*(x2+x1)
    xl=0.5D0*(x2-x1)
    DO i=1,m
       z1=0.0D0
       z=COS(pi*(i - 0.25D0)/(n + 0.5D0))
       pp=n*z/(z*z-1.0D0)
       DO while ( abs(z-z1) > EPS )
          p1=1.0D0
          p2=0.0D0
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2.0D0*j-1.0D0)*z*p2-(j-1.0D0)*p3)/j
          end DO
          pp=n*(z*p1-p2)/(z*z-1.0D0)
          z1=z
          z=z-p1/pp
       end DO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.0D0*xl/((1.0D0-z*z)*pp*pp)
       w(n+1-i)=w(i)
    end DO
    ! set return values
    mesh%pnt(:)%x = x(:); mesh%pnt(:)%w = w(:); mesh%pnt(:)%xxw = x(:) * x(:) * w(:)
    DEALLOCATE(w)
    DEALLOCATE(x)
  END SUBROUTINE chp_setup_gauleg_mesh
  
  ! The expressions for some of the loop diagrams contains integrals that must be solved numerically
  ! This function sets up the integration meshs used in the numerical integrations
  SUBROUTINE chp_setup_delta_loop_int_data(nof_mu_points, mu_low_in, mu_high_in)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: nof_mu_points
    REAL*8, INTENT(IN) :: mu_high_in, mu_low_in
    chp_delta_2PE_loop_int_mesh_mu%info = chp_mesh_info(nof_mu_points, mu_low_in, mu_high_in)
    CALL chp_setup_gauleg_mesh(chp_delta_2PE_loop_int_mesh_mu)
    if ( nof_mu_points == 25 ) then 
       chp_delta_2PE_loop_int_mesh_mu%pnt(13)%x = 488.038986666667
       chp_delta_2PE_loop_int_mesh_mu%pnt(13)%w = 25.5194921127220
    end if
  END SUBROUTINE chp_setup_delta_loop_int_data
  
  
  
  FUNCTION chp_sigma_dot_sigma(ms1,ms2,ms3,ms4) result(res)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    COMPLEX*16 :: res, res1, spin, jx1, jx2, jy1, jy2, jz1, jz2 
    REAL*8 :: m1,m2,m3,m4
    
    m1 = 0.5d0*dble(ms1); m2 = 0.5d0*dble(ms2)
    m3 = 0.5d0*dble(ms3); m4 = 0.5d0*dble(ms4)    
    jx1 = 0.d0
    jx2 = 0.d0 
    jy1 = 0.d0
    jy2 = 0.d0 
    jz1 = 0.d0
    jz2 = 0.d0 
    
    spin = 0.5d0 
    res = 0.0D0
    IF ( ms1 == ms3 ) THEN
       jx1 = 0.d0
       jy1 = 0.d0 
       jz1 = 2.d0*m1
    end if
    IF ( ms1 == ms3 + 2 ) then 
       jx1 = sqrt( (spin+m3+1.d0)*(spin-m3) )
       jy1 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m3+1.d0)*(spin-m3) )
       jz1 = 0.d0
    end if
    IF ( ms1 == ms3 - 2 ) then 
       jx1 = sqrt( (spin-m3+1.d0)*(spin+m3) )
       jy1 = dcmplx(0.d0,1.d0)*sqrt( (spin-m3+1.d0)*(spin+m3) )
       jz1 = 0.d0
    end if
    
    IF ( ms2 == ms4 ) THEN
       jx2 = 0.d0
       jy2 = 0.d0 
       jz2 = 2.d0*m2 
    end if
    IF ( ms2 == ms4 + 2 ) then 
       jx2 = sqrt( (spin+m4+1.d0)*(spin-m4) )
       jy2 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m4+1.d0)*(spin-m4) )
       jz2 = 0.d0
    end if
    IF ( ms2 == ms4 - 2 ) then 
       jx2 = sqrt( (spin-m4+1.d0)*(spin+m4) )
       jy2 = dcmplx(0.d0,1.d0)*sqrt( (spin-m4+1.d0)*(spin+m4) )
       jz2 = 0.d0
    end if
    
    res1 = jx1*jx2 + jy1*jy2 + jz1*jz2
    res = res1
    
  end FUNCTION chp_sigma_dot_sigma

  FUNCTION chp_tau_dot_tau(ms1,ms2,ms3,ms4) result(res)    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    COMPLEX*16 :: res, res1, spin, jx1, jx2, jy1, jy2, jz1, jz2 
    REAL*8 :: m1,m2,m3,m4
    
    m1 = 0.5d0*dble(ms1); m2 = 0.5d0*dble(ms2)
    m3 = 0.5d0*dble(ms3); m4 = 0.5d0*dble(ms4)
    jx1 = 0.d0
    jx2 = 0.d0 
    jy1 = 0.d0
    jy2 = 0.d0 
    jz1 = 0.d0
    jz2 = 0.d0 
    
    spin = 0.5d0 
    res = 0.0D0
    IF ( ms1 == ms3 ) THEN
       jx1 = 0.d0
       jy1 = 0.d0 
       jz1 = 2.d0*m1
    end if
    IF ( ms1 == ms3 + 2 ) then 
       jx1 = sqrt( (spin+m3+1.d0)*(spin-m3) )
       jy1 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m3+1.d0)*(spin-m3) )
       jz1 = 0.d0
    end if
    IF ( ms1 == ms3 - 2 ) then 
       jx1 = sqrt( (spin-m3+1.d0)*(spin+m3) )
       jy1 = dcmplx(0.d0,1.d0)*sqrt( (spin-m3+1.d0)*(spin+m3) )
       jz1 = 0.d0
    end if
    
    IF ( ms2 == ms4 ) THEN
       jx2 = 0.d0
       jy2 = 0.d0 
       jz2 = 2.d0*m2 
    end if
    IF ( ms2 == ms4 + 2 ) then 
       jx2 = sqrt( (spin+m4+1.d0)*(spin-m4) )
       jy2 = -dcmplx(0.d0,1.d0)*sqrt( (spin+m4+1.d0)*(spin-m4) )
       jz2 = 0.d0
    end if
    IF ( ms2 == ms4 - 2 ) then 
       jx2 = sqrt( (spin-m4+1.d0)*(spin+m4) )
       jy2 = dcmplx(0.d0,1.d0)*sqrt( (spin-m4+1.d0)*(spin+m4) )
       jz2 = 0.d0
    end if
            
    res1 = jx1*jx2 + jy1*jy2 + jz1*jz2
    res = res1
    
  end FUNCTION chp_tau_dot_tau

  FUNCTION chp_sigma_dot_q_mtx(ms1,ms3,q) result(res)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ms1, ms3
    REAL*8, INTENT(IN) :: q(3)
    COMPLEX*16 :: res, res1
    COMPLEX*16 :: chi1(2), chi3(2), mat(2,2)
    INTEGER :: i1
    
    chi1 = 0.d0; chi3 = 0.d0;
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    
    mat = q(1)*sigma_x+q(2)*sigma_y+q(3)*sigma_z
    res1 = dot_product(chi1, matmul( mat, chi3))
    res = res1
  end FUNCTION chp_sigma_dot_q_mtx

  FUNCTION chp_tau_dot_tau_mtx(ms1,ms2,ms3,ms4) result(res)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1,ms2,ms3,ms4
    COMPLEX*16 :: res, res1
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2)
    INTEGER :: i1 
    
    res = 0.0D0
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms4) 
    chi4(i1) = 1.d0 
    
    res1 = dot_product(chi1, matmul( sigma_x, chi3)) * dot_product(chi2, matmul( sigma_x, chi4)) &
         + dot_product(chi1, matmul( sigma_y, chi3)) * dot_product(chi2, matmul( sigma_y, chi4)) &
         + dot_product(chi1, matmul( sigma_z, chi3)) * dot_product(chi2, matmul( sigma_z, chi4)) 
    res = res1    
  end FUNCTION chp_tau_dot_tau_mtx
  
  FUNCTION chp_sigma1_dot_q_sigma2_dot_q_mtx(ms1,ms2,ms3,ms4,q) result(res)    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ms1, ms2, ms3, ms4
    REAL*8, INTENT(IN) :: q(3)
    COMPLEX*16 :: res, res1, res2
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2), mat(2,2)
    INTEGER :: i1
    
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0;
    i1 = nint(1.5-0.5*ms1)
    chi1(i1) = 1.d0
    i1 = nint(1.5-0.5*ms2)
    chi2(i1) = 1.d0
    i1 = nint(1.5-0.5*ms3)
    chi3(i1) = 1.d0
    i1 = nint(1.5-0.5*ms4)
    chi4(i1) = 1.d0
    
    mat = q(1)*sigma_x+q(2)*sigma_y+q(3)*sigma_z
    res1 = dot_product(chi1, matmul( mat, chi3))
    res2 = dot_product(chi2, matmul( mat, chi4))
    res = res1 * res2
  end FUNCTION chp_sigma1_dot_q_sigma2_dot_q_mtx

  FUNCTION chp_spin_dot_qxk_mtx(ms1,ms2,ms3,ms4,qxk) result(res)
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: ms1, ms2, ms3, ms4
    REAL*8, INTENT(IN) :: qxk(3)
    COMPLEX*16 :: res, res1, res2 
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2), mat(2,2) 
    REAL*8 :: delta 
    INTEGER :: i1 
    
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms4) 
    chi4(i1) = 1.d0 
    
    mat =-dcmplx(0.d0,1.d0) * 0.5D0*(qxk(1)*sigma_x+qxk(2)*sigma_y+qxk(3)*sigma_z)    
    res1 = dot_product(chi1, matmul( mat, chi3))*delta(ms2,ms4)
    res2 = dot_product(chi2, matmul( mat, chi4))*delta(ms1,ms3)
    res = (res1+res2)    
  end FUNCTION chp_spin_dot_qxk_mtx

  FUNCTION chp_tau1_dot_tauXtau(t1,t2,t3,t4,t5,t6) result(res)    
    IMPLICIT NONE 
    INTEGER, INTENT(IN) :: t1, t2, t3, t4, t5, t6
    COMPLEX*16 :: facx, facy, facz,res
    COMPLEX*16 :: chi1(2), chi2(2), chi3(2), chi4(2), chi5(2), chi6(2), mat(2,2)
    INTEGER :: i1
    
    chi1 = 0.d0; chi2 = 0.d0; chi3 = 0.d0; chi4 = 0.d0; chi5 = 0.d0; chi6 = 0.d0
    i1 = nint(1.5-0.5*t1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*t2) 
    chi2(i1) = 1.d0 
    i1 = nint(1.5-0.5*t3) 
    chi3(i1) = 1.d0 
    i1 = nint(1.5-0.5*t4) 
    chi4(i1) = 1.d0 
    i1 = nint(1.5-0.5*t5) 
    chi5(i1) = 1.d0 
    i1 = nint(1.5-0.5*t6) 
    chi6(i1) = 1.d0 
    
    facx = dot_product(chi2,matmul(sigma_y,chi5))*dot_product(chi3,matmul(sigma_z,chi6))  - &
         dot_product(chi2,matmul(sigma_z,chi5))*dot_product(chi3,matmul(sigma_y,chi6))  

    facy = dot_product(chi2,matmul(sigma_z,chi5))*dot_product(chi3,matmul(sigma_x,chi6))  - &
         dot_product(chi2,matmul(sigma_x,chi5))*dot_product(chi3,matmul(sigma_z,chi6))  

    facz = dot_product(chi2,matmul(sigma_x,chi5))*dot_product(chi3,matmul(sigma_y,chi6))  - &
         dot_product(chi2,matmul(sigma_y,chi5))*dot_product(chi3,matmul(sigma_x,chi6))  
    
    mat = facx*sigma_x+facy*sigma_y+facz*sigma_z 
    res = dot_product(chi1, matmul( mat, chi4))    
  end FUNCTION chp_tau1_dot_tauXtau

  
  FUNCTION chp_sigma1_dot_q_ope(ms1,ms2,q) result(res)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ms1, ms2
    REAL*8, INTENT(IN) :: q(3)
    REAL*8 :: q2
    COMPLEX*16 :: res, res1 
    COMPLEX*16 :: chi1(2), chi2(2), mat(2,2) 
    INTEGER :: i1 
    
    chi1 = 0.d0; chi2 = 0.d0
    i1 = nint(1.5-0.5*ms1) 
    chi1(i1) = 1.d0 
    i1 = nint(1.5-0.5*ms2) 
    chi2(i1) = 1.d0 
    
    mat = q(1)*sigma_x+q(2)*sigma_y+q(3)*sigma_z 
    res1 = dot_product(chi1, matmul( mat, chi2))    
    q2 = sum( q*q )
    res = res1 /( q2 + mpi2(2))
  end FUNCTION chp_sigma1_dot_q_ope
  
  
  
  ! NLO loop FUNCTION w [1] Eq 4.12 (DR and SFR)
  FUNCTION chp_NLO_two_pion_exchange_loop_w(q2, impi) result(res)
    IMPLICIT none
    REAL*8, INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = DSQRT(fourmpi2(impi) + q2)
  end FUNCTION chp_NLO_two_pion_exchange_loop_w
  
  ! NLO SFR loop FUNCTION s [2] Eq 2.16
  FUNCTION chp_NLO_sfr_two_pion_exchange_loop_s(impi) result(res)
    IMPLICIT none
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = DSQRT(sfr2 - fourmpi2(impi))
  end FUNCTION chp_NLO_sfr_two_pion_exchange_loop_s
  
  ! NLO SFR loop FUNCTION L [2] Eq 2.16
  ! Note: s is a constant, while w is FUNCTION of q
  FUNCTION chp_NLO_sfr_two_pion_exchange_loop_L(q, q2, w, s, impi) result(res)
    IMPLICIT none
    ! REAL*8, INTENT(INOUT) :: q, q2
    REAL*8, INTENT(IN) :: q, q2
    REAL*8, INTENT(IN) :: w, s
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res, eps, q0, q20
    
    res = 0.D0
    eps = 1.D-8 
    IF ( sfr_heavyside(impi) == 0.D0 ) return
    q0 = q
    q20 = q2
    IF ( q == 0.d0 ) then  
       q0 = q0 + eps 
       q20 = q0*q0 
    end IF
    ! res = w * dlog( (sfr2*w*w+q2*s*s+2.0D0*sfr*q*w*s)/( fourmpi2(impi)*(sfr2+q2) ) )/(2.0D0*q)
    res = w * dlog( (sfr2*w*w+q20*s*s+2.0D0*sfr*q0*w*s)/( fourmpi2(impi)*(sfr2+q20) ) )/(2.0D0*q0)    
  end FUNCTION chp_NLO_sfr_two_pion_exchange_loop_L

  ! NNLO loop FUNCTION wtilde SQUARED [1] Eq 4.20 (DR)
  FUNCTION chp_NNLO_two_pion_exchange_loop_wtilde2(q2, impi) result(res)
    IMPLICIT none
    REAL*8, INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = 2.0D0*mpi2(impi) + q2    
  end FUNCTION chp_NNLO_two_pion_exchange_loop_wtilde2
  
  ! NNLO loop FUNCTION wtilde [2] Eq 2.17 (SFR), [6] Eq. 49
  FUNCTION chp_NNLO_sfr_two_pion_exchange_loop_A(q, q2, impi) result(res)
    IMPLICIT none    
    ! REAL*8, INTENT(INOUT) :: q, q2
    REAL*8, INTENT(IN) :: q, q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res, eps, q0, q20
    
    res = 0.0D0
    eps = 1.0D-8
    IF ( sfr_heavyside(impi) == 0.D0 ) return
    q0 = q
    q20 = q2
    IF ( q == 0.d0 ) then
       q0 = q0 + eps
       q20 = q0*q0
    end IF
    ! res = datan( q*(sfr-twompi(impi) )/(q2 + sfr*twompi(impi) ) )/(2.0D0*q)
    res = datan( q0*(sfr-twompi(impi) )/(q20 + sfr*twompi(impi) ) )/(2.0D0*q0)
  end FUNCTION chp_NNLO_sfr_two_pion_exchange_loop_A


  
  ! NNLO FUNCTION [1] Eq 4.13 (ci part) OR [2] Eq 2.15 (V_C part) OR [3] Eq C1
  FUNCTION chp_two_pion_exchange_1loop_d_Vc(q2, A, wt2, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2, A, wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res    
    res = -const(9)*(mpi2(impi)*const(11) - q2*c3)*wt2*A
  end FUNCTION chp_two_pion_exchange_1loop_d_Vc

  ! NNLO FUNCTION Eq [1] 4.16
  FUNCTION chp_two_pion_exchange_1loop_d_WT(w, A) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: w, A
    REAL*8 :: res
    res = -1.0D0*const(16)*A*( c4*w*w )    
  end FUNCTION chp_two_pion_exchange_1loop_d_WT

  ! NNLO FUNCTION Eq [1] 4.16
  FUNCTION chp_two_pion_exchange_1loop_d_Ws(q2, w, A) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2, w, A
    REAL*8 :: res    
    res = -q2 * chp_two_pion_exchange_1loop_d_WT(w, A)
  end FUNCTION chp_two_pion_exchange_1loop_d_Ws
  

  
  ! static one pion exchange, [1] eq 4.5
  FUNCTION chp_one_pion_exchange(q2, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = -1.0D0 * const(2) / (q2 + mpi2(impi))
  end FUNCTION chp_one_pion_exchange  

  
  ! NLO FUNCTION [1] 4.9 OR [2] 2.14 (W_C part) OR [3] B1 OR [4] 6
  FUNCTION chp_NLO_two_pion_exchange_Wc(q2, L, w, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2, L, w
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = -L *(fourmpi2(impi)*const(4) + q2*const(5) + const(6)*mpi4(impi)/(w*w))/const(3) 
  end FUNCTION chp_NLO_two_pion_exchange_Wc
  
  ! NLO FUNCTION Eq 4.10
  FUNCTION chp_NLO_two_pion_exchange_Vs(q2,L,impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2, L
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = const(7)*L*q2/const(8)    
  end FUNCTION chp_NLO_two_pion_exchange_Vs
  
  ! NLO Vt FUNCTION Eq 4.10
  FUNCTION chp_NLO_two_pion_exchange_VT(L,impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: L
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = -const(7)*L/const(8)    
  end FUNCTION chp_NLO_two_pion_exchange_VT


  ! NNLO FUNCTION Eq [1] 4.13 and 4.21
  FUNCTION chp_NNLO_two_pion_exchange_Vc(q2, w, A, wt2, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2, w, A, wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = 0.0D0
    IF ( power_counting == pc_machleidt2011 ) then
       res = const(9)*(const(10)*mpi5(impi)/(mnuc(0)*w*w) - &
            (mpi2(impi)*const(11) - q2*c3 - q2*3.0D0*const(10)/mnuc(0) ) *wt2*A)
       ! IF 2PE in EM format add correction to the NNLO central term
       IF ( chp_itope == 'EM' ) THEN
          res = res - const(12)*(mpi(impi)*w*w+wt2*wt2*A)/mnuc(0)
       end IF
    else IF ( power_counting == pc_epelbaum2005 ) then
       !from 1loop_d_Vc
       res = const(9)*(-(mpi2(impi)*const(11) - q2*c3) *wt2*A)
       !from 1loop_r_Vc
       res = res + 3*const(14)*(mpi5(impi)/(2*w*w) + (2*mpi2(impi) + q2) * (q2 - mpi2(impi))*A) / mnuc(0)
    end IF
  end FUNCTION chp_NNLO_two_pion_exchange_Vc

  ! NNLO FUNCTION Eq [1] 4.14 and 4.22
  FUNCTION chp_NNLO_two_pion_exchange_Wc(q2, w, A, wt2, impi) result(res)
    IMPLICIT NONE    
    REAL*8, INTENT(IN) :: q2, w, A, wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = 0.0D0
    IF ( power_counting == pc_machleidt2011 ) then ! [1] 4.14      
       res = const(13) * (3.0D0*gA2*mpi5(impi)/(w*w) - & 
            (fourmpi2(impi) + 2.0D0*q2 - gA2*(fourmpi2(impi)+3.0D0*q2))*wt2*A)/mnuc(0)
       ! IF 2PE in EM format add correction to the NNLO central term
       IF ( chp_itope == 'EM' ) THEN   ! [1] 4.22
          res = res + const(14)*(mpi(impi)*w*w + wt2*wt2*A)/mnuc(0)
       end IF
    else IF ( power_counting == pc_epelbaum2005 ) then
       !from loop_r_Wc 
       res = 2.0D0*const(13) * (3.0D0*gA2*mpi5(impi)/(2.0D0*w*w) + &
            (gA2*(3.0D0*mpi2(impi) + 2.0D0*q2) - 2.0D0*mpi2(impi) - q2)*&
            (2.0D0*mpi2(impi) + q2)*A)/mnuc(0)
    end IF
  end FUNCTION chp_NNLO_two_pion_exchange_Wc
  
  ! NNLO FUNCTION Eq [1] 4.15 and 4.23
  FUNCTION chp_NNLO_two_pion_exchange_VT(w, A, wt2, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: w, A, wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = 0.0D0
    IF ( power_counting == pc_machleidt2011 ) then
       res = 3.0D0*const(15)*wt2*A/mnuc(0)
       ! IF 2PE in EM format add correction to the NNLO central term
       IF ( chp_itope == 'EM' ) THEN
          res = res + const(15)*(mpi(impi) + w*w*A )/mnuc(0)
       end IF
       ! else IF ( power_counting == pc_epelbaum2005 ) then
       !    ! else IF ( power_counting == 1 ) then ! sam only used in pc_machleidt2011
       !    ![3] D9 = [4] 19 + 22
       !    res = 3.0D0*const(18)*(5.0D0*mpi2(impi) + 2.0D0*q2)*A/mnuc(0)
    end IF
  end FUNCTION chp_NNLO_two_pion_exchange_VT

  ! NNLO FUNCTION Eq [1] 4.15 and 4.23
  FUNCTION chp_NNLO_two_pion_exchange_Vs(q2, w, A, wt2, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2, w, A, wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = 0.0D0
    IF ( power_counting == pc_machleidt2011 ) then
       res = -3.0D0*q2*const(15)*wt2*A/mnuc(0)
       ! IF 2PE in EM format add correction to the NNLO central term
       IF (chp_itope == 'EM') THEN
          res = res - 1.0D0*const(15)*q2*(mpi(impi) + w*w*A)/mnuc(0)
       end IF
    else IF ( power_counting == pc_epelbaum2005 ) then
       res = -q2*3.0D0*const(18)*(5.0D0*mpi2(impi) + 2.0D0*q2)*A/mnuc(0)
    end IF    
  end FUNCTION chp_NNLO_two_pion_exchange_Vs
  
  ! NNLO FUNCTION Eq [1] 4.16 and 4.24
  FUNCTION chp_NNLO_two_pion_exchange_WT(q2, w, A, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2, w, A
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = 0.0D0
    IF ( power_counting == pc_machleidt2011 ) then
       res = -1.0D0*const(16)*A*( (c4 + 1.0D0/(4.0D0*mnuc(0)))*w*w - &
            const(17)*(10.0D0*mpi2(impi) + 3.0D0*q2)/mnuc(0))
       ! IF 2PE in EM format add correction to the NNLO central term
       IF (chp_itope == 'EM') THEN
          res = res - const(18)*(mpi(impi) + w*w*A)/mnuc(0)
       end IF
    else IF ( power_counting == pc_epelbaum2005 ) then
       res =  -1.0D0*const(16)*A*( c4*w*w )
       ! [3] D10 = [4] 19 + 22
       res = res + const(13)*(gA2*(3.0D0*mpi2(impi) + q2) - w**2)*A/mnuc(0)
    end IF
  end FUNCTION chp_NNLO_two_pion_exchange_WT

  ! NNLO FUNCTION Eq [1] 4.16 and 4.24
  FUNCTION chp_NNLO_two_pion_exchange_Ws(q2, w, A, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2, w, A
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = 0.0D0
    IF ( power_counting == pc_machleidt2011 ) then
       res = q2*const(16)*A*( (c4+1.0D0/(4.0D0*mnuc(0)))*w*w - & 
            const(17)*(10.0D0*mpi2(impi) +3.0D0*q2)/mnuc(0))
       ! IF 2PE in EM format add correction to the NNLO central term
       IF (chp_itope == 'EM') THEN
          res = res + const(18)*q2*(mpi(impi) + w*w*A)/mnuc(0)
       end IF
    else IF ( power_counting == pc_epelbaum2005 ) then
       res = q2*1.0D0*const(16)*A*( c4*w*w )
       ! [3] D10 = [4] 19 + 22
       res = res - q2*const(13)*(gA2*(3.0D0*mpi2(impi) + q2) - w**2)*A/mnuc(0)
    end IF
  end FUNCTION chp_NNLO_two_pion_exchange_Ws

  ! NNLO FUNCTION Eq [1] 4.17
  FUNCTION chp_NNLO_two_pion_exchange_VLS(A, wt2, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: A, wt2
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = const(19) * wt2*A/mnuc(0)             
  end FUNCTION chp_NNLO_two_pion_exchange_VLS

  ! NNLO FUNCTION Eq [1] 4.18
  FUNCTION chp_NNLO_two_pion_exchange_WLS(w, A, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: w, A
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = const(20)*w*w*A/mnuc(0)             
  end FUNCTION chp_NNLO_two_pion_exchange_WLS

  

  ! Dloop FUNCTION [5] Eq 2.8
  FUNCTION chp_delta_2PE_sfr_loop_D(q2, delta_in, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2, delta_in
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: integral, res
    REAL*8 :: mu, mu2, mut
    INTEGER :: j, mu_nr
    
    mu_nr = chp_delta_2PE_loop_int_mesh_mu%info%amount
    integral = 0.0D0
    DO j = 1, mu_nr
       mu = chp_delta_2PE_loop_int_mesh_mu%pnt(j)%x
       mu2 = mu*mu
       mut = sqrt( mu2 - fourmpi2(impi) )/(2.0D0*delta_in)
       integral = integral + chp_delta_2PE_loop_int_mesh_mu%pnt(j)%w * ( atan(mut)/ (mu2 + q2) )
    END DO
    res = integral/delta_in
    
  end FUNCTION chp_delta_2PE_sfr_loop_D
  
  ! H loop functions [5] Eq 2.8
  FUNCTION chp_delta_loop_H(q2, delta, w, sigma, L, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2, delta, w, sigma, L
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res, qdel, qdel2, delL, wdel, sdel
    
    qdel  = 2.0D0*sqrt(delta**2 - mpi2(impi))
    qdel2 = qdel*qdel
    wdel  = sqrt(qdel2 + fourmpi2(impi))
    sdel  = sqrt(sfr2 - fourmpi2(impi))
    delL  = 0.0D0
    IF ( sfr_heavyside(impi) == 1.0D0 ) then    
       delL = wdel * dlog( (sfr2*wdel*wdel+qdel2*sdel*sdel+2.0D0*sfr*qdel*wdel*sdel)/( fourmpi2(impi)*(sfr2+qdel2) ) )/(2.0D0*qdel)
    end IF
    res = 2.0D0*sigma*(L - delL)/(w**2-4.0D0*delta**2)
    
  end FUNCTION chp_delta_loop_H

  
  ! Vc contributions
  FUNCTION chp_NLO_delta_2PE_single_box_Vc_SFR(q2, A, delta_in, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2, A, delta_in
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = -gA2*hA2*const(22)*( 2.0D0*mpi2(impi) + q2 )*( 2.0D0*mpi2(impi) + q2 )*A/(12.0D0*delta_in)
  end FUNCTION chp_NLO_delta_2PE_single_box_Vc_SFR
  
  FUNCTION chp_NLO_delta_2PE_double_box_Vc_SFR(L, S, D, H, delta_in, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: L, S, D, H, delta_in
    INTEGER, INTENT(IN) :: impi    
    REAL*8 :: res
    res = -hA4*const(21)*(-4.0D0*delta_in*delta_in*L + S*(H + (S + 8.0D0*delta_in*delta_in)*D))/(27.0D0)    
  end FUNCTION chp_NLO_delta_2PE_double_box_Vc_SFR

  FUNCTION chp_NNLO_delta_2PE_triangle_Vc_SFR(w, L, S, D, delta_in, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: w, L, S, D, delta_in
    INTEGER, INTENT(IN) :: impi    
    REAL*8 :: res
    res = -hA2*const(21)*delta_in*( 6.0D0*S*(4.0D0*c1*mpi2(impi)-2.0D0*c2*delta_in*delta_in - &
         c3*(2.0D0*delta_in*delta_in + S))*D + (-24.0D0*c1*mpi2(impi) + c2*(w*w-6.0D0*S) + &
         6.0D0*c3*(2.0D0*delta_in*delta_in+S))*L)/(18.0D0)
  end FUNCTION chp_NNLO_delta_2PE_triangle_Vc_SFR


  ! Wc contributions
  FUNCTION chp_NLO_delta_2PE_single_box_Wc_SFR(q2, L, S, D, delta_in, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2, L, S, D, delta_in
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = -gA2*hA2*const(21)*( (12.0D0*delta_in*delta_in - 20.0D0*mpi2(impi) - 11.0D0*q2 )*L + 6.0D0*S*S*D )/(216.0D0)
  end FUNCTION chp_NLO_delta_2PE_single_box_Wc_SFR
  
  FUNCTION chp_NLO_delta_2PE_double_box_Wc_SFR(w, L, S, D, H, delta_in, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: w, L, S, D, H, delta_in
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = -hA4*const(21)*((12.0D0*S-w*w)*L + 3.0D0*S*(H + (8.0*delta_in*delta_in-S)*D))/(486.0D0)    
  end FUNCTION chp_NLO_delta_2PE_double_box_Wc_SFR
  
  FUNCTION chp_NLO_delta_2PE_single_triangle_Wc_SFR(w, L, S, D, delta_in) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: w, L, S, D, delta_in
    REAL*8 :: res
    res = -hA2*const(21)*( (6.0D0*S-w*w)*L + 12.0D0*delta_in*delta_in*S*D)/216.0D0
  end FUNCTION chp_NLO_delta_2PE_single_triangle_Wc_SFR

  
  ! VT contributions
  FUNCTION chp_NLO_delta_2PE_single_box_Vt_SFR(w, L, D, delta_in, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: w, L, D, delta_in
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = -gA2*hA2*const(21)*(-2.0D0*L + (w*w - 4.0D0*delta_in**2)*D )/(48.0D0)    
  end FUNCTION chp_NLO_delta_2PE_single_box_Vt_SFR
  
  FUNCTION chp_NLO_delta_2PE_double_box_Vt_SFR(w, L, D, delta_in, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: w, L, D, delta_in
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = -hA4*const(21)*(6.0D0*L + (12.0D0*delta_in*delta_in - w*w)*D)/(216.0D0)    
  end FUNCTION chp_NLO_delta_2PE_double_box_Vt_SFR
  

  ! WT contributions
  FUNCTION chp_NLO_delta_2PE_single_box_Wt_SFR(w, A, delta_in, impi) result(res)
    REAL*8, INTENT(IN) :: w, A, delta_in
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = -gA2*hA2*const(22)*(w*w*A )/(144.0D0*delta_in)    
  end FUNCTION chp_NLO_delta_2PE_single_box_Wt_SFR
  
  FUNCTION chp_NLO_delta_2PE_double_box_Wt_SFR(w, L, D, delta_in, impi) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: w, L, D, delta_in
    INTEGER, INTENT(IN) :: impi
    REAL*8 :: res
    res = -hA4*const(21)*(2.0D0*L + (4.0D0*delta_in*delta_in + w*w)*D)/(1296.0D0)    
  end FUNCTION chp_NLO_delta_2PE_double_box_Wt_SFR

  FUNCTION chp_NNLO_delta_2PE_triangle_Wt_SFR(w, L, D, delta_in) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: w, L, D, delta_in    
    REAL*8 :: res
    res = -c4*hA2*const(21)*delta_in*( (w*w-4.0D0*delta_in*delta_in)*D-2.0D0*L )/(72.0D0)    
  end FUNCTION chp_NNLO_delta_2PE_triangle_Wt_SFR

  

  ! regulator, eq 4.63 
  ! pfinal : final momentum
  ! pinit  : initial momentum
  ! n      : cutoff order 
  FUNCTION freg(pfinal, pinit, n) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: pfinal(3), pinit(3)
    INTEGER, INTENT(IN) :: n
    REAL*8 :: res,exponent, p2, pp2 
    ! res = 1.e-4
    ! return
    p2 = sum( pfinal*pfinal )
    pp2 = sum( pinit*pinit )
    exponent = (p2**n/lambda**(2*n) + pp2**n/lambda**(2*n) )
    res = dexp(-exponent)
  end FUNCTION freg
  
  ! local 3NF regulator, eq 11 in P. Navratil's paper
  ! pfinal : final momentum
  ! pinit  : initial momentum
  ! n      : cutoff order 
  FUNCTION freg_3NFlocal(q2) result(res)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: q2
    REAL*8 :: res,exponent
    ! res = 1.e-4
    ! return
    res = 1.d0   
    IF ( chp_3NF_regulator_type .ne. chp_3NF_reg_LOCAL ) return
    exponent = (q2/lambda**2)**2
    res = dexp(-exponent)
  end FUNCTION freg_3NFlocal
  
  ! non-local 3NF regulator, eq 3 in http://arxiv.org/abs/1304.2212
  ! pfinal : final momentum
  ! pinit  : initial momentum
  ! n      : cutoff order 
  ! LAMBDA is accessed from the chp constant chp_lambda
  FUNCTION freg_3NFnonlocal(k1,k2,k3) result(res)
    USE constants
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: k1(3), k2(3), k3(3)
    REAL*8 :: p(3), q(3) 
    REAL*8 :: res,exponent, q2, llambda
    ! res = 1.e-4
    ! return
    res = 1.d0 
    IF ( chp_3NF_regulator_type .ne. chp_3NF_reg_NONLOCAL ) return ! if local 3NF is used, return 1
    p = 0.5d0*(k1-k2)
    q = (2./3.)*(k3-0.5d0*(k1+k2))
    llambda = 4.d0 * lambda_3NF_nonloc**2 
    q2 = 4.d0*dot_product(p,p)+ 3.d0*dot_product(q,q)
    exponent = (q2/llambda)**nreg_3NF_nonloc
    res = dexp(-exponent)
  end FUNCTION freg_3NFnonlocal
  
end MODULE chiral_constants


MODULE chiral_potentials
  IMPLICIT none
contains 
  
  FUNCTION chiral_NN(p,q,r,s) result(res) 
    USE single_particle_orbits
    USE constants
    USE chiral_constants
    USE chiral_tables
    USE ang_mom_functions, only : tjs 
    
    IMPLICIT none 
    INTEGER, INTENT(IN) :: p,q,r,s
    COMPLEX*16 :: res
    REAL*8 :: k1(3), k2(3), k3(3), k4(3)
    REAL*8 :: qtrans(3), prel(3), pprel(3), qxk(3), kmean(3) 
    INTEGER :: nx1,ny1,nz1, nx2,ny2,nz2, nx3,ny3,nz3, nx4,ny4,nz4
    INTEGER :: m1,m2,m3,m4, t1,t2,t3,t4, Tiso, iph
    REAL*8 :: nucleon_mass, relativity_factor, freg_nnlo
    REAL*8 :: q2, p2, kmean2, qabs, pp2, cg1, cg2
    REAL*8 :: delta_spin, delta_isospin
    COMPLEX*16 :: sigma1_dot_q_sigma2_dot_qtrans
    COMPLEX*16 :: t_dot_t, s_dot_s, LS_struct
    REAL*8 :: loop_A, loop_D, loop_H, loop_L
    REAL*8 :: loop_s, loop_Sigma, loop_w, loop_wtilde2
    COMPLEX*16 :: term_CS, term_CT, term_C(1:7)
    COMPLEX*16 :: Vc, Wc, Vs, Ws, VLS, WLS, VT, WT
    COMPLEX*16 :: Vc_delta, Wc_delta, Vt_delta, Wt_delta, Vs_delta, Ws_delta    
    COMPLEX*16 :: vlo, vnlo, vnnlo, vdir
    COMPLEX*16 :: cont_lo, cont_nlo
   
    res = 0.d0
    
    ! Conservation of isospin
    IF ( all_orbit%tz(p) + all_orbit%tz(q) /= all_orbit%tz(r) + all_orbit%tz(s) ) return
    
    Vc      = 0.0D0
    Wc      = 0.0D0
    Vs      = 0.0D0
    Ws      = 0.0D0
    VLS     = 0.0D0
    WLS     = 0.0D0
    VT      = 0.0D0
    WT      = 0.0D0
    term_C  = 0.0D0
    
    vlo = 0.0D0; vnlo = 0.0D0; vnnlo = 0.0D0;
    cont_lo = 0.0D0; cont_nlo = 0.0D0;
    
    ! Contribution from the Delta isobar are described in EPJ A 32, 127-137 (2007)
    vc_delta    = 0.0D0
    wc_delta    = 0.0D0
    vt_delta    = 0.0D0
    wt_delta    = 0.0D0
    vs_delta    = 0.0D0
    ws_delta    = 0.0D0
    loop_H      = 0.0D0 
    loop_Sigma  = 0.0D0

    nx1 = all_orbit%nx(p)
    ny1 = all_orbit%ny(p)
    nz1 = all_orbit%nz(p)
    nx2 = all_orbit%nx(q)
    ny2 = all_orbit%ny(q)
    nz2 = all_orbit%nz(q)
    nx3 = all_orbit%nx(r)
    ny3 = all_orbit%ny(r)
    nz3 = all_orbit%nz(r)
    nx4 = all_orbit%nx(s)
    ny4 = all_orbit%ny(s)
    nz4 = all_orbit%nz(s)
    ! Conservation of linear momentum
    IF ( nx1 + nx2 /= nx3 + nx4 ) return 
    IF ( ny1 + ny2 /= ny3 + ny4 ) return 
    IF ( nz1 + nz2 /= nz3 + nz4 ) return 
    
    k1(1) = all_orbit%kx(p)
    k1(2) = all_orbit%ky(p)
    k1(3) = all_orbit%kz(p)
    k2(1) = all_orbit%kx(q)
    k2(2) = all_orbit%ky(q)
    k2(3) = all_orbit%kz(q)
    k3(1) = all_orbit%kx(r)
    k3(2) = all_orbit%ky(r)
    k3(3) = all_orbit%kz(r)
    k4(1) = all_orbit%kx(s)
    k4(2) = all_orbit%ky(s)
    k4(3) = all_orbit%kz(s)
    ! momenta in MeV
    k1 = k1*hbarc
    k2 = k2*hbarc
    k3 = k3*hbarc
    k4 = k4*hbarc
    
    m1 = all_orbit%sz(p) 
    m2 = all_orbit%sz(q) 
    m3 = all_orbit%sz(r) 
    m4 = all_orbit%sz(s) 
  
    t1 = all_orbit%tz(p) 
    t2 = all_orbit%tz(q) 
    t3 = all_orbit%tz(r) 
    t4 = all_orbit%tz(s) 
  
    ! RELATIVE MOMENTA <prel |v| pprel>
    prel  = 0.5d0*(k1-k2)
    pprel = 0.5d0*(k3-k4)
    p2  = sum(prel*prel) 
    pp2 = sum(pprel*pprel)
    ! AVERAGE MOMENTA kav = 1/2(prel+pprel)
    kmean = 0.5d0*( prel + pprel ) 
    kmean2 = sum(kmean*kmean)
    ! momentum transfer
    qtrans = prel - pprel
    q2 = sum(qtrans*qtrans)
    qabs = dsqrt(q2)
    
    ! [1], p. 28
    nucleon_mass = mnuc((t1+t2)/2)
    relativity_factor = nucleon_mass / sqrt( sqrt( (nucleon_mass**2 + p2) * (nucleon_mass**2 + pp2) ) )
    freg_nnlo = freg(prel, pprel, nexp_NN)
    IF ( freg_nnlo < 1.e-10 ) return

    ! cross product between momentum transfer and average momenta q X k
    qxk(1) = qtrans(2)*kmean(3)-qtrans(3)*kmean(2)
    qxk(2) = qtrans(3)*kmean(1)-qtrans(1)*kmean(3)
    qxk(3) = qtrans(1)*kmean(2)-qtrans(2)*kmean(1)
    ! define operator structures
    delta_spin    = delta_tab(m1,m3)*delta_tab(m2,m4)
    delta_isospin = delta_tab(t1,t3)*delta_tab(t2,t4)
    s_dot_s       = sigma_dot_sigma_tab(m1,m2,m3,m4)
    t_dot_t       = tau_dot_tau_tab(t1,t2,t3,t4)
    sigma1_dot_q_sigma2_dot_qtrans = chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,qtrans) ! used several times
    LS_struct     = chp_spin_dot_qxk_mtx(m1,m2,m3,m4,qxk)
    
    ! define loop functions
    IF ( regulatization_scheme == chp_SFR ) then ! SFR
       ! NLO loop functions
       loop_w = chp_NLO_two_pion_exchange_loop_w(q2, 2)
       loop_s = chp_NLO_sfr_two_pion_exchange_loop_s(2) ! [2], Eq. 2.16
       loop_L = chp_NLO_sfr_two_pion_exchange_loop_L(qabs, q2, loop_w, loop_s, 2)
       ! NNLO
       loop_wtilde2 = chp_NNLO_two_pion_exchange_loop_wtilde2(q2, 2)
       loop_A = chp_NNLO_sfr_two_pion_exchange_loop_A(qabs, q2, 2)
       ! Delta-full
       IF ( chiral_delta_flag ) then
          ! [5], Eq. 2.8
          loop_Sigma  = 2.d0*mpi2(2) + q2 - 2.d0*delta_split2
          loop_D = chp_delta_2PE_sfr_loop_D(q2, delta_split, 2)
          loop_H = chp_delta_loop_H(q2, delta_split, loop_w, loop_Sigma, loop_L, 2)
       end IF
    end IF

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !
    ! NEXT TO LEADING ORDER
    !
    vlo = 0.d0
    cont_lo = 0.d0

    WT = 0.d0 
    IF ( t1 + t2 == 0 ) THEN ! np, [1] Eq. 4.78 (T=1) and 4.79 (T=0)           
       WT = 0.d0 
       DO Tiso = 0, 2, 2
          cg1 = tjs(1,1,Tiso,t1,t2,-(t1+t2))*iph( (t1+t2)/2 )*sqrt(Tiso+1.d0)
          cg2 = tjs(1,1,Tiso,t3,t4,-(t3+t4))*iph( (t3+t4)/2 )*sqrt(Tiso+1.d0)     
          WT = WT + cg1*cg2*(-chp_one_pion_exchange(q2, 0) + iph(Tiso/2+1)*2.0D0*chp_one_pion_exchange(q2,1)) 
       end DO
    end IF
    IF ( t1 + t2 /= 0 ) THEN ! pp or nn, [1] Eq. 4.77
       WT = chp_one_pion_exchange(q2, 0) 
    end IF
    vlo = WT*sigma1_dot_q_sigma2_dot_qtrans

    ! leading order CIB contacts 
    term_CS = CS((t1+t2)/2) * delta_spin * delta_isospin
    term_CT = CT((t1+t2)/2) * s_dot_s    * delta_isospin
    cont_lo = term_CS + term_CT
    
    ! NLO contributions
    Vc      = 0.0D0
    Wc      = 0.0D0
    Vs      = 0.0D0
    Ws      = 0.0D0
    VT      = 0.0D0
    WT      = 0.0D0

    Vc_delta = 0.0D0
    Wc_delta = 0.0D0
    Vt_delta = 0.0D0
    Wt_delta = 0.0D0
    Vs_delta = 0.0D0
    Ws_delta = 0.0D0

    ! Delta-less terms (same in both power countings)
    Wc = chp_NLO_two_pion_exchange_Wc(q2, loop_L, loop_w, 2) ! [1], Eq. 4.9
    Vs = chp_NLO_two_pion_exchange_Vs(q2, loop_L, 2)         ! [1], Eq. 4.10
    VT = chp_NLO_two_pion_exchange_VT(loop_L, 2)             ! [1], Eq. 4.10
    
    ! Delta-full terms
    IF ( chiral_delta_flag ) then
       ! Vc 
       Vc_delta = &
            chp_NLO_delta_2PE_single_box_Vc_SFR(q2, loop_A, delta_split,2) +   &                     ! leading, single box, eq. 2.6
            chp_NLO_delta_2PE_double_box_Vc_SFR(loop_L, loop_Sigma, loop_D, loop_H, delta_split, 2)  ! leading, double box, eq. 2.7
       ! Wc (3 terms at NLO)
       Wc_delta = &
            chp_NLO_delta_2PE_single_triangle_Wc_SFR(loop_w, loop_L, loop_Sigma, loop_D, delta_split) + &      ! leading, triangle, eq. 2.5
            chp_NLO_delta_2PE_single_box_Wc_SFR(q2, loop_L, loop_Sigma, loop_D, delta_split, 2)  + &           ! leading, single box, eq. 2.6
            chp_NLO_delta_2PE_double_box_Wc_SFR(loop_w, loop_L, loop_Sigma, loop_D, loop_H, delta_split, 2)    ! leading, double box, eq. 2.7
       ! VT
       Vt_delta = &
            chp_NLO_delta_2PE_single_box_Vt_SFR(loop_w, loop_L, loop_D, delta_split, 2) + &   ! leading, single box, eq. 2.6
            chp_NLO_delta_2PE_double_box_Vt_SFR(loop_w, loop_L, loop_D, delta_split, 2)       ! leading, double box, eq. 2.7
       ! Vs
       Vs_delta = -q2 * Vt_delta
       ! WT
       Wt_delta = &
            chp_NLO_delta_2PE_single_box_Wt_SFR(loop_w, loop_A, delta_split, 2) + &      ! leading, single box, eq. 2.6
            chp_NLO_delta_2PE_double_box_Wt_SFR(loop_w, loop_L, loop_D, delta_split, 2)  ! leading, double box, eq. 2.7
       ! Ws
       Ws_delta = -q2 * Wt_delta
       !
       Vc = Vc + Vc_delta
       Wc = Wc + Wc_delta
       Vs = Vs + Vs_delta
       Ws = Ws + Ws_delta
       VT = VT + Vt_delta
       WT = WT + Wt_delta
    end IF
    
    ! multiply by operator structures
    Vc = Vc * delta_spin * delta_isospin
    Wc = Wc * delta_spin * t_dot_t
    Vs = Vs * s_dot_s    * delta_isospin
    Ws = Ws * s_dot_s    * t_dot_t
    VT = VT * sigma1_dot_q_sigma2_dot_qtrans * delta_isospin
    WT = WT * sigma1_dot_q_sigma2_dot_qtrans * t_dot_t
    
    IF ( chiral_delta_flag ) then
       vnlo = ( Vc + Wc + Vs + Ws + VT + WT )
    else
       vnlo = ( WC + VS + VT )
    end IF

    ! NLO contacts 
    term_C(1) = cnlo(1)*q2*delta_spin*delta_isospin 
    term_C(2) = cnlo(2)*kmean2*delta_spin*delta_isospin
    term_C(3) = cnlo(3)*q2*s_dot_s*delta_isospin
    term_C(4) = cnlo(4)*kmean2*s_dot_s*delta_isospin
    term_C(5) = cnlo(5)*LS_struct*delta_isospin
    term_C(6) = cnlo(6)*sigma1_dot_q_sigma2_dot_qtrans*delta_isospin
    term_C(7) = cnlo(7)*chp_sigma1_dot_q_sigma2_dot_q_mtx(m1,m2,m3,m4,kmean)*delta_isospin 
    cont_nlo = SUM(term_C)
    ! 
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    !
    ! NEXT TO NEXT TO LEADING ORDER 
    !
    Vc      = 0.0D0
    Wc      = 0.0D0
    Vs      = 0.0D0
    Ws      = 0.0D0
    VT      = 0.0D0
    WT      = 0.0D0
    VLS     = 0.0D0
    WLS     = 0.0D0
    
    Vc_delta = 0.0D0
    Wc_delta = 0.0D0
    Vt_delta = 0.0D0
    Wt_delta = 0.0D0
    Vs_delta = 0.0D0
    Ws_delta = 0.0D0
    
    ! The expressions of the NNLO functions differs according to the adopted power counting
    ! scheme. Machleidt 2011 [1] and Epelbaum [2,4]/Machleidt [3,6]  use different power countings.
    ! In particular, the way relativistic corrections proportional to 1/MN are included is
    ! different.
    ! The relevant expressions can be found e.g. in the following papers:
    ! Power counting Machleidt 2011 (e.g. N2LOopt,NNLOsat)
    ! [1] Eqs. 4.13-4.18 (+ corrections 4.21-4.24)    
    ! Power counting Epelbaum
    ! [6] Eqs. 47,48 (and 59-64 for the relativistic corrections)
    
    ! Machleidt power counting
    IF ( power_counting == pc_machleidt2011 ) then
       ! Delta-less contributions
       Vc  = chp_NNLO_two_pion_exchange_Vc (q2, loop_w, loop_A, loop_wtilde2, 2)
       Wc  = chp_NNLO_two_pion_exchange_Wc (q2, loop_w, loop_A, loop_wtilde2, 2)
       Vs  = chp_NNLO_two_pion_exchange_Vs (q2, loop_w, loop_A, loop_wtilde2, 2)
       Ws  = chp_NNLO_two_pion_exchange_Ws (q2, loop_w, loop_A, 2)
       VT  = chp_NNLO_two_pion_exchange_VT (loop_w, loop_A, loop_wtilde2, 2)
       WT  = chp_NNLO_two_pion_exchange_WT (q2, loop_w, loop_A, 2)
       VLS = chp_NNLO_two_pion_exchange_VLS(loop_A, loop_wtilde2, 2)
       WLS = chp_NNLO_two_pion_exchange_WLS(loop_w, loop_A, 2)
    ! Epelbaum power counting
    else IF ( power_counting==pc_epelbaum2005 ) then
       ! Delta-less contributions -> loop_d 2PE diagrams
       Vc = chp_two_pion_exchange_1loop_d_Vc (q2, loop_A, loop_wtilde2, 2)
       WT = chp_two_pion_exchange_1loop_d_WT (    loop_w, loop_A )
       Ws = chp_two_pion_exchange_1loop_d_Ws (q2, loop_w, loop_A )
       ! Delta-full contributions
       IF ( chiral_delta_flag ) then
          Vc_delta = chp_NNLO_delta_2PE_triangle_Vc_SFR(loop_w, loop_L, loop_Sigma, loop_D, delta_split, 2) ! subleading, triangle, [5], eq.2.9
          Wt_delta = chp_NNLO_delta_2PE_triangle_Wt_SFR(loop_w, loop_L, loop_D, delta_split) ! subleading, triangle, [5], eq.2.9
          Ws_delta = -q2 * Wt_delta
          Vc = Vc + Vc_delta
          WT = WT + Wt_delta
          Ws = Ws + Ws_delta
       end IF
    end IF
    
    ! multiply by operator structures
    Vc  = Vc * delta_spin * delta_isospin
    Wc  = Wc * delta_spin * t_dot_t
    Vs  = Vs * s_dot_s    * delta_isospin
    Ws  = Ws * s_dot_s    * t_dot_t
    VT  = VT * sigma1_dot_q_sigma2_dot_qtrans * delta_isospin
    WT  = WT * sigma1_dot_q_sigma2_dot_qtrans * t_dot_t
    VLS = VLS* LS_struct  * delta_isospin
    WLS = WLS* LS_struct  * t_dot_t

    IF ( power_counting == pc_machleidt2011 ) then 
       vnnlo = (Vc + Wc + VLS + WLS + VT + WT + Vs + Ws)
    elseIF ( power_counting == pc_epelbaum2005 ) then
       vnnlo = (Vc + WT + Ws)
    end IF
    ! 
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    ! sum up all orders
    ! regulator and relativity factor
    vdir = (vlo + cont_lo + vnlo + cont_nlo + vnnlo ) * & 
         relativity_factor * freg_nnlo
    res = hbarc**3 * (vdir)/volume
    
  end FUNCTION chiral_NN
  
  
  FUNCTION chiral_3N(p,q,r,s,t,u) result(res)
    USE single_particle_orbits
    USE constants
    USE chiral_constants
    USE chiral_tables 
    USE ang_mom_functions, only : tjs 
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: p,q,r,s,t,u
    COMPLEX*16 :: res
    INTEGER :: nx1,ny1,nz1, nx2,ny2,nz2, nx3,ny3,nz3
    INTEGER :: nx4,ny4,nz4, nx5,ny5,nz5, nx6,ny6,nz6
    INTEGER :: n1(3),n2(3),n3(3),n4(3),n5(3),n6(3)
    INTEGER :: ntrans1(3), ntrans2(3), ntrans3(3)
    INTEGER :: m1,m2,m3,m4,m5,m6, t1,t2,t3,t4,t5,t6
    REAL*8 :: k1(3),k2(3),k3(3),k4(3),k5(3),k6(3)
    REAL*8 :: q1xq2(3), q2xq3(3), q1xq3(3), q2sq, q3sq, q1sq
    REAL*8 :: qtrans1(3), qtrans2(3), qtrans3(3)
    REAL*8 :: freg1, freg2, freg3, freg_nl1, freg_nl2
    COMPLEX*16 :: opstruct_1, opstruct_2, opstruct_3
    COMPLEX*16 :: sigma_dot_q141, sigma_dot_q142, sigma_dot_q143
    COMPLEX*16 :: sigma_dot_q251, sigma_dot_q252, sigma_dot_q253
    COMPLEX*16 :: sigma_dot_q361, sigma_dot_q362, sigma_dot_q363
    COMPLEX*16 :: sigma_dot_q_ope141, sigma_dot_q_ope252, sigma_dot_q_ope363
    COMPLEX*16 :: tauXtau_dot_tau_1, tauXtau_dot_tau_2, tauXtau_dot_tau_3
    COMPLEX*16 :: V3NF_sum, V3NF_k3, V3NF_k2, V3NF_k1
    
    res = 0.d0

    ! Conservation of isospin 
    IF ( all_orbit%tz(p) + all_orbit%tz(q)+ all_orbit%tz(r) /= & 
         all_orbit%tz(s) + all_orbit%tz(t)+ all_orbit%tz(u) ) return
    
    tauXtau_dot_tau_1 = 0.d0
    tauXtau_dot_tau_2 = 0.d0
    tauXtau_dot_tau_3 = 0.d0
    
    nx1 = all_orbit%nx(p)
    ny1 = all_orbit%ny(p)
    nz1 = all_orbit%nz(p)
    nx2 = all_orbit%nx(q)
    ny2 = all_orbit%ny(q)
    nz2 = all_orbit%nz(q)
    nx3 = all_orbit%nx(r)
    ny3 = all_orbit%ny(r)
    nz3 = all_orbit%nz(r)
    nx4 = all_orbit%nx(s)
    ny4 = all_orbit%ny(s)
    nz4 = all_orbit%nz(s)
    nx5 = all_orbit%nx(t)
    ny5 = all_orbit%ny(t)
    nz5 = all_orbit%nz(t)
    nx6 = all_orbit%nx(u)
    ny6 = all_orbit%ny(u)
    nz6 = all_orbit%nz(u)
    ! Conservation of linear momentum
    IF ( nx1 + nx2 + nx3 /= nx4 + nx5 + nx6 ) return 
    IF ( ny1 + ny2 + ny3 /= ny4 + ny5 + ny6 ) return 
    IF ( nz1 + nz2 + nz3 /= nz4 + nz5 + nz6 ) return 
    
    k1(1) = all_orbit%kx(p)
    k1(2) = all_orbit%ky(p)
    k1(3) = all_orbit%kz(p)
    k2(1) = all_orbit%kx(q)
    k2(2) = all_orbit%ky(q)
    k2(3) = all_orbit%kz(q)
    k3(1) = all_orbit%kx(r)
    k3(2) = all_orbit%ky(r)
    k3(3) = all_orbit%kz(r)
    k4(1) = all_orbit%kx(s)
    k4(2) = all_orbit%ky(s)
    k4(3) = all_orbit%kz(s)
    k5(1) = all_orbit%kx(t)
    k5(2) = all_orbit%ky(t)
    k5(3) = all_orbit%kz(t)
    k6(1) = all_orbit%kx(u)
    k6(2) = all_orbit%ky(u)
    k6(3) = all_orbit%kz(u)    
    ! momenta in MeV
    k1 = k1*hbarc
    k2 = k2*hbarc
    k3 = k3*hbarc
    k4 = k4*hbarc
    k5 = k5*hbarc
    k6 = k6*hbarc
    
    m1 = all_orbit%sz(p) 
    m2 = all_orbit%sz(q) 
    m3 = all_orbit%sz(r) 
    m4 = all_orbit%sz(s) 
    m5 = all_orbit%sz(t) 
    m6 = all_orbit%sz(u) 
  
    t1 = all_orbit%tz(p) 
    t2 = all_orbit%tz(q) 
    t3 = all_orbit%tz(r) 
    t4 = all_orbit%tz(s) 
    t5 = all_orbit%tz(t) 
    t6 = all_orbit%tz(u) 

    n1(1) = nx1
    n1(2) = ny1
    n1(3) = nz1
    n2(1) = nx2
    n2(2) = ny2
    n2(3) = nz2
    n3(1) = nx3
    n3(2) = ny3
    n3(3) = nz3
    n4(1) = nx4
    n4(2) = ny4
    n4(3) = nz4
    n5(1) = nx5
    n5(2) = ny5
    n5(3) = nz5
    n6(1) = nx6
    n6(2) = ny6
    n6(3) = nz6
    
    ! Momentum transfers
    qtrans1 = k4-k1
    qtrans2 = k5-k2
    qtrans3 = k6-k3    
    ntrans1 = n4 - n1
    ntrans2 = n5 - n2
    ntrans3 = n6 - n3
    q1sq = sum(qtrans1*qtrans1)
    q2sq = sum(qtrans2*qtrans2)
    q3sq = sum(qtrans3*qtrans3)
    
    ! local regulators
    freg1 = freg_3NFlocal(q1sq)
    freg2 = freg_3NFlocal(q2sq) 
    freg3 = freg_3NFlocal(q3sq)
    freg_nl1 = freg_3NFnonlocal(k1,k2,k3)
    freg_nl2 = freg_3NFnonlocal(k4,k5,k6)
    IF ( max(freg1*freg2, freg1*freg3, freg2*freg3) < 1.e-12 ) return
    IF ( abs(freg_3NFnonlocal(k1,k2,k3)*freg_3NFnonlocal(k4,k5,k6)) < 1.e-12 ) return
    
    !  cross product between momentum transfer and average momenta q X k 
    q1xq2(1) = qtrans1(2)*qtrans2(3)-qtrans1(3)*qtrans2(2) 
    q1xq2(2) = qtrans1(3)*qtrans2(1)-qtrans1(1)*qtrans2(3) 
    q1xq2(3) = qtrans1(1)*qtrans2(2)-qtrans1(2)*qtrans2(1) 
    
    q1xq3(1) = qtrans1(2)*qtrans3(3)-qtrans1(3)*qtrans3(2) 
    q1xq3(2) = qtrans1(3)*qtrans3(1)-qtrans1(1)*qtrans3(3) 
    q1xq3(3) = qtrans1(1)*qtrans3(2)-qtrans1(2)*qtrans3(1) 
    
    q2xq3(1) = qtrans2(2)*qtrans3(3)-qtrans2(3)*qtrans3(2) 
    q2xq3(2) = qtrans2(3)*qtrans3(1)-qtrans2(1)*qtrans3(3) 
    q2xq3(3) = qtrans2(1)*qtrans3(2)-qtrans2(2)*qtrans3(1) 
    
    ! Ex. 251 = m2, m5, ntrans1 (qtrans1);    143 = m1, m4, ntrans3 (qtrans3)
    sigma_dot_q251 = sigma_dot_q_tab(m2,m5,ntrans1(1),ntrans1(2),ntrans1(3))
    sigma_dot_q141 = sigma_dot_q_tab(m1,m4,ntrans1(1),ntrans1(2),ntrans1(3))
    sigma_dot_q361 = sigma_dot_q_tab(m3,m6,ntrans1(1),ntrans1(2),ntrans1(3))
    sigma_dot_q252 = sigma_dot_q_tab(m2,m5,ntrans2(1),ntrans2(2),ntrans2(3))
    sigma_dot_q142 = sigma_dot_q_tab(m1,m4,ntrans2(1),ntrans2(2),ntrans2(3))
    sigma_dot_q362 = sigma_dot_q_tab(m3,m6,ntrans2(1),ntrans2(2),ntrans2(3))  
    sigma_dot_q363 = sigma_dot_q_tab(m3,m6,ntrans3(1),ntrans3(2),ntrans3(3))
    sigma_dot_q143 = sigma_dot_q_tab(m1,m4,ntrans3(1),ntrans3(2),ntrans3(3))
    sigma_dot_q253 = sigma_dot_q_tab(m2,m5,ntrans3(1),ntrans3(2),ntrans3(3))

    sigma_dot_q_ope141 = sigma_dot_q_ope_tab(m1,m4,ntrans1(1),ntrans1(2),ntrans1(3))
    sigma_dot_q_ope252 = sigma_dot_q_ope_tab(m2,m5,ntrans2(1),ntrans2(2),ntrans2(3))
    sigma_dot_q_ope363 = sigma_dot_q_ope_tab(m3,m6,ntrans3(1),ntrans3(2),ntrans3(3))

    ! tau dot tau terms
    opstruct_1 = tau_dot_tau_tab(t2,t3,t5,t6)*delta_tab(m1,m4)*delta_tab(t1,t4)
    opstruct_2 = tau_dot_tau_tab(t1,t3,t4,t6)*delta_tab(m2,m5)*delta_tab(t2,t5)
    opstruct_3 = tau_dot_tau_tab(t1,t2,t4,t5)*delta_tab(m3,m6)*delta_tab(t3,t6)
    ! tauXtau dot tau terms
    tauXtau_dot_tau_1 = tau1_dot_tauXtau_tab(t1,t2,t3,t4,t5,t6)
    tauXtau_dot_tau_2 = tau1_dot_tauXtau_tab(t2,t1,t3,t5,t4,t6)
    tauXtau_dot_tau_3 = tau1_dot_tauXtau_tab(t3,t1,t2,t6,t4,t5)
    
    ! permutations (ijk) = (123) + (213)
    ! factor 1/2 in cD are omitted because only three permutations are summed (123 and 213 are equal contributions etc.)
    ! const(23) = gA/(8*fpi^2)

    V3NF_k3 = 0.d0
    IF ( abs(opstruct_3) .ge. 1e-12 ) then
       V3NF_k3 = &
            ! factor out [tau.tau * delta(m3,m6) ]
            opstruct_3 * ( &
            ! vce terms
            Econst * delta_tab(m1,m4)*delta_tab(m2,m5) + &
            ! vcd terms
            - Dconst * const(23) * ( &
            ! 123
            sigma_dot_q252 * sigma_dot_q142/( q2sq + mpi2(2) ) + &
            ! 213
            sigma_dot_q141 * sigma_dot_q251/( q1sq + mpi2(2) ) ) + &
            ! V2pi terms (c1,c3), 123
            sigma_dot_q_ope141 * sigma_dot_q_ope252 * const(2) * ( -c1_3NF * (4.d0*mpi2(2)/fpi2) + &
            c3_3NF * (2.d0/fpi2) * sum(qtrans1*qtrans2) ) &
            )
    end IF
    ! V2pi terms (c4), 123
    IF ( abs(tauXtau_dot_tau_3) .ge. 1e-12 ) then
       V3NF_k3 = V3NF_k3 + &
            c4_3NF * (1.d0/fpi2) * const(2) * sigma_dot_q_ope141 * sigma_dot_q_ope252 * &
            chp_sigma_dot_q_mtx(m3,m6,q1xq2) * tauXtau_dot_tau_3
    end IF
    
    ! permutations (312) + (132)
    V3NF_k2 = 0.d0
    IF ( abs(opstruct_2) .ge. 1e-12 ) then
       V3NF_k2 = &
            ! factor out [tau.tau * delta_tab(r,u)] from vce and vcd
            opstruct_2 * ( &
            ! vce terms
            Econst * delta_tab(m1,m4)*delta_tab(m3,m6) + &
            ! vcd terms, 312
            - Dconst * const(23) * ( sigma_dot_q363 * sigma_dot_q143/( q3sq + mpi2(2) ) + &
            ! 132
            sigma_dot_q141 * sigma_dot_q361/( q1sq + mpi2(2) ) ) + &
            ! V2pi terms (c1,c3), 312
            sigma_dot_q_ope141 * sigma_dot_q_ope363 * const(2) * ( -c1_3NF * (4.d0*mpi2(2)/fpi2) + &
            c3_3NF * (2.d0/fpi2) * sum(qtrans1*qtrans3) ) &
            )
    end IF
    ! V2pi terms (c4), 312
    IF ( abs(tauXtau_dot_tau_2) .ge. 1e-12 ) then
       V3NF_k2 = V3NF_k2 +     &
            c4_3NF * (1.d0/fpi2) * const(2) * sigma_dot_q_ope141 * sigma_dot_q_ope363 * &
            chp_sigma_dot_q_mtx(m2,m5,q1xq3) * tauXtau_dot_tau_2
    end IF
    
    ! permutations (321) + (231)
    V3NF_k1 = 0.d0
    IF ( abs(opstruct_1) .ge. 1e-12 ) then
       V3NF_k1 = &
            ! factor out [tau.tau * delta_tab(r,u)] from vce and vcd
            opstruct_1 * ( &
            ! vce terms
            Econst * delta_tab(m2,m5)*delta_tab(m3,m6) + &
            ! vcd terms, 321
            - Dconst * const(23) * ( sigma_dot_q363 * sigma_dot_q253/( q3sq + mpi2(2) ) + &
            ! 231
            sigma_dot_q252 * sigma_dot_q362/( q2sq + mpi2(2) ) ) + &
            ! V2pi terms (c1,c3), 321
            sigma_dot_q_ope252 * sigma_dot_q_ope363 * const(2) * ( -c1_3NF * (4.d0*mpi2(2)/fpi2) + &
            c3_3NF * (2.d0/fpi2) * sum(qtrans2*qtrans3) ) &
            )
    end IF
    ! V2pi terms (c4), 231
    IF ( abs(tauXtau_dot_tau_1) .ge. 1e-12 ) then
       V3NF_k1 = V3NF_k1 + &
            c4_3NF * (1.d0/fpi2) * const(2) * sigma_dot_q_ope252 * sigma_dot_q_ope363 * &
            chp_sigma_dot_q_mtx(m1,m4,q2xq3) * tauXtau_dot_tau_1
    end IF

    V3NF_sum = &
         freg1*freg2 * V3NF_k3 + & 
         freg1*freg3 * V3NF_k2 + & 
         freg2*freg3 * V3NF_k1
    
    res = freg_nl1 * freg_nl2 * V3NF_sum
    res = res * hbarc**6/volume**2
    
  end FUNCTION chiral_3N

  
  FUNCTION chiral_3N_asym(i,j,p,k,l,q) result(res)
    USE single_particle_orbits
    USE constants
    USE chiral_constants
    IMPLICIT none 
    INTEGER, intent(in) :: i,j,p,k,l,q
    COMPLEX*16 :: res
    
    res = ( chiral_3N(i,j,p,k,l,q) - chiral_3N(i,j,p,l,k,q) &
         - chiral_3N(i,j,p,q,l,k) - chiral_3N(i,j,p,k,q,l) &
         + chiral_3N(i,j,p,l,q,k) + chiral_3N(i,j,p,q,k,l) )
    
  end FUNCTION chiral_3N_asym
  
  
  FUNCTION vmom_minnesota(p,q,r,s) result(res)
    USE single_particle_orbits
    USE constants
    USE chiral_constants
    
    IMPLICIT none
    REAL*8 :: res
    INTEGER :: p,q,r,s, m1,m2,m3,m4, t1,t2,t3,t4
    INTEGER :: nx1, ny1, nz1, nx2, ny2, nz2, nx3, ny3, nz3, nx4, ny4, nz4 
    REAL*8 :: k1(3), k2(3), k3(3), k4(3)
    REAL*8 :: qtrans(3), prel(3), pprel(3)
    REAL*8 :: q2, vdir, delta
    REAL*8 :: v0r, v0t, v0s, kr,kt,ks, vr,vt,vs, vcentral, vsigma

    res = 0.d0
    
    nx1 = all_orbit%nx(p)
    ny1 = all_orbit%ny(p)
    nz1 = all_orbit%nz(p)
    nx2 = all_orbit%nx(q)
    ny2 = all_orbit%ny(q)
    nz2 = all_orbit%nz(q)
    nx3 = all_orbit%nx(r)
    ny3 = all_orbit%ny(r)
    nz3 = all_orbit%nz(r)
    nx4 = all_orbit%nx(s)
    ny4 = all_orbit%ny(s)
    nz4 = all_orbit%nz(s)
    
    ! Conservation of linear momentum
    IF ( nx1 + nx2 /= nx3 + nx4 ) return 
    IF ( ny1 + ny2 /= ny3 + ny4 ) return 
    IF ( nz1 + nz2 /= nz3 + nz4 ) return 
    ! conservation of spin and isospin 
    IF ( all_orbit%tz(p) + all_orbit%tz(q) /= all_orbit%tz(r) + all_orbit%tz(s) ) return 
  
    k1(1) = all_orbit%kx(p)
    k1(2) = all_orbit%ky(p)
    k1(3) = all_orbit%kz(p)
    k2(1) = all_orbit%kx(q)
    k2(2) = all_orbit%ky(q)
    k2(3) = all_orbit%kz(q)
    k3(1) = all_orbit%kx(r)
    k3(2) = all_orbit%ky(r)
    k3(3) = all_orbit%kz(r)
    k4(1) = all_orbit%kx(s)
    k4(2) = all_orbit%ky(s)
    k4(3) = all_orbit%kz(s)
  
    m1 = all_orbit%sz(p) 
    m2 = all_orbit%sz(q) 
    m3 = all_orbit%sz(r) 
    m4 = all_orbit%sz(s) 
  
    t1 = all_orbit%tz(p) 
    t2 = all_orbit%tz(q) 
    t3 = all_orbit%tz(r) 
    t4 = all_orbit%tz(s) 
  
    ! RELATIVE MOMENTA <prel |v| pprel > 
    prel  = 0.5d0*(k1-k2)
    pprel = 0.5d0*(k3-k4)
    ! momentum transfer 
    qtrans = prel - pprel
    q2 = sum(qtrans*qtrans)     
    
    v0r=200.0  ! MeV
    v0t=178.0  ! MeV
    v0s=91.85  ! MeV
    kr=1.487  ! fm**-2
    kt=0.639  ! fm**-2
    ks=0.465  ! fm**-2
    
    vr =  v0r * pi**1.5d0 * exp(-q2/(4.d0*kr) )/ (kr**1.5d0) 
    vt = -v0t * pi**1.5d0 * exp(-q2/(4.d0*kt) )/ (kt**1.5d0)
    vs = -v0s * pi**1.5d0 * exp(-q2/(4.d0*ks) )/ (ks**1.5d0)
    
    vcentral = 0.25d0*(vr+vs)*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) 
    vsigma = -0.25d0*(vr+vs)*real(chp_sigma_dot_sigma(m1,m2,m3,m4))*delta(t1,t3)*delta(t2,t4)
    
    vr = vr * (3.d0*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) - & 
         real(chp_tau_dot_tau(t1,t2,t3,t4))*delta(m1,m3)*delta(m2,m4) - & 
         real(chp_sigma_dot_sigma(m1,m2,m3,m4))*delta(t1,t3)*delta(t2,t4) - & 
         real(chp_tau_dot_tau(t1,t2,t3,t4))*real(chp_sigma_dot_sigma(m1,m2,m3,m4)) )/8.d0
    
    vs = vs * (3.d0*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) + & 
         real(chp_tau_dot_tau(t1,t2,t3,t4))*delta(m1,m3)*delta(m2,m4) - & 
         3.d0*real(chp_sigma_dot_sigma(m1,m2,m3,m4))*delta(t1,t3)*delta(t2,t4) - & 
         real(chp_tau_dot_tau(t1,t2,t3,t4))*real(chp_sigma_dot_sigma(m1,m2,m3,m4)) )/16.d0
    
    vt = vt * (3.d0*delta(m1,m3)*delta(m2,m4)*delta(t1,t3)*delta(t2,t4) - & 
         3.d0*real(chp_tau_dot_tau(t1,t2,t3,t4))*delta(m1,m3)*delta(m2,m4) + & 
         real(chp_sigma_dot_sigma(m1,m2,m3,m4))*delta(t1,t3)*delta(t2,t4) - & 
         real(chp_tau_dot_tau(t1,t2,t3,t4))*real(chp_sigma_dot_sigma(m1,m2,m3,m4)) )/16.d0
    
    vdir = vs+vr+vt
    res = vdir/volume
    
  end FUNCTION vmom_minnesota

  
  FUNCTION v2int(p,q,r,s) result(v2b)
    USE chiral_constants
    USE single_particle_orbits
    ! USE NNForceCartesian
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: p,q,r,s
    INTEGER :: pp(5), qq(5), rr(5), ss(5)
    COMPLEX*16 :: v2b
    v2b = 0.d0
    
    IF ( which_pot > 0 ) then
       IF ( which_pot < 7 ) then
          v2b = chiral_NN(p,q,r,s) - chiral_NN(p,q,s,r)
       ! else
       !    pp = [ all_orbit%nx(p), all_orbit%ny(p), all_orbit%nz(p), all_orbit%sz(p), all_orbit%tz(p) ]
       !    qq = [ all_orbit%nx(q), all_orbit%ny(q), all_orbit%nz(q), all_orbit%sz(q), all_orbit%tz(q) ]
       !    rr = [ all_orbit%nx(r), all_orbit%ny(r), all_orbit%nz(r), all_orbit%sz(r), all_orbit%tz(r) ]
       !    ss = [ all_orbit%nx(s), all_orbit%ny(s), all_orbit%nz(s), all_orbit%sz(s), all_orbit%tz(s) ]
       !    v2b = vrel%GetME(pp,qq,rr,ss)
       end IF
    else
       v2b = vmom_minnesota(p,q,r,s) - vmom_minnesota(p,q,s,r)
    end IF
    
  end FUNCTION v2int
  
  FUNCTION v3int(p,q,r,s,t,u) result(v3b)
    USE chiral_constants
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: p,q,r,s,t,u
    COMPLEX*16 :: v3b
    
    IF ( which_pot > 0 ) then
       v3b = chiral_3N_asym(p,q,r,s,t,u)
    else
       v3b = 0.d0
    end IF
    
  end FUNCTION v3int


end MODULE chiral_potentials


  ! subroutine set_nonlocal_int(vint, v)
  !   type(NNForceMom), intent(inout) :: vint
  !   type(v_inter), intent(in) :: v
  !   type(TwoBodyRelSpaceSpinMBasis), pointer :: two
  !   type(TwoBodyRelChanSpinMBasis), pointer :: ch_s0, ch_s1, ch_cp
  !   integer :: ich_s0, ich_s1, ich_cp
  !   integer :: z, j, bra, ket, n1, n2
  !   real(8) :: v6(6)

  !   two => vint%ms
  !   do z = -1, 1
  !     do j = 0, two%GetJmax()

  !       !$omp parallel
  !       !$omp do private(ket, bra, v6, ich_s0, ich_s1, ich_cp, &
  !       !$omp &  ch_s0, ch_s1, ch_cp, n1, n2)
  !       do ket = 1, two%GetNMesh()
  !         do bra = 1, two%GetNMesh()

  !           v6(:) = 0.d0
  !           v6 = v%v(bra,ket,j,z,:)
  !           ich_s0 = two%GetIndex(j, (-1)**j, 0, z)
  !           ich_s1 = two%GetIndex(j, (-1)**j, 1, z)
  !           ich_cp = two%GetIndex(j, (-1)**(j+1), 1, z)
  !           ch_s0 => two%GetChannel(j, (-1)**j, 0, z)
  !           ch_s1 => two%GetChannel(j, (-1)**j, 1, z)
  !           ch_cp => two%GetChannel(j, (-1)**(j+1), 1, z)
  !           if(associated( ch_s0 )) then
  !             n1 = ch_s0%GetIndex(bra, j)
  !             n2 = ch_s0%GetIndex(ket, j)
  !             vint%MatCh(ich_s0)%m(n1,n2) = v6(1)
  !           end if

  !           if(associated( ch_s1 )) then
  !             n1 = ch_s1%GetIndex(bra, j)
  !             n2 = ch_s1%GetIndex(ket, j)
  !             vint%MatCh(ich_s1)%m(n1,n2) = v6(2)
  !           end if

  !           if(associated( ch_cp )) then
  !             n1 = ch_cp%GetIndex(bra, j+1)
  !             n2 = ch_cp%GetIndex(ket, j+1)
  !             vint%MatCh(ich_cp)%m(n1,n2) = v6(3)
  !             if(j-1 < 0) cycle
  !             n1 = ch_cp%GetIndex(bra, j-1)
  !             n2 = ch_cp%GetIndex(ket, j-1)
  !             vint%MatCh(ich_cp)%m(n1,n2) = v6(4)
  !             n1 = ch_cp%GetIndex(bra, j+1)
  !             n2 = ch_cp%GetIndex(ket, j-1)
  !             vint%MatCh(ich_cp)%m(n1,n2) = v6(5)
  !             n1 = ch_cp%GetIndex(bra, j-1)
  !             n2 = ch_cp%GetIndex(ket, j+1)
  !             vint%MatCh(ich_cp)%m(n1,n2) = v6(6)

  !           end if
  !         end do
  !       end do
  !       !$omp end do
  !       !$omp end parallel

  !     end do
  !   end do
  ! end subroutine set_nonlocal_int


  ! subroutine vnn_n2lo_emn500(vint)
  !   use MyLibrary, only: hc
  !   type(NNForceMom), intent(inout) :: vint
  !   type(TwoBodyRelSpaceSpinMBasis), pointer :: two
  !   integer :: nmesh, jmax
  !   real(8)  :: v, xmev, ymev
  !   real(8) :: vpp(6), vnp(6), vnn(6)
  !   integer :: j, inn, ij, bra, ket, i
  !   real(8), allocatable :: pmom(:)
  !   type(v_inter) :: vmom
  !   real(8) :: hbc3
  !   character (len=4) :: label
  !   common /cnn/ inn
  !   common /cpot/ v(6),xmev,ymev
  !   common /cstate/ j, heform, sing, trip, coup, endep, label
  !   logical :: sing, trip, coup, heform, endep
  !   two => vint%ms
  !   nmesh = two%NMesh
  !   jmax = two%jmax

  !   heform=.false.
  !   sing=.true.
  !   trip=.true.
  !   coup=.true.
  !   hbc3=hc**3
  !   j=0
  !   xmev=0.d0
  !   ymev=0.d0
  !   v(:)=0.d0
  !   allocate(vmom%v(nmesh, nmesh, 0:jmax, -1:1, 6))
  !   allocate(pmom(NMesh))
  !   do i = 1, nmesh
  !     pmom(i) = two%jpsz(1)%points(i)%GetP()
  !   end do
  !   vmom%v = 0.d0

  !   do ij = 0, jmax
  !     do bra = 1, nmesh
  !       do ket = 1, nmesh
  !         xmev = pmom(bra) * hc
  !         ymev = pmom(ket) * hc
  !         j = ij

  !         inn = 1
  !         call n2lo500
  !         vpp = v * hbc3

  !         inn = 2
  !         call n2lo500
  !         vnp = v * hbc3

  !         inn = 3
  !         call n2lo500
  !         vnn = v * hbc3
  !         vmom%v(bra, ket, ij, -1, :) = vpp
  !         vmom%v(bra, ket, ij,  0, :) = vnp
  !         vmom%v(bra, ket, ij,  1, :) = vnn
  !       end do
  !     end do
  !   end do
  !   call set_nonlocal_int(vint, vmom)
  !   deallocate(vmom%v,pmom)
  ! end subroutine vnn_n2lo_emn500
