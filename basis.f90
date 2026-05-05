

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE allocate_sp_data
  USE parallel
  USE single_particle_orbits
  USE constants
  USE one_body_operators
  USE chiral_potentials
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: shell, nx,ny,nz,ix,iy,iz, count
  INTEGER, DIMENSION (39) :: Nmaxlist
  Nmaxlist = (/ 0,1,2,3,4,5,6,8,9,10,11,12,13,14,16,&
                17,18,19,20,21,22,24,25,26,27,29,30,&
                32,33,34,35,36,37,38,40,41,42,43,44 /)

  Nmax0 = Nmaxlist(Nmax)
  
  szmax =  1
  szmin = -1
  tzmin =  1
  tzmax = -1
  IF ( below_ef_n > 0 ) tzmax =  1
  IF ( below_ef_p > 0 ) tzmin = -1
  
  lx = (real(below_ef)/rho)**(1.d0/3.d0)
  ly = lx
  lz = lx

  rho_n = real(below_ef_n)/(lx*ly*lz)
  rho_p = real(below_ef_p)/(lx*ly*lz)
  kf_n = (3.d0*pi**2*rho_n)**(1.d0/3.d0)
  kf_p = (3.d0*pi**2*rho_p)**(1.d0/3.d0)
  kf   = (3.d0*pi**2*rho)**(1.d0/3.d0)

  IF ( iam == 0 ) then
     write(6,'(A12,9x,F7.4,3x,A2)') 'Box Length:', lx, 'fm'
     IF ( below_ef_n > 0 ) write(6,'(A11,10x,F7.4,3x,A5)') 'N Density:', rho_n, 'fm^-3'
     IF ( below_ef_p > 0 ) write(6,'(A11,10x,F7.4,3x,A5)') 'P Density:', rho_p, 'fm^-3'
     write(6,'(A9,12x,F7.4,3x,A5)') 'k_Fermi:', kf, 'fm^-1'
     IF ( below_ef_n > 0 ) write(6,'(A11,10x,F7.4,3x,A5)') 'N k_Fermi:', kf_n, 'fm^-1'
     IF ( below_ef_p > 0 ) write(6,'(A11,10x,F7.4,3x,A5)') 'P k_Fermi:', kf_p, 'fm^-1'
     write(6,*)
  end IF
  volume = lx**3
  
  ! count states so that nx^2 + ny^2 + nz^2 <= Nmax0
  count = 0
  n1max = 0
  DO shell = 0, Nmax0
     DO nx = 0, Nmax0
        DO ny = 0, Nmax0
           DO nz = 0, Nmax0
              IF ( nx**2 + ny**2 + nz**2 /= shell ) cycle
              IF ( nx > n1max ) n1max = nx
              DO ix = -min(nx,1), min(nx,1), 2
                 DO iy = -min(ny,1), min(ny,1), 2
                    DO iz = -min(nz,1), min(nz,1), 2
                       ! IF ( iam == 0 ) write(6,*) 'state:', ix*nx, iy*ny, iz*nz, nx**2 + ny**2 + nz**2
                       count = count + 1
                    end DO
                 end DO
              end DO
           end DO
        end DO
     end DO
  end DO

  all_orbit%total_orbits = (szmax-szmin+2) * (tzmax-tzmin+2) * count/4
  CALL allocate_sp_array(all_orbit,all_orbit%total_orbits)
  CALL mem_register('base', REAL(all_orbit%total_orbits * (5 * 4.d0 + 4 * 8.d0 + 10 * 1.d0), dp))
  
  tot_orbs = all_orbit%total_orbits
  IF ( iam == 0 ) write(6,*) 'Total number of states', tot_orbs

  tot_orbs_n = tot_orbs
  tot_orbs_p = 0
  IF ( below_ef_p > 0 ) then
     tot_orbs_n = tot_orbs/2
     tot_orbs_p = tot_orbs/2
  end IF

end SUBROUTINE allocate_sp_data


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE setup_sp_data
  USE parallel
  USE single_particle_orbits
  USE constants
  USE one_body_operators
  USE chiral_potentials
  USE build_status
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: i, kk, shell, maxn, N2
  INTEGER :: sum_sz,sum_tz, sz,tzp
  INTEGER :: nx,ny,nz,ix,iy,iz
  INTEGER :: count_n, count_p, ih, ip
  REAL(dp) :: k2
  INTEGER, ALLOCATABLE :: check_orbital(:), order(:)
  TYPE (single_particle_descript) :: tmp_orbit
  INTEGER, DIMENSION (45) :: Nmaxlist0
  Nmaxlist0 = (/ 0,1,2,3,4,5,6,6,7,8,9,10,11,12,13,13,14,&
                15,16,17,18,19,20,20,21,22,23,24,24,25,26,&
                26,27,28,29,30,30,32,33,34,35,36,37,38,39 /)

  kk = 0
  DO shell = 0, Nmax0
     DO maxn = 0, Nmax0
        DO nx = 0, Nmax0
           DO ny = 0, Nmax0
              DO nz = 0, Nmax0
                 IF ( max(nx,ny,nz) /= maxn ) cycle
                 IF ( nx**2 + ny**2 + nz**2 /= shell ) cycle
                 DO ix = -min(nx,1), min(nx,1), 2
                    DO iy = -min(ny,1), min(ny,1), 2
                       DO iz = -min(nz,1), min(nz,1), 2
                          
                          DO sz = szmin, szmax, 2
                             DO tzp = tzmin, tzmax, 2
                                
                                kk = kk + 1
                                all_orbit%kx(kk) = (2.d0*ix*nx*pi)/lx
                                all_orbit%ky(kk) = (2.d0*iy*ny*pi)/ly
                                all_orbit%kz(kk) = (2.d0*iz*nz*pi)/lz
                                
                                all_orbit%nx(kk) = ix*nx
                                all_orbit%ny(kk) = iy*ny
                                all_orbit%nz(kk) = iz*nz
                                
                                all_orbit%sz(kk) = sz
                                all_orbit%tz(kk) = tzp
                                
                                k2 = all_orbit%kx(kk)**2 + all_orbit%ky(kk)**2 + all_orbit%kz(kk)**2
                                all_orbit%e(kk) = k2 * hbarc**2/nuc_mass/2.
                                ! IF ( tzp == -1 ) all_orbit%e(kk) = k2 * hbarc**2/p_mass/2.
                                ! IF ( tzp ==  1 ) all_orbit%e(kk) = k2 * hbarc**2/n_mass/2.
                                ! write(6,*) kk, (all_orbit%kx(kk)**2 + all_orbit%ky(kk)**2 + all_orbit%kz(kk)**2), hbarc, p_mass
                                ! IF ( iam == 0 ) write(6,*) 'state:', ix*nx, iy*ny, iz*nz, sz, tzp, nx**2 + ny**2 + nz**2

                             end DO
                          end DO
                       end DO
                    end DO
                 end DO
              end DO
           end DO
        end DO
     end DO
  end DO

  allocate( check_orbital(all_orbit%total_orbits) )
  DO i = 1, all_orbit%total_orbits
     check_orbital(i) = 1 
  end DO
  
  count_n = 0 
  count_p = 0
  DO i = 1, all_orbit%total_orbits
     N2 = all_orbit%nx(i)**2 + all_orbit%ny(i)**2 + all_orbit%nz(i)**2
     ! neutron
     IF ( all_orbit%tz(i) == 1 ) then
        count_n = count_n + 1
        IF ( count_n <= below_ef_n ) then
           check_orbital(i) = 0
           NF2(1) = max(NF2(1), N2)
        end IF
     end IF
     ! proton
     IF ( all_orbit%tz(i) == -1 ) then
        count_p = count_p + 1
        IF ( count_p <= below_ef_p ) then
           check_orbital(i) = 0
           NF2(-1) = max(NF2(-1), N2)
        end IF
     end IF
  end DO

  DO i = 1, all_orbit%total_orbits
     IF ( check_orbital(i) == 0 ) all_orbit%orbit_status(i) = 'hole'
     IF ( check_orbital(i) == 1 ) all_orbit%orbit_status(i) = 'particle'
  end DO

  ! Sort orbits: holes first (1..below_ef), particles after (below_ef+1..tot_orbs).
  ! This is required because the entire code assumes this index partitioning.
  ! O(N) two-pass sort via permutation array.
  ALLOCATE( order(all_orbit%total_orbits) )
  ih = 0
  ip = below_ef
  DO i = 1, all_orbit%total_orbits
     IF ( check_orbital(i) == 0 ) THEN
        ih = ih + 1
        order(ih) = i
     ELSE
        ip = ip + 1
        order(ip) = i
     END IF
  END DO
  deallocate( check_orbital )

  IF ( ih /= below_ef ) THEN
     IF ( iam == 0 ) write(6,*) 'FATAL: hole count mismatch in orbit sort:', ih, below_ef
     CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierror)
  END IF

  ! Apply permutation via temporary copy
  NULLIFY(tmp_orbit%nx, tmp_orbit%ny, tmp_orbit%nz)
  NULLIFY(tmp_orbit%kx, tmp_orbit%ky, tmp_orbit%kz)
  NULLIFY(tmp_orbit%sz, tmp_orbit%tz, tmp_orbit%e)
  NULLIFY(tmp_orbit%orbit_status, tmp_orbit%orb_type, tmp_orbit%model_space)
  CALL allocate_sp_array(tmp_orbit, all_orbit%total_orbits)
  DO i = 1, all_orbit%total_orbits
     kk = order(i)
     tmp_orbit%kx(i) = all_orbit%kx(kk)
     tmp_orbit%ky(i) = all_orbit%ky(kk)
     tmp_orbit%kz(i) = all_orbit%kz(kk)
     tmp_orbit%nx(i) = all_orbit%nx(kk)
     tmp_orbit%ny(i) = all_orbit%ny(kk)
     tmp_orbit%nz(i) = all_orbit%nz(kk)
     tmp_orbit%sz(i) = all_orbit%sz(kk)
     tmp_orbit%tz(i) = all_orbit%tz(kk)
     tmp_orbit%e(i)  = all_orbit%e(kk)
     tmp_orbit%orbit_status(i) = all_orbit%orbit_status(kk)
  END DO
  DO i = 1, all_orbit%total_orbits
     all_orbit%kx(i) = tmp_orbit%kx(i)
     all_orbit%ky(i) = tmp_orbit%ky(i)
     all_orbit%kz(i) = tmp_orbit%kz(i)
     all_orbit%nx(i) = tmp_orbit%nx(i)
     all_orbit%ny(i) = tmp_orbit%ny(i)
     all_orbit%nz(i) = tmp_orbit%nz(i)
     all_orbit%sz(i) = tmp_orbit%sz(i)
     all_orbit%tz(i) = tmp_orbit%tz(i)
     all_orbit%e(i)  = tmp_orbit%e(i)
     all_orbit%orbit_status(i) = tmp_orbit%orbit_status(i)
  END DO
  CALL deallocate_sp_array(tmp_orbit)
  DEALLOCATE( order )

  ! Validate: orbits 1..below_ef must be holes
  DO i = 1, below_ef
     IF ( all_orbit%orbit_status(i) /= 'hole' ) THEN
        IF ( iam == 0 ) write(6,*) 'FATAL: orbit', i, 'is not a hole but has index <= below_ef'
        CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierror)
     END IF
  END DO
    
  IF ( iam == 0 ) write(6,*) 'Occupied orbits'
  sum_tz = 0
  sum_sz = 0
  DO i = 1, below_ef
     IF ( iam == 0 ) write(6,'(1(i6,1x),1x,1(g20.10,1x),1x,5(I3,1x),a)')  &
          i, all_orbit%e(i), (all_orbit%nx(i)), (all_orbit%ny(i)), (all_orbit%nz(i)), &
          (all_orbit%sz(i)), (all_orbit%tz(i)), all_orbit%orbit_status(i)
     sum_sz = sum_sz + all_orbit%sz(i)
     sum_tz = sum_tz + all_orbit%tz(i)
  end DO
  IF ( iam == 0 ) write(6,*) 'Sz and Tz for vacuum state', sum_sz, sum_tz

  IF ( iam == 0 ) write(6,*) 'Unoccupied orbits'
  DO i = below_ef+1, tot_orbs
     IF ( iam == 0 ) write(6,'(1(i6,1x),1x,1(g20.10,1x),1x,5(I3,1x),a)')  &
          i, all_orbit%e(i), (all_orbit%nx(i)), (all_orbit%ny(i)), (all_orbit%nz(i)), &
          (all_orbit%sz(i)), (all_orbit%tz(i)), all_orbit%orbit_status(i)
  end DO
  IF ( iam == 0 ) write(6,*)

  basis_built = .TRUE.

END SUBROUTINE setup_sp_data


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE setup_structures
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE build_status
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: configs_2b
  INTEGER :: number_channels, ch, ii,jj, nconfs, n_pr_proc
  INTEGER :: Nx,Ny,Nz,Tz
  INTEGER :: p,q, k1,k2,k3,k4
  INTEGER, ALLOCATABLE :: channel_2b(:,:,:,:)
  
  ! get number of 2b channels
  NNmin = 1000
  NNmax = -1000
  DO p = 1, tot_orbs
     DO q = 1, tot_orbs
        IF ( all_orbit%nx(p) + all_orbit%nx(q) > NNmax ) NNmax = all_orbit%nx(p) + all_orbit%nx(q)
        IF ( all_orbit%nx(p) + all_orbit%nx(q) < NNmin ) NNmin = all_orbit%nx(p) + all_orbit%nx(q)
     end DO
  end DO

  ALLOCATE( channel_2b(NNmin:NNmax, NNmin:NNmax, NNmin:NNmax, tzmin:tzmax) )
  channel_2b = 0
  
  DO p = 1, tot_orbs
     DO q = 1, tot_orbs
        Nx = all_orbit%nx(p) + all_orbit%nx(q)
        Ny = all_orbit%ny(p) + all_orbit%ny(q)
        Nz = all_orbit%nz(p) + all_orbit%nz(q)
        Tz = (all_orbit%tz(p) + all_orbit%tz(q))/2
        channel_2b(Nx,Ny,Nz,Tz) = 1
     end DO
  end DO
  
  number_channels = 0
  DO Nx = NNmin, NNmax
     DO Ny = NNmin, NNmax
        DO Nz = NNmin, NNmax
           DO Tz = tzmin, tzmax, 1
              IF ( channel_2b(Nx,Ny,Nz,Tz) == 0 ) cycle
              number_channels = number_channels + 1
           end DO
        end DO
     end DO
  end DO
  channels_2b%number_confs = number_channels

  
  ! set up 2b channels
  ALLOCATE( channels_2b%config_NxNyNz_Tz(4*number_channels) )  
  ALLOCATE( number_2b(1:3,number_channels) )
  channels_2b%config_NxNyNz_Tz = 0
  number_2b = 0
  CALL mem_register('channels_2b', REAL(number_channels * 7 * 4.d0, dp))
  
  n_pr_proc = int(channels_2b%number_confs / num_procs) + 1
  ii = 0
  DO Nx = NNmin, NNmax
     DO Ny = NNmin, NNmax
        DO Nz = NNmin, NNmax
           DO Tz = tzmin, tzmax, 1
              IF ( channel_2b(Nx,Ny,Nz,Tz) == 0 ) cycle
              ii = ii + 1
              
              IF ( ii <= iam*n_pr_proc .or. (iam+1)*n_pr_proc < ii ) cycle
              ch = ii
              k1 = ch*4 - 3
              k2 = ch*4 - 2
              k3 = ch*4 - 1
              k4 = ch*4
              channels_2b%config_NxNyNz_Tz(k1) = Nx
              channels_2b%config_NxNyNz_Tz(k2) = Ny
              channels_2b%config_NxNyNz_Tz(k3) = Nz
              channels_2b%config_NxNyNz_Tz(k4) = Tz
              DO jj = 1, 3
                 CALL number_2b_confs(Nx,Ny,Nz,Tz,jj,configs_2b)
                 nconfs = configs_2b%number_confs
                 IF ( nconfs <= 0 ) cycle
                 number_2b(jj,ch) = nconfs
              end DO
           end DO
        end DO
     end DO
  end DO
  CALL mpi_allreduce(mpi_in_place,channels_2b%config_NxNyNz_Tz,size(channels_2b%config_NxNyNz_Tz),&
       mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_structures', 'allreduce')
  CALL mpi_allreduce(mpi_in_place,number_2b,size(number_2b),mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_structures', 'allreduce')

  
  ALLOCATE( lookup_2b_configs(1:3,number_channels) )
  
  ALLOCATE( hh_config_2b%ival2(1:below_ef, 1:below_ef) )
  ALLOCATE( hh_channel_2b%ival2(1:below_ef, 1:below_ef) )
  hh_config_2b%ival2 = 0
  hh_channel_2b%ival2 = 0
  CALL mem_register('channels_2b', REAL(2 * (below_ef**2) * 4.d0, dp))

  ALLOCATE( hp_config_2b%ival2(1:below_ef, below_ef+1:tot_orbs) )
  ALLOCATE( hp_channel_2b%ival2(1:below_ef, below_ef+1:tot_orbs) )
  hp_config_2b%ival2 = 0
  hp_channel_2b%ival2 = 0
  CALL mem_register('channels_2b', REAL(2 * (below_ef*(tot_orbs-below_ef)) * 4.d0, dp))

  ALLOCATE( pp_config_2b%ival2(below_ef+1:tot_orbs, below_ef+1:tot_orbs) )
  ALLOCATE( pp_channel_2b%ival2(below_ef+1:tot_orbs, below_ef+1:tot_orbs) )
  pp_config_2b%ival2 = 0
  pp_channel_2b%ival2 = 0
  CALL mem_register('channels_2b', REAL(2 * ((tot_orbs-below_ef)**2) * 4.d0, dp))

  
  ii = 0
  DO Nx = NNmin, NNmax
     DO Ny = NNmin, NNmax
        DO Nz = NNmin, NNmax
           DO Tz = tzmin, tzmax, 1
              IF ( channel_2b(Nx,Ny,Nz,Tz) == 0 ) cycle
              ii = ii + 1
              ch = ii

              DO jj = 1, 3
                 nconfs = number_2b(jj,ch)
                 IF ( nconfs <= 0 ) cycle
                 ALLOCATE( lookup_2b_configs(jj,ch)%ival2(2,nconfs) )
                 lookup_2b_configs(jj,ch)%ival2 = 0
                 CALL mem_register('channels_2b', REAL(2 * nconfs * 4.d0, dp))
              end DO
              
              IF ( ii <= iam*n_pr_proc .or. (iam+1)*n_pr_proc < ii ) cycle
              DO jj = 1, 3
                 nconfs = number_2b(jj,ch)
                 IF ( nconfs <= 0 ) cycle
                 ALLOCATE( configs_2b%config(2*nconfs) )
                 configs_2b%number_confs = nconfs
                 CALL setup_2b_confs(Nx,Ny,Nz,Tz,jj,configs_2b)
                 DO nconfs = 1, configs_2b%number_confs
                    p = configs_2b%config(nconfs*2-1)
                    q = configs_2b%config(nconfs*2)
                    lookup_2b_configs(jj,ch)%ival2(1,nconfs) = p
                    lookup_2b_configs(jj,ch)%ival2(2,nconfs) = q
                    IF ( jj == 1 ) then
                       hh_config_2b%ival2(p,q) = nconfs
                       hh_config_2b%ival2(q,p) = nconfs
                       hh_channel_2b%ival2(p,q) = ch
                       hh_channel_2b%ival2(q,p) = ch
                    else IF ( jj == 2 ) then
                       hp_config_2b%ival2(p,q) = nconfs
                       hp_channel_2b%ival2(p,q) = ch
                    else IF ( jj == 3 ) then
                       pp_config_2b%ival2(p,q) = nconfs
                       pp_config_2b%ival2(q,p) = nconfs
                       pp_channel_2b%ival2(p,q) = ch
                       pp_channel_2b%ival2(q,p) = ch
                    end IF
                 end DO
                 DEALLOCATE( configs_2b%config )
              end DO

           end DO
        end DO
     end DO
  end DO
  DEALLOCATE( channel_2b )
  CALL mpi_allreduce(mpi_in_place,hh_config_2b%ival2,size(hh_config_2b%ival2),mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_structures', 'allreduce')
  CALL mpi_allreduce(mpi_in_place,hp_config_2b%ival2,size(hp_config_2b%ival2),mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_structures', 'allreduce')
  CALL mpi_allreduce(mpi_in_place,pp_config_2b%ival2,size(pp_config_2b%ival2),mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_structures', 'allreduce')
  CALL mpi_allreduce(mpi_in_place,hh_channel_2b%ival2,size(hh_channel_2b%ival2),mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_structures', 'allreduce')
  CALL mpi_allreduce(mpi_in_place,hp_channel_2b%ival2,size(hp_channel_2b%ival2),mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_structures', 'allreduce')
  CALL mpi_allreduce(mpi_in_place,pp_channel_2b%ival2,size(pp_channel_2b%ival2),mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_structures', 'allreduce')
  DO ch = 1, channels_2b%number_confs
     DO jj = 1, 3
        nconfs = number_2b(jj,ch)
        IF ( nconfs <= 0 ) cycle
        CALL mpi_allreduce(mpi_in_place,lookup_2b_configs(jj,ch)%ival2,&
             size(lookup_2b_configs(jj,ch)%ival2),mpi_integer,mpi_sum,mpi_comm_world,ierror)
        CALL check_mpi(ierror, 'setup_structures', 'allreduce')
     end DO
  end DO
  
  ! [old memory print removed - replaced by mem_report]

  CALL setup_proc_mappings_v2b

  channels_2b_built = .TRUE.
  mappings_built = .TRUE.

end SUBROUTINE setup_structures


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE setup_cross_structures
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: configs_2b
  INTEGER :: number_channels, ch, ii,jj, nconfs, n_pr_proc
  INTEGER :: Nx,Ny,Nz,Tz, tzmin2,tzmax2
  INTEGER :: p,q, k1,k2,k3,k4
  INTEGER, ALLOCATABLE :: channel_2b(:,:,:,:)
  
  ! get number of 2bcross channels
  tzmin2 = -1
  tzmax2 = 1
  IF ( tzmin == tzmax ) then
     tzmin2 = 0
     tzmax2 = 0
  end IF
  
  ALLOCATE( channel_2b(NNmin:NNmax, NNmin:NNmax, NNmin:NNmax, tzmin2:tzmax2) )
  channel_2b = 0

  N2max = 0
  DO p = 1, tot_orbs
     DO q = 1, tot_orbs
        Nx = all_orbit%nx(p) - all_orbit%nx(q)
        Ny = all_orbit%ny(p) - all_orbit%ny(q)
        Nz = all_orbit%nz(p) - all_orbit%nz(q)
        Tz = (all_orbit%tz(p) - all_orbit%tz(q))/2
        channel_2b(Nx,Ny,Nz,Tz) = 1
        IF ( Nx**2 + Ny**2 + Nz**2 > N2max ) N2max = Nx**2 + Ny**2 + Nz**2
     end DO
  end DO
  
  number_channels = 0
  DO Nx = NNmin, NNmax
     DO Ny = NNmin, NNmax
        DO Nz = NNmin, NNmax
           DO Tz = tzmin2, tzmax2, 1
              IF ( channel_2b(Nx,Ny,Nz,Tz) == 0 ) cycle
              number_channels = number_channels + 1
           end DO
        end DO
     end DO
  end DO
  channels_2bcross%number_confs = number_channels

  
  ! set up 2bcross channels
  ALLOCATE( channels_2bcross%config_NxNyNz_Tz(4*number_channels) )  
  ALLOCATE( number_2bcross(2:3,number_channels) )
  channels_2bcross%config_NxNyNz_Tz = 0
  number_2bcross = 0
  CALL mem_register('channels_cross', REAL(number_channels * 6 * 4.d0, dp))
  
  n_pr_proc = int(channels_2bcross%number_confs / num_procs) + 1
  ii = 0
  DO Nx = NNmin, NNmax
     DO Ny = NNmin, NNmax
        DO Nz = NNmin, NNmax
           DO Tz = tzmin2, tzmax2, 1
              IF ( channel_2b(Nx,Ny,Nz,Tz) == 0 ) cycle
              ii = ii + 1
              
              IF ( ii <= iam*n_pr_proc .or. (iam+1)*n_pr_proc < ii ) cycle
              ch = ii
              k1 = ch*4 - 3
              k2 = ch*4 - 2
              k3 = ch*4 - 1
              k4 = ch*4
              channels_2bcross%config_NxNyNz_Tz(k1) = Nx
              channels_2bcross%config_NxNyNz_Tz(k2) = Ny
              channels_2bcross%config_NxNyNz_Tz(k3) = Nz
              channels_2bcross%config_NxNyNz_Tz(k4) = Tz
              DO jj = 2, 3
                 CALL number_2bcross_confs(Nx,Ny,Nz,Tz,jj,configs_2b)
                 nconfs = configs_2b%number_confs
                 IF ( nconfs <= 0 ) cycle
                 number_2bcross(jj,ch) = nconfs
              end DO
           end DO
        end DO
     end DO
  end DO
  CALL mpi_allreduce(mpi_in_place,channels_2bcross%config_NxNyNz_Tz,size(channels_2bcross%config_NxNyNz_Tz),&
       mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_cross_structures', 'allreduce')
  CALL mpi_allreduce(mpi_in_place,number_2bcross,size(number_2bcross),mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_cross_structures', 'allreduce')

  
  ALLOCATE( lookup_2bcross_configs(2:3,number_channels) )
  
  ALLOCATE( hp_config_2bcross%ival2(1:below_ef, below_ef+1:tot_orbs) )
  ALLOCATE( hp_channel_2bcross%ival2(1:below_ef, below_ef+1:tot_orbs) )
  hp_config_2bcross%ival2 = 0
  hp_channel_2bcross%ival2 = 0
  CALL mem_register('channels_cross', REAL(2 * (below_ef*(tot_orbs-below_ef)) * 4.d0, dp))

  ALLOCATE( ph_config_2bcross%ival2(below_ef+1:tot_orbs, 1:below_ef) )
  ALLOCATE( ph_channel_2bcross%ival2(below_ef+1:tot_orbs, 1:below_ef) )
  ph_config_2bcross%ival2 = 0
  ph_channel_2bcross%ival2 = 0
  CALL mem_register('channels_cross', REAL(2 * ((tot_orbs-below_ef)*below_ef) * 4.d0, dp))

  
  ii = 0
  DO Nx = NNmin, NNmax
     DO Ny = NNmin, NNmax
        DO Nz = NNmin, NNmax
           DO Tz = tzmin2, tzmax2, 1
              IF ( channel_2b(Nx,Ny,Nz,Tz) == 0 ) cycle
              ii = ii + 1
              ch = ii

              DO jj = 2, 3
                 nconfs = number_2bcross(jj,ch)
                 IF ( nconfs <= 0 ) cycle
                 ALLOCATE( lookup_2bcross_configs(jj,ch)%ival2(2,nconfs) )
                 lookup_2bcross_configs(jj,ch)%ival2 = 0
                 CALL mem_register('channels_cross', REAL(2 * nconfs * 4.d0, dp))
              end DO
              
              IF ( ii <= iam*n_pr_proc .or. (iam+1)*n_pr_proc < ii ) cycle
              DO jj = 2, 3
                 nconfs = number_2bcross(jj,ch)
                 IF ( nconfs <= 0 ) cycle
                 ALLOCATE( configs_2b%config(2*nconfs) )
                 configs_2b%number_confs = nconfs
                 CALL setup_2bcross_confs(Nx,Ny,Nz,Tz,jj,configs_2b)
                 DO nconfs = 1, configs_2b%number_confs
                    p = configs_2b%config(nconfs*2-1)
                    q = configs_2b%config(nconfs*2)
                    lookup_2bcross_configs(jj,ch)%ival2(1,nconfs) = p
                    lookup_2bcross_configs(jj,ch)%ival2(2,nconfs) = q
                    IF ( jj == 2 ) then
                       hp_config_2bcross%ival2(p,q) = nconfs
                       hp_channel_2bcross%ival2(p,q) = ch
                    else IF ( jj == 3 ) then
                       ph_config_2bcross%ival2(p,q) = nconfs
                       ph_channel_2bcross%ival2(p,q) = ch
                    end IF
                 end DO
                 DEALLOCATE( configs_2b%config )
              end DO

           end DO
        end DO
     end DO
  end DO
  DEALLOCATE( channel_2b )
  CALL mpi_allreduce(mpi_in_place,hp_config_2bcross%ival2,size(hp_config_2bcross%ival2),mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_cross_structures', 'allreduce')
  CALL mpi_allreduce(mpi_in_place,ph_config_2bcross%ival2,size(ph_config_2bcross%ival2),mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_cross_structures', 'allreduce')
  CALL mpi_allreduce(mpi_in_place,hp_channel_2bcross%ival2,size(hp_channel_2bcross%ival2),mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_cross_structures', 'allreduce')
  CALL mpi_allreduce(mpi_in_place,ph_channel_2bcross%ival2,size(ph_channel_2bcross%ival2),mpi_integer,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'setup_cross_structures', 'allreduce')
  DO ch = 1, channels_2bcross%number_confs
     DO jj = 2, 3
        nconfs = number_2bcross(jj,ch)
        IF ( nconfs <= 0 ) cycle
        CALL mpi_allreduce(mpi_in_place,lookup_2bcross_configs(jj,ch)%ival2,&
             size(lookup_2bcross_configs(jj,ch)%ival2),mpi_integer,mpi_sum,mpi_comm_world,ierror)
        CALL check_mpi(ierror, 'setup_cross_structures', 'allreduce')
     end DO
  end DO

  ! [old memory print removed - replaced by mem_report]
  
  CALL setup_proc_mappings_v2b_cross
     
end SUBROUTINE setup_cross_structures


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
