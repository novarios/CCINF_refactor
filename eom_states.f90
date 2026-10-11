

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE get_eom_ind
  USE single_particle_orbits
  USE constants
  USE mpi_check

  IMPLICIT NONE
  INTEGER :: i, ind, ind0, add0
  INTEGER :: count, count0
  INTEGER :: nx,ny,nz, tz0
  INTEGER, DIMENSION (45) :: Nmaxlist0
  Nmaxlist0 = (/ 0,1,2,3,4,5,6,6,7,8,9,10,11,12,13,13,14,&
                15,16,17,18,19,20,20,21,22,23,24,24,25,26,&
                26,27,28,29,30,30,32,33,34,35,36,37,38,39 /)

  IF ( add_n > 0 ) then
     tz0 = 1
     count0 = below_ef_n
     add0 = add_n
  else IF ( add_p > 0 ) then
     tz0 = -1
     count0 = below_ef_p
     add0 = add_p
  else IF ( add_n < 0 ) then
     tz0 = 1
     count0 = tot_orbs_n-below_ef_n
     add0 = add_n
  else
     tz0 = -1
     count0 = tot_orbs_p-below_ef_p
     add0 = add_p
  end IF

  ! default
  nx_paeom = 0
  ny_paeom = 0
  nz_paeom = 0
  tz_paeom = 0

  ! neutron/proton attachement
  IF ( add0 > 0 ) then
     count = 0
     ind0 = 0
     DO i = 1, all_orbit%total_orbits
        nx = all_orbit%nx(i)
        ny = all_orbit%ny(i)
        nz = all_orbit%nz(i)
        ind = nx**2 + ny**2 + nz**2
        IF ( all_orbit%tz(i) == tz0 ) then
           count = count + 1
           IF ( count <= count0 ) then
              ind0 = nx**2 + ny**2 + nz**2
           else IF ( Nmaxlist0(ind+1) - Nmaxlist0(ind0+1) >= add0 ) then
              nx_paeom = nx
              ny_paeom = ny
              nz_paeom = nz
              tz_paeom = tz0
              return
           end IF
        end IF
     end DO
  end IF

  ! neutron/proton removal
  IF ( add0 < 0 ) then
     count = 0
     ind0 = 0
     DO i = all_orbit%total_orbits, 1, -1
        nx = all_orbit%nx(i)
        ny = all_orbit%ny(i)
        nz = all_orbit%nz(i)
        ind = nx**2 + ny**2 + nz**2
        IF ( all_orbit%tz(i) == tz0 ) then
           count = count + 1
           IF ( count <= count0 ) then
              ind0 = nx**2 + ny**2 + nz**2
           else IF ( Nmaxlist0(ind+1) - Nmaxlist0(ind0+1) <= add0 ) then
              nx_preom = nx
              ny_preom = ny
              nz_preom = nz
              tz_preom = tz0
              return
           end IF
        end IF
     end DO
  end IF

end SUBROUTINE get_eom_ind


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE compute_eom_states
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mem_tracker
  USE mpi_check

  IMPLICIT NONE
  INTEGER :: ii, i0, ind, isel
  INTEGER(i8) :: ndim, ndim2, ndim3
  INTEGER :: eom_nums(8,4), eom_sel(8)
  REAL(dp) :: eom_results(8,6)

  ! EOM parameters
  eom_iterations = 550

  eom_nums = 0
  eom_sel = 0
  eom_results = 0.d0
  ii = 0
  DO ind = -4, 4
     IF ( ind == 0 ) cycle
     ii = ii + 1
     IF ( add_n /= 0 ) add_n = ind
     IF ( add_p /= 0 ) add_p = ind
     CALL get_eom_ind

     IF ( iam == 0 ) write(6,*)
     IF ( iam == 0 .and. ind < 0 ) then
        write(6,*) '...Computing PR-EOM states...'
        write(6,'(A15,28x,A5,I3,3x,A5,I3,3x,A5,I3,3x,A5,I3,3x,A5,I3)') &
             'PR-EOM state: ','nx =',nx_preom,'ny =',ny_preom,'nz =',nz_preom,'tz =',tz_preom
     end IF
     IF ( iam == 0 .and. ind > 0 ) then
        write(6,*) '...Computing PA-EOM states...'
        write(6,'(A15,28x,A5,I3,3x,A5,I3,3x,A5,I3,3x,A5,I3,3x,A5,I3)') &
             'PA-EOM state: ','nx =',nx_paeom,'ny =',ny_paeom,'nz =',nz_paeom,'tz =',tz_paeom
     end IF
     IF ( iam == 0 ) write(6,*)

     IF ( (ind < 0 .and. nx_preom == 0 .and. ny_preom == 0 .and. nz_preom == 0 .and. tz_preom == 0) .or. &
          (ind > 0 .and. nx_paeom == 0 .and. ny_paeom == 0 .and. nz_paeom == 0 .and. tz_paeom == 0) ) cycle

     ! Setup amplitudes and get dimensions
     IF ( iam == 0 .and. ind < 0 ) write(6,*) '...Setting up PR-EOM amplitudes...'
     IF ( iam == 0 .and. ind > 0 ) write(6,*) '...Setting up PA-EOM amplitudes...'

     ndim2 = 0
     ndim3 = 0
     IF ( ind < 0 ) CALL setup_preom_amplitudes(ndim2)
     IF ( ind > 0 ) CALL setup_paeom_amplitudes(ndim2)
     ndim = ndim2
     eom_ndim2 = ndim2
     IF ( eom_approx > 0 ) then
        IF ( ind < 0 ) CALL setup_preom3_amplitudes(ndim3)
        IF ( ind > 0 ) CALL setup_paeom3_amplitudes(ndim3)
        ndim = ndim + ndim3
        eom_ndim3 = ndim3
     end IF
     eom_ndim = ndim
     ! Split up EOM vector
     IF ( ind < 0 ) CALL setup_proc_mappings_preomvec
     IF ( ind > 0 ) CALL setup_proc_mappings_paeomvec

     ALLOCATE( eom_eigs(eom_iterations), eom_wts(eom_iterations), eom_res(eom_iterations), eom_str(eom_iterations) )
     IF ( ind < 0 ) then
        CALL mem_register('preom2', REAL(2 * eom_states * 4.d0, dp))
        CALL mem_register('preom2', REAL(eom_states * 16.d0, dp))
        CALL mem_register('preom2', REAL(eom_iterations * 16.d0, dp))
        CALL mem_register('preom2', REAL(3 * eom_iterations**2 * 16.d0, dp))
        CALL mem_register('preom3', REAL((eom%my_stop-eom%my_start+1)*(eom_iterations+2) * 16.d0, dp))
        IF ( iam == 0 ) write(6,*) '...Setting up PR-EOM amplitudes done!'
        CALL mem_report('PR-EOM total')
     else IF ( ind > 0 ) then
        CALL mem_register('paeom2', REAL(2 * eom_states * 4.d0, dp))
        CALL mem_register('paeom2', REAL(eom_states * 16.d0, dp))
        CALL mem_register('paeom2', REAL(eom_iterations * 16.d0, dp))
        CALL mem_register('paeom2', REAL(3 * eom_iterations**2 * 16.d0, dp))
        CALL mem_register('paeom3', REAL((eom%my_stop-eom%my_start+1)*(eom_iterations+2) * 16.d0, dp))
        IF ( iam == 0 ) write(6,*) '...Setting up PA-EOM amplitudes done!'
        CALL mem_report('PA-EOM total')
     end IF

     IF ( iam == 0 ) write(6,*) '...Setting up Hbar structures...'
     IF ( ind < 0 ) then
        CALL allocate_hbar_preom
        CALL build_hbar_preom
     else IF ( ind > 0 ) then
        CALL allocate_hbar_paeom
        CALL build_hbar_paeom
     end IF
     IF ( iam == 0 ) write(6,*) '...Setting up Hbar structures done!'
     CALL mem_report('interactions + Hbar')

     ! Total Memory
     IF ( iam == 0 .and. ind < 0 ) write(6,*) '...Solving PR-EOM Amplitudes Ready!'
     IF ( iam == 0 .and. ind > 0 ) write(6,*) '...Solving PA-EOM Amplitudes Ready!'
     CALL mem_report('total (EOM)')

     IF ( iam == 0 ) write(6,*)
     IF ( iam == 0 .and. ind < 0 ) write(6,*) 'Solving PR-EOM equations...'
     IF ( iam == 0 .and. ind > 0 ) write(6,*) 'Solving PA-EOM equations...'
     IF ( iam == 0 ) write(6,*)
     CALL compute_right_eom_states(ind > 0)
     CALL mpi_barrier(mpi_comm_world,ierror)
     CALL select_eom_root(ind > 0, isel)
     IF ( eom_sf > 0 ) CALL finish_eom_sf(ind > 0)

     eom_results(ii,1) = kf**2
     eom_sel(ii) = isel
     IF ( ind < 0 ) then
        i0 = r1_preom_ind(1)
        eom_nums(ii,1) = nx_preom
        eom_nums(ii,2) = ny_preom
        eom_nums(ii,3) = nz_preom
        eom_nums(ii,4) = tz_preom
        eom_results(ii,2) = (all_orbit%kx(i0)**2 + all_orbit%ky(i0)**2 + all_orbit%kz(i0)**2)
        eom_results(ii,3) = all_orbit%e(i0)
        eom_results(ii,4) = real(fock_mtx(i0,i0))
        IF ( isel > 0 ) eom_results(ii,5) = -1.d0 * real(eom_eigs(isel))
        CALL deallocate_preom
     else IF ( ind > 0 ) then
        i0 = r1_paeom_ind(1)
        eom_nums(ii,1) = nx_paeom
        eom_nums(ii,2) = ny_paeom
        eom_nums(ii,3) = nz_paeom
        eom_nums(ii,4) = tz_paeom
        eom_results(ii,2) = (all_orbit%kx(i0)**2 + all_orbit%ky(i0)**2 + all_orbit%kz(i0)**2)
        eom_results(ii,3) = all_orbit%e(i0)
        eom_results(ii,4) = real(fock_mtx(i0,i0))
        IF ( isel > 0 ) eom_results(ii,5) = real(eom_eigs(isel))
        CALL deallocate_paeom
     end IF
     IF ( isel > 0 ) eom_results(ii,6) = eom_wts(isel)
     IF ( isel > 0 .and. eom_sf > 0 ) eom_results(ii,6) = eom_str(isel)
     DEALLOCATE( eom_eigs, eom_wts, eom_res, eom_str )
  end DO

  ! Print Results: kF^2, k^2, e_k, f_kk, omega_k, selected root, its 1p (1h) weight
  ! (its strength s_i with eom_sf > 0)
  IF ( iam == 0 ) write(6,*) ' -- EOM Energy Band --'
  ii = 0
  DO ind = -4, 4
     IF ( ind == 0 ) cycle
     ii = ii + 1
     IF ( iam == 0 ) then
        write(6,'(i4,4x,f4.2,4x,5f16.8,i6,f12.6)') below_ef, rho, eom_results(ii,1), eom_results(ii,2), &
             eom_results(ii,3), eom_results(ii,4), eom_results(ii,5), eom_sel(ii), eom_results(ii,6)
     end IF
  end DO
  IF ( iam == 0 ) write(6,*) ' -------------------'

end SUBROUTINE compute_eom_states


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE compute_paeom_states
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mem_tracker
  USE mpi_check

  IMPLICIT NONE
  INTEGER :: isel, i1,i2
  INTEGER(i8) :: ndim, ndim2, ndim3

  CALL get_eom_ind
  IF ( iam == 0 ) write(6,*)
  IF ( iam == 0 ) write(6,*) '...Computing PA-EOM states...'
  IF ( iam == 0 ) write(6,'(A15,28x,A5,I3,3x,A5,I3,3x,A5,I3,3x,A5,I3,3x,A5,I3)') &
       'PA-EOM state: ','nx =',nx_paeom,'ny =',ny_paeom,'nz =',nz_paeom,'tz =',tz_paeom
  IF ( iam == 0 ) write(6,*)

  ! PA-EOM parameters
  eom_iterations = 550

  ! Setup amplitudes and get dimensions
  IF ( iam == 0 ) write(6,*) '...Setting up PA-EOM amplitudes...'
  ndim2 = 0
  ndim3 = 0
  CALL setup_paeom_amplitudes(ndim2)
  ndim = ndim2
  eom_ndim2 = ndim2
  IF ( eom_approx > 0 ) then
     CALL setup_paeom3_amplitudes(ndim3)
     ndim = ndim + ndim3
     eom_ndim3 = ndim3
  end IF
  eom_ndim = ndim
  ! Split up PA-EOM vector
  CALL setup_proc_mappings_paeomvec

  ALLOCATE( eom_eigs(eom_iterations), eom_wts(eom_iterations), eom_res(eom_iterations), eom_str(eom_iterations) )
  CALL mem_register('paeom2', REAL(2 * eom_states * 4.d0, dp))
  CALL mem_register('paeom2', REAL(eom_states * 16.d0, dp))
  CALL mem_register('paeom2', REAL(eom_iterations * 16.d0, dp))
  CALL mem_register('paeom2', REAL(3 * eom_iterations**2 * 16.d0, dp))
  CALL mem_register('paeom3', REAL((eom%my_stop-eom%my_start+1)*(eom_iterations+2) * 16.d0, dp))
  IF ( iam == 0 ) write(6,*) '...Setting up PA-EOM amplitudes done!'
  CALL mem_report('PA-EOM total')


  IF ( iam == 0 ) write(6,*) '...Setting up Hbar structures...'
  CALL allocate_hbar_paeom
  CALL build_hbar_paeom
  IF ( iam == 0 ) write(6,*) '...Setting up Hbar structures done!'
  CALL mem_report('interactions + Hbar')


  ! Total Memory
  IF ( iam == 0 ) write(6,*) '...Solving PA-EOM Amplitudes Ready!'
  CALL mem_report('total (EOM)')


  IF ( iam == 0 ) write(6,*)
  IF ( iam == 0 ) write(6,*) 'Solving PA-EOM equations...'
  IF ( iam == 0 ) write(6,*)
  CALL compute_right_eom_states(.true.)
  CALL mpi_barrier(mpi_comm_world,ierror)
  CALL select_eom_root(.true., isel)
  IF ( eom_sf > 0 ) CALL finish_eom_sf(.true.)

  ! Print selected state
  IF ( iam == 0 ) write(6,*) ' -- PA-EOM States --'
  i1 = r1_paeom_ind(1)
  i2 = r1_paeom_ind(2)
  IF ( iam == 0 .and. isel > 0 ) then
     write(6,'(5(i4,a2),2f16.8,a2,6f16.8,f12.6)') isel, ':', nx_paeom, ',', ny_paeom, ',', nz_paeom, ',', tz_paeom, ':', &
          kf**2, (all_orbit%kx(i1)**2 + all_orbit%ky(i1)**2 + all_orbit%kz(i1)**2), ':', &
          all_orbit%e(i1), real(fock_mtx(i1,i1)), real(fock_mtx(i1,i2)), real(hbar1b_I2(i1,i1)), &
          real(hbar1b_I2(i1,i2)), real(eom_eigs(isel)), eom_wts(isel)
  end IF
  IF ( iam == 0 ) write(6,*) ' -------------------'
  DEALLOCATE( eom_eigs, eom_wts, eom_res, eom_str )

end SUBROUTINE compute_paeom_states


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE compute_preom_states
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mem_tracker
  USE mpi_check

  IMPLICIT NONE
  INTEGER :: isel, i1,i2
  INTEGER(i8) :: ndim, ndim2, ndim3

  CALL get_eom_ind
  IF ( iam == 0 ) write(6,*)
  IF ( iam == 0 ) write(6,*) '...Computing PR-EOM states...'
  IF ( iam == 0 ) write(6,'(A15,28x,A5,I3,3x,A5,I3,3x,A5,I3,3x,A5,I3,3x,A5,I3)') &
       'PR-EOM state: ','nx =',nx_preom,'ny =',ny_preom,'nz =',nz_preom,'tz =',tz_preom
  IF ( iam == 0 ) write(6,*)

  ! PR-EOM parameters
  eom_iterations = 550

  ! Setup amplitudes and get dimensions
  IF ( iam == 0 ) write(6,*) '...Setting up PR-EOM amplitudes...'
  ndim2 = 0
  ndim3 = 0
  CALL setup_preom_amplitudes(ndim2)
  ndim = ndim2
  eom_ndim2 = ndim2
  IF ( eom_approx > 0 ) then
     CALL setup_preom3_amplitudes(ndim3)
     ndim = ndim + ndim3
     eom_ndim3 = ndim3
  end IF
  eom_ndim = ndim
  ! Split up PR-EOM vector
  CALL setup_proc_mappings_preomvec

  ALLOCATE( eom_eigs(eom_iterations), eom_wts(eom_iterations), eom_res(eom_iterations), eom_str(eom_iterations) )
  CALL mem_register('preom2', REAL(2 * eom_states * 4.d0, dp))
  CALL mem_register('preom2', REAL(eom_states * 16.d0, dp))
  CALL mem_register('preom2', REAL(eom_iterations * 16.d0, dp))
  CALL mem_register('preom2', REAL(3 * eom_iterations**2 * 16.d0, dp))
  CALL mem_register('preom3', REAL((eom%my_stop-eom%my_start+1)*(eom_iterations+2) * 16.d0, dp))
  IF ( iam == 0 ) write(6,*) '...Setting up PR-EOM amplitudes done!'
  CALL mem_report('PR-EOM total')


  IF ( iam == 0 ) write(6,*) '...Setting up Hbar structures...'
  CALL allocate_hbar_preom
  CALL build_hbar_preom
  IF ( iam == 0 ) write(6,*) '...Setting up Hbar structures done!'
  CALL mem_report('interactions + Hbar')


  ! Total Memory
  IF ( iam == 0 ) write(6,*) '...Solving PR-EOM Amplitudes Ready!'
  CALL mem_report('total (EOM)')


  IF ( iam == 0 ) write(6,*)
  IF ( iam == 0 ) write(6,*) 'Solving PR-EOM equations...'
  IF ( iam == 0 ) write(6,*)
  CALL compute_right_eom_states(.false.)
  CALL mpi_barrier(mpi_comm_world,ierror)
  CALL select_eom_root(.false., isel)
  IF ( eom_sf > 0 ) CALL finish_eom_sf(.false.)

  ! Print selected state
  IF ( iam == 0 ) write(6,*) ' -- PR-EOM States --'
  i1 = r1_preom_ind(1)
  i2 = r1_preom_ind(2)
  IF ( iam == 0 .and. isel > 0 ) then
     write(6,'(5(i4,a2),2f16.8,a2,6f16.8,f12.6)') isel, ':', nx_preom, ',', ny_preom, ',', nz_preom, ',', tz_preom, ':', &
          kf**2, (all_orbit%kx(i1)**2 + all_orbit%ky(i1)**2 + all_orbit%kz(i1)**2), ':', &
          all_orbit%e(i1), real(fock_mtx(i1,i1)), real(fock_mtx(i1,i2)), real(hbar1b_I3(i1,i1)), &
          real(hbar1b_I3(i1,i2)), -1.d0 * real(eom_eigs(isel)), eom_wts(isel)
  end IF
  IF ( iam == 0 ) write(6,*) ' -------------------'
  DEALLOCATE( eom_eigs, eom_wts, eom_res, eom_str )

end SUBROUTINE compute_preom_states


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Arnoldi solver for the right PA-EOM (pa = .true.) or PR-EOM (pa = .false.) states.
! EOM vector entries 1:2 are the 1p (1h) amplitudes for the two spins. Every Ritz pair
! (theta, y = V c) gets its 1p (1h) weight and its residual |Hbar y - theta y| =
! |h(m+1,m) c(m)|/|c|, exact for an orthonormal Krylov basis V. Candidates for
! select_eom_root go to eom_eigs, eom_wts and eom_res.
!  eom_qp = 0: random start vector; stop when the sum of the eom_states lowest roots is
!              stable; the candidates are those roots.
!  eom_qp = 1: start from the 1p (1h) state; stop when the Ritz pair with the largest
!              weight has converged; the candidates are all converged Ritz pairs.
! eom_sf > 0 overrides eom_qp: start from the 1p (1h) state v and give every Ritz pair its
! strength s_i = Y(1,i) x(i), Y x = e1, from the eigenvectors Y of the Hessenberg matrix H_m.
! The (theta_i, s_i) reproduce the moments <v|Hbar^k|v>, k < m, and sum to 1.
!  eom_sf = 1 (sum): stop when the roots converged to res_tol hold sf_frac of the strength;
!              the candidates are those roots.
!  eom_sf = 2 (lit): stop when A(E) = -Im e1^T (E + i*eom_gamma - H_m)^-1 e1 / pi on a fixed
!              grid changes by less than lit_tol (relative to its maximum).
SUBROUTINE compute_right_eom_states(pa)
  USE parallel
  USE constants
  USE operator_storage
  USE mpi_check

  IMPLICIT NONE
  LOGICAL, INTENT(in) :: pa
  INTEGER :: iterations, nroots
  INTEGER(i8) :: ii
  REAL(dp) :: startwtime, endwtime
  REAL(dp) :: sigma, tolerance, min_val
  COMPLEX(dpc) :: norm, e_temp, e1, c1, c2
  COMPLEX(dpc), ALLOCATABLE :: start_vector(:)
  COMPLEX(dpc), ALLOCATABLE :: h_new(:,:)
  COMPLEX(dpc), ALLOCATABLE :: lanc_r(:,:), lanc_r1(:,:)
  COMPLEX(dpc), ALLOCATABLE :: eigs_temp(:), h_temp(:,:), vec_temp(:,:)
  INTEGER, ALLOCATABLE :: vec_no(:)
  LOGICAL :: selfconsistency, exhausted
  INTEGER :: p, i,j, pass, nloc, imax
  COMPLEX(dpc), ALLOCATABLE :: dots(:)
  REAL(dp), ALLOCATABLE :: wts_temp(:), res_temp(:)
  REAL(dp) :: rand1, vnorm
  REAL(dp), PARAMETER :: res_tol = 1.d-5    ! residual (MeV) of the target root, eom_qp = 1
  REAL(dp), PARAMETER :: report_tol = 1.d-3 ! residual (MeV) for reporting other roots
  REAL(dp), PARAMETER :: sf_frac = 0.95d0   ! converged strength fraction, eom_sf = 1
  REAL(dp), PARAMETER :: lit_tol = 1.d-6    ! relative change of A(E), eom_sf = 2
  REAL(dp), PARAMETER :: lit_half = 250.d0  ! grid half-width (MeV) around <v|Hbar|v>, eom_sf = 2
  COMPLEX(dpc), ALLOCATABLE :: str_temp(:), ymat(:,:), xvec(:)
  INTEGER, ALLOCATABLE :: ipiv(:)
  REAL(dp), ALLOCATABLE :: lit_old(:)
  REAL(dp) :: sconv
  INTEGER :: info, ig

  ! EOM parameters
  tolerance = 1.e-6

  eom_eigs = 0.d0
  eom_wts = 0.d0
  eom_res = 0.d0
  eom_nfound = 0

  ALLOCATE( vec_no(eom_states) )
  ALLOCATE( h_new(eom_iterations, eom_iterations) )
  ALLOCATE( start_vector(eom%my_start:eom%my_stop) )
  h_new  = 0.d0
  vec_no = 0
  start_vector = 0.d0

  IF ( eom_sf > 0 ) then
     ALLOCATE( eom_ritz(eom_iterations), eom_ritz_s(eom_iterations), eom_ritz_res(eom_iterations) )
     eom_ritz = 0.d0; eom_ritz_s = 0.d0; eom_ritz_res = 0.d0
     eom_m = 0
     eom_ngrid = 0
  end IF

  ! Initialize EOM vec: the 1p (1h) state with both spins, or random
  DO ii = eom%my_start, eom%my_stop
     IF ( eom_qp == 1 .or. eom_sf > 0 ) then
        IF ( ii <= 2 ) start_vector(ii) = 1.d0
     else
        CALL random_number(rand1)
        start_vector(ii) = ( (-0.5d0 + rand1) )/real(ii)**2
        IF ( ii > eom_ndim2 ) start_vector(ii) = 0.d0 ! initialize r3 to zero
     end IF
  end DO

  ! Normalize vector
  norm = sum( start_vector(:)*start_vector(:) )
  CALL mpi_allreduce(mpi_in_place,norm,1,mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'compute_right_eom_states', 'allreduce')
  start_vector(:) = start_vector(:)/sqrt(norm)

  ALLOCATE( lanc_r(eom%my_start:eom%my_stop, 0:eom_iterations) )
  ALLOCATE( lanc_r1(2, 0:eom_iterations) ) ! 1p (1h) entries of each Krylov vector, on every rank
  ALLOCATE( dots(0:eom_iterations) )
  lanc_r = 0.d0
  lanc_r1 = 0.d0
  dots = 0.d0
  nloc = int(eom%my_stop - eom%my_start + 1)
  lanc_r(:,0) = start_vector(:)
  CALL gather_r1(0)

  sigma = 1000.d0
  e1 = 100.d0
  min_val = 1.e-4
  selfconsistency = .FALSE.
  exhausted = .FALSE.
  CALL mpi_barrier(mpi_comm_world,ierror)

  iterations = 1
  DO WHILE ( .not.selfconsistency .and. iterations <= eom_iterations-1 )
     p = iterations-1

     ! Calculate new Krylov vector by Hbar|lanc> -> |lanc>
     CALL mpi_barrier(mpi_comm_world,ierror)
     startwtime = MPI_WTIME()
     CALL apply_hbar(start_vector)
     CALL mpi_barrier(mpi_comm_world,ierror)
     endwtime = MPI_WTIME()
     IF ( iam == master ) write(6,*) 'Total execution for right Hbar.Veom', endwtime - startwtime

     ! Orthogonalize the work vector against the Krylov basis and fill the Hessenberg
     ! column. Classical Gram-Schmidt applied twice: a single pass loses orthogonality,
     ! and the residual estimate below assumes an orthonormal basis.
     h_new(1:p+1,p+1) = 0.d0
     DO pass = 1, 2
        CALL zgemv('C', nloc, p+1, (1.d0,0.d0), lanc_r, max(1,nloc), start_vector, 1, (0.d0,0.d0), dots, 1)
        CALL mpi_allreduce(mpi_in_place, dots, p+1, mpi_complex16, mpi_sum, mpi_comm_world, ierror)
        CALL check_mpi(ierror, 'compute_right_eom_states', 'allreduce')
        CALL zgemv('N', nloc, p+1, (-1.d0,0.d0), lanc_r, max(1,nloc), dots, 1, (1.d0,0.d0), start_vector, 1)
        h_new(1:p+1,p+1) = h_new(1:p+1,p+1) + dots(0:p)
     end DO

     norm = sum( conjg(start_vector(:))*start_vector(:) )
     CALL mpi_allreduce(mpi_in_place, norm, 1, mpi_complex16, mpi_sum, mpi_comm_world, ierror)
     CALL check_mpi(ierror, 'compute_right_eom_states', 'allreduce')

     h_new(p+2,p+1) = sqrt(norm)

     ! A vanishing residual means the Krylov space is exhausted (only in tiny spaces)
     IF ( abs(h_new(p+2,p+1)) < 1.d-10 ) then
        exhausted = .TRUE.
     else
        start_vector(:) = start_vector(:)/sqrt(norm)
        lanc_r(:,p+1) = start_vector(:)
        CALL gather_r1(p+1)
     end IF

     ALLOCATE( vec_temp(iterations,iterations) )
     ALLOCATE( eigs_temp(iterations) )
     ALLOCATE( h_temp(iterations,iterations) )
     vec_temp = 0.d0; eigs_temp = 0.d0
     h_temp = h_new(1:iterations,1:iterations)

     CALL lapack_diag_cmplx(h_temp, vec_temp, eigs_temp, iterations)

     ! 1p (1h) weight, |sum_j c_j lanc_r(1:2,j-1)|^2/|c|^2, and residual of every Ritz pair
     ALLOCATE( wts_temp(iterations), res_temp(iterations) )
     DO i = 1, iterations
        vnorm = sqrt( sum( abs(vec_temp(:,i))**2 ) )
        c1 = sum( vec_temp(:,i)*lanc_r1(1,0:iterations-1) )
        c2 = sum( vec_temp(:,i)*lanc_r1(2,0:iterations-1) )
        wts_temp(i) = ( abs(c1)**2 + abs(c2)**2 )/vnorm**2
        res_temp(i) = abs( h_new(p+2,p+1)*vec_temp(iterations,i) )/vnorm
     end DO

     ! Ritz strengths s_i = Y(1,i) x(i) with Y x = e1 (start vector = first Krylov vector)
     IF ( eom_sf > 0 ) then
        ALLOCATE( str_temp(iterations), ymat(iterations,iterations), xvec(iterations), ipiv(iterations) )
        ymat = vec_temp
        xvec = 0.d0
        xvec(1) = 1.d0
        CALL zgesv(iterations, 1, ymat, iterations, ipiv, xvec, iterations, info)
        DO i = 1, iterations
           str_temp(i) = vec_temp(1,i)*xvec(i)
        end DO
        eom_m = iterations
        eom_ritz(1:iterations) = eigs_temp
        eom_ritz_s(1:iterations) = str_temp
        eom_ritz_res(1:iterations) = res_temp
        DEALLOCATE( ymat, xvec, ipiv )
     end IF

     IF ( eom_sf == 1 ) then
        ! Spectral function from the converged states
        eom_nfound = 0
        sconv = 0.d0
        DO i = 1, iterations
           IF ( abs(eigs_temp(i)) < min_val ) cycle
           IF ( res_temp(i) > res_tol ) cycle
           sconv = sconv + real(str_temp(i))
           eom_nfound = eom_nfound + 1
           eom_eigs(eom_nfound) = eigs_temp(i)
           eom_wts(eom_nfound)  = wts_temp(i)
           eom_res(eom_nfound)  = res_temp(i)
           eom_str(eom_nfound)  = real(str_temp(i))
        end DO
        sigma = sconv
        IF ( iam == 0 ) write(6,'(a,i5,f12.6,i6)') ' Iterations, converged strength, roots: ', &
             iterations, sconv, eom_nfound
     else IF ( eom_sf == 2 ) then
        ! Lorentzian-smoothed A(E) on a grid fixed at the first iteration
        IF ( p == 0 ) then
           eom_ngrid = 2*nint(lit_half/(0.25d0*eom_gamma)) + 1
           ALLOCATE( eom_grid(eom_ngrid), eom_lit(eom_ngrid), lit_old(eom_ngrid) )
           DO ig = 1, eom_ngrid
              eom_grid(ig) = real(h_new(1,1)) - lit_half + (ig-1)*0.25d0*eom_gamma
           end DO
           lit_old = 0.d0
        end IF
        DO ig = 1, eom_ngrid
           eom_lit(ig) = -aimag( sum( str_temp(:)/(cmplx(eom_grid(ig), eom_gamma, dpc) - eigs_temp(:)) ) )/pi
        end DO
        sigma = maxval( abs(eom_lit - lit_old) )/max( maxval(abs(eom_lit)), 1.d-300 )
        lit_old = eom_lit

        ! Candidates: every converged Ritz pair, plus the one with the largest strength
        imax = 0
        DO i = 1, iterations
           IF ( abs(eigs_temp(i)) < min_val ) cycle
           IF ( imax == 0 ) then
              imax = i
           else IF ( real(str_temp(i)) > real(str_temp(imax)) + 1.d-6 ) then
              imax = i
           end IF
        end DO
        eom_nfound = 0
        DO i = 1, iterations
           IF ( abs(eigs_temp(i)) < min_val ) cycle
           IF ( res_temp(i) > report_tol .and. i /= imax ) cycle
           eom_nfound = eom_nfound + 1
           eom_eigs(eom_nfound) = eigs_temp(i)
           eom_wts(eom_nfound)  = wts_temp(i)
           eom_res(eom_nfound)  = res_temp(i)
           eom_str(eom_nfound)  = real(str_temp(i))
        end DO
        IF ( iam == 0 ) write(6,'(a,i5,e12.4)') ' Iterations, relative change of A(E): ', iterations, sigma
     else IF ( eom_qp == 0 ) then
        ! Keep the eom_states lowest non-zero roots
        nroots = 0
        e_temp = 0.d0
        DO i = 1, iterations
           IF ( abs(eigs_temp(i)) < min_val ) cycle
           IF ( nroots == eom_states ) exit
           nroots = nroots + 1
           vec_no(nroots) = i
           e_temp = e_temp + eigs_temp(i)
        end DO
        sigma = abs( e1 - e_temp )
        e1 = e_temp

        eom_nfound = nroots
        DO i = 1, nroots
           j = vec_no(i)
           eom_eigs(i) = eigs_temp(j)
           eom_wts(i)  = wts_temp(j)
           eom_res(i)  = res_temp(j)
           IF ( iam == 0 ) write(6,'(i4,2x,2f16.8,f12.6,e12.4)') i, eom_eigs(i), eom_wts(i), eom_res(i)
        end DO
        IF ( iam == 0 ) write(6,'(a,i5,e16.8)') ' Iterations, Delta_E: ', iterations, sigma
     else
        ! Target the Ritz pair with the largest weight (lowest one on noise-level ties)
        imax = 0
        DO i = 1, iterations
           IF ( abs(eigs_temp(i)) < min_val ) cycle
           IF ( imax == 0 ) then
              imax = i
           else IF ( wts_temp(i) > wts_temp(imax) + 1.d-6 ) then
              imax = i
           end IF
        end DO

        ! Candidates: every converged Ritz pair, plus the target
        eom_nfound = 0
        DO i = 1, iterations
           IF ( abs(eigs_temp(i)) < min_val ) cycle
           IF ( res_temp(i) > report_tol .and. i /= imax ) cycle
           eom_nfound = eom_nfound + 1
           eom_eigs(eom_nfound) = eigs_temp(i)
           eom_wts(eom_nfound)  = wts_temp(i)
           eom_res(eom_nfound)  = res_temp(i)
        end DO
        sigma = 1000.d0
        IF ( imax > 0 ) sigma = res_temp(imax)
        IF ( iam == 0 .and. imax > 0 ) write(6,'(a,i5,2f16.8,f12.6,e12.4)') &
             ' Iterations, target (E, weight, residual): ', iterations, eigs_temp(imax), wts_temp(imax), sigma
     end IF
     IF ( iam == 0 ) write(6,*)
     DEALLOCATE( vec_temp, eigs_temp, h_temp, wts_temp, res_temp )
     IF ( ALLOCATED(str_temp) ) DEALLOCATE( str_temp )

     iterations = iterations + 1
     IF ( eom_sf == 1 ) then
        selfconsistency = ( iterations == eom_iterations .or. sigma >= sf_frac ) .and. &
             iterations > min(eom_iterations-1,12)
     else IF ( eom_sf == 2 ) then
        selfconsistency = ( iterations == eom_iterations .or. sigma < lit_tol ) .and. &
             iterations > min(eom_iterations-1,12)
     else IF ( eom_qp == 0 ) then
        selfconsistency = ( iterations == eom_iterations .or. abs(sigma) < tolerance ) .and. &
             iterations > min(eom_iterations-1,12)
     else
        selfconsistency = ( iterations == eom_iterations .or. sigma < res_tol ) .and. &
             iterations > min(eom_iterations-1,12)
     end IF
     IF ( exhausted ) selfconsistency = .TRUE.
  end DO

  ! Keep the Hessenberg matrix, with its last subdiagonal element, for the strength output
  IF ( eom_sf > 0 ) then
     ALLOCATE( eom_hess(eom_m+1, eom_m) )
     eom_hess = h_new(1:eom_m+1, 1:eom_m)
  end IF
  IF ( ALLOCATED(lit_old) ) DEALLOCATE( lit_old )

  DEALLOCATE( vec_no )
  DEALLOCATE( h_new )
  DEALLOCATE( start_vector, lanc_r, lanc_r1, dots )

CONTAINS

  ! Hbar applied to an EOM vector, in place
  SUBROUTINE apply_hbar(vec)
    COMPLEX(dpc), INTENT(inout) :: vec(eom%my_start:eom%my_stop)
    IF ( pa ) then
       CALL populate_paeom(eom%my_start, eom%my_stop, vec)
       CALL build_hbar3b_paeom
       CALL paeom_1p
       CALL paeom_2p1h
       IF ( eom_approx > 0 ) then
          CALL paeom_1p_3p2h
          CALL paeom_2p1h_3p2h
          CALL paeom_3p2h
       end IF
       CALL populate_paeom_vec(eom%my_start, eom%my_stop, vec)
    else
       CALL populate_preom(eom%my_start, eom%my_stop, vec)
       CALL build_hbar3b_preom
       CALL preom_1h
       CALL preom_1p2h
       IF ( eom_approx > 0 ) then
          CALL preom_1h_2p3h
          CALL preom_1p2h_2p3h
          CALL preom_2p3h
       end IF
       CALL populate_preom_vec(eom%my_start, eom%my_stop, vec)
    end IF
  END SUBROUTINE apply_hbar


  ! Copy the 1p (1h) entries of Krylov vector pp to every rank
  SUBROUTINE gather_r1(pp)
    INTEGER, INTENT(in) :: pp
    INTEGER(i8) :: s
    DO s = 1, 2
       IF ( eom%my_start <= s .and. s <= eom%my_stop ) lanc_r1(s,pp) = lanc_r(s,pp)
    end DO
    CALL mpi_allreduce(mpi_in_place, lanc_r1(:,pp), 2, mpi_complex16, mpi_sum, mpi_comm_world, ierror)
    CALL check_mpi(ierror, 'compute_right_eom_states', 'allreduce')
  END SUBROUTINE gather_r1

end SUBROUTINE compute_right_eom_states


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Pick the most single-particle-like candidate root (largest 1p/1h weight, the lowest
! root on ties) and print all candidates. PR energies are printed as omega = -E.
SUBROUTINE select_eom_root(pa, isel)
  USE parallel
  USE constants

  IMPLICIT NONE
  LOGICAL, INTENT(in) :: pa
  INTEGER, INTENT(out) :: isel
  INTEGER :: i
  REAL(dp) :: sgn

  ! roots come in degenerate spin pairs, so ignore weight differences at the noise level;
  ! with eom_sf > 0 select on the strength s_i instead of the right-vector weight
  isel = 0
  DO i = 1, eom_nfound
     IF ( isel == 0 ) then
        isel = i
     else IF ( eom_sf > 0 ) then
        IF ( eom_str(i) > eom_str(isel) + 1.d-6 ) isel = i
     else IF ( eom_wts(i) > eom_wts(isel) + 1.d-6 ) then
        isel = i
     end IF
  end DO

  sgn = 1.d0
  IF ( .not. pa ) sgn = -1.d0
  IF ( iam == 0 .and. eom_sf == 0 ) then
     IF ( pa ) write(6,*) ' -- PA-EOM roots: root, omega (re, im), 1p weight, residual --'
     IF ( .not. pa ) write(6,*) ' -- PR-EOM roots: root, omega (re, im), 1h weight, residual --'
     DO i = 1, eom_nfound
        IF ( i == isel ) then
           write(6,'(i4,2x,2f16.8,f12.6,e12.4,a)') i, sgn*eom_eigs(i), eom_wts(i), eom_res(i), '   <- selected'
        else
           write(6,'(i4,2x,2f16.8,f12.6,e12.4)') i, sgn*eom_eigs(i), eom_wts(i), eom_res(i)
        end IF
     end DO
  else IF ( iam == 0 ) then
     IF ( pa ) write(6,*) ' -- PA-EOM roots: root, omega (re, im), strength, 1p weight, residual --'
     IF ( .not. pa ) write(6,*) ' -- PR-EOM roots: root, omega (re, im), strength, 1h weight, residual --'
     DO i = 1, eom_nfound
        IF ( i == isel ) then
           write(6,'(i4,2x,2f16.8,2f12.6,e12.4,a)') i, sgn*eom_eigs(i), eom_str(i), eom_wts(i), eom_res(i), &
                '   <- selected'
        else
           write(6,'(i4,2x,2f16.8,2f12.6,e12.4)') i, sgn*eom_eigs(i), eom_str(i), eom_wts(i), eom_res(i)
        end IF
     end DO
     write(6,'(a,f12.6)') ' Strength in listed roots: ', sum(eom_str(1:eom_nfound))
  end IF

end SUBROUTINE select_eom_root


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Strength-function output for the current target (eom_sf > 0), then free its arrays.
! Prints the LIT peak and FWHM (eom_sf = 2, in omega: E for PA, -E for PR) and writes to
! <cc_file>_eom<2|3>[_c<eom3_cut>]_<sum|lit>[_g<eom_gamma>]_<n|p><add>.h5, created on the first
! call of the run. Datasets per target, tag = PA<n>/PR<n>:
!  tag_info   (int)  m, nx, ny, nz, tz, eom_sf, eom_ngrid
!  tag_params        eom_gamma, sign (+1 PA, -1 PR: omega = sign*E)
!  tag_H_real/imag   Hessenberg matrix h(1:m+1,1:m) in E (Hbar eigenvalue) units
!  tag_ritz          m x 5: Re E, Im E, Re s, Im s, residual
!  tag_lit           ngrid x 2: E, A(E)   (eom_sf = 2)
SUBROUTINE finish_eom_sf(pa)
  USE parallel
  USE constants
  USE hdf5_wrapper
  USE hdf5

  IMPLICIT NONE
  LOGICAL, INTENT(in) :: pa
  INTEGER :: ind, ipk, il, ir, error
  INTEGER(HID_T) :: file_id
  INTEGER :: info(7)
  REAL(dp) :: sgn, amax, wl, wr
  REAL(dp), ALLOCATABLE :: buf(:,:)
  CHARACTER(len=8) :: tag
  CHARACTER(len=200) :: sf_file
  CHARACTER(len=20) :: str

  sgn = 1.d0
  IF ( .not. pa ) sgn = -1.d0
  ind = add_n
  IF ( add_n == 0 ) ind = add_p

  ! LIT peak and full width at half maximum, by linear interpolation on the grid
  IF ( iam == 0 .and. eom_sf == 2 ) then
     ipk = maxloc(eom_lit, 1)
     amax = eom_lit(ipk)
     il = ipk
     DO WHILE ( il > 1 .and. eom_lit(il) > 0.5d0*amax )
        il = il - 1
     end DO
     ir = ipk
     DO WHILE ( ir < eom_ngrid .and. eom_lit(ir) > 0.5d0*amax )
        ir = ir + 1
     end DO
     wl = eom_grid(il)
     IF ( eom_lit(il+1) /= eom_lit(il) ) wl = eom_grid(il) + (0.5d0*amax - eom_lit(il)) * &
          (eom_grid(il+1) - eom_grid(il))/(eom_lit(il+1) - eom_lit(il))
     wr = eom_grid(ir)
     IF ( eom_lit(ir-1) /= eom_lit(ir) ) wr = eom_grid(ir-1) + (0.5d0*amax - eom_lit(ir-1)) * &
          (eom_grid(ir) - eom_grid(ir-1))/(eom_lit(ir) - eom_lit(ir-1))
     write(6,'(a,f8.3,a,3f16.8)') ' LIT (gamma =', eom_gamma, ' MeV) peak omega, height, FWHM: ', &
          sgn*eom_grid(ipk), amax, abs(wr - wl)
  end IF

  IF ( iam == 0 ) then
     ! file name: ground-state name plus every EOM setting that changes the output
     sf_file = cc_file(1:len_trim(cc_file)-3) // '_eom2'
     IF ( eom_approx > 0 ) then
        sf_file = cc_file(1:len_trim(cc_file)-3) // '_eom3'
        IF ( eom3_cut < 1000 ) then
           write(str,'(I0)') eom3_cut
           sf_file = trim(sf_file) // '_c' // trim(str)
        end IF
     end IF
     IF ( eom_sf == 1 ) sf_file = trim(sf_file) // '_sum'
     IF ( eom_sf == 2 ) then
        write(str,'(F0.2)') eom_gamma
        IF ( str(1:1) == '.' ) str = '0' // str
        sf_file = trim(sf_file) // '_lit_g' // trim(str)
     end IF
     write(str,'(I0)') eom_add0
     IF ( add_n /= 0 ) sf_file = trim(sf_file) // '_n' // trim(str) // '.h5'
     IF ( add_n == 0 ) sf_file = trim(sf_file) // '_p' // trim(str) // '.h5'
     IF ( eom_sf_new ) then
        CALL h5fcreate_f(trim(sf_file), H5F_ACC_TRUNC_F, file_id, error)
        eom_sf_new = .FALSE.
     else
        CALL h5fopen_f(trim(sf_file), H5F_ACC_RDWR_F, file_id, error)
     end IF
     IF ( pa ) write(tag,'(a2,i1)') 'PA', abs(ind)
     IF ( .not. pa ) write(tag,'(a2,i1)') 'PR', abs(ind)

     info = 0
     info(1) = eom_m
     IF ( pa ) info(2:5) = (/ nx_paeom, ny_paeom, nz_paeom, tz_paeom /)
     IF ( .not. pa ) info(2:5) = (/ nx_preom, ny_preom, nz_preom, tz_preom /)
     info(6) = eom_sf
     info(7) = eom_ngrid
     error = write_integer_vector(file_id, info, trim(tag)//'_info')
     error = write_double_vector(file_id, (/ eom_gamma, sgn /), trim(tag)//'_params')
     error = write_double_matrix(file_id, real(eom_hess, dp), trim(tag)//'_H_real')
     error = write_double_matrix(file_id, aimag(eom_hess), trim(tag)//'_H_imag')

     ALLOCATE( buf(eom_m, 5) )
     buf(:,1) = real(eom_ritz(1:eom_m))
     buf(:,2) = aimag(eom_ritz(1:eom_m))
     buf(:,3) = real(eom_ritz_s(1:eom_m))
     buf(:,4) = aimag(eom_ritz_s(1:eom_m))
     buf(:,5) = eom_ritz_res(1:eom_m)
     error = write_double_matrix(file_id, buf, trim(tag)//'_ritz')
     DEALLOCATE( buf )
     IF ( eom_sf == 2 ) then
        ALLOCATE( buf(eom_ngrid, 2) )
        buf(:,1) = eom_grid
        buf(:,2) = eom_lit
        error = write_double_matrix(file_id, buf, trim(tag)//'_lit')
        DEALLOCATE( buf )
     end IF
     CALL h5fclose_f(file_id, error)
  end IF

  DEALLOCATE( eom_hess, eom_ritz, eom_ritz_s, eom_ritz_res )
  IF ( ALLOCATED(eom_grid) ) DEALLOCATE( eom_grid, eom_lit )

end SUBROUTINE finish_eom_sf


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
