

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
  INTEGER :: ii, i0, ind
  INTEGER(i8) :: ndim, ndim2, ndim3
  INTEGER :: eom_nums(8,4)
  REAL(dp) :: eom_results(8,5)
  
  ! EOM parameters
  eom_states = 1
  eom_iterations = 550

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
          (ind > 0 .and. nx_paeom == 0 .and. ny_paeom == 0 .and. nz_paeom == 0 .and. tz_paeom == 0) ) then
        eom_nums(ii,1) = 0
        eom_nums(ii,2) = 0
        eom_nums(ii,3) = 0
        eom_nums(ii,4) = 0
        eom_results(ii,1) = 0.d0
        eom_results(ii,2) = 0.d0
        eom_results(ii,3) = 0.d0
        eom_results(ii,4) = 0.d0
        eom_results(ii,5) = 0.d0
        cycle
     end IF
     
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
     
     IF ( ind < 0 ) then
        ALLOCATE( preom_eigs(eom_states) )
        preom_eigs = 0.d0
        CALL mem_register('preom2', REAL(2 * eom_states * 4.d0, dp))
        CALL mem_register('preom2', REAL(eom_states * 16.d0, dp))
        CALL mem_register('preom2', REAL(eom_iterations * 16.d0, dp))
        CALL mem_register('preom2', REAL(3 * eom_iterations**2 * 16.d0, dp))
        CALL mem_register('preom3', REAL((eom%my_stop-eom%my_start+1)*(eom_iterations+2) * 16.d0, dp))
        IF ( iam == 0 ) write(6,*) '...Setting up PR-EOM amplitudes done!'
        CALL mem_report('PR-EOM total')
     else IF ( ind > 0 ) then
        ALLOCATE( paeom_eigs(eom_states) )
        paeom_eigs = 0.d0
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
     IF ( ind < 0 ) CALL compute_right_preom_states
     IF ( ind > 0 ) CALL compute_right_paeom_states
     CALL mpi_barrier(mpi_comm_world,ierror)

     eom_results(ii,1) = kf**2
     IF ( ind < 0 ) then
        i0 = r1_preom_ind(1)
        eom_nums(ii,1) = nx_preom
        eom_nums(ii,2) = ny_preom
        eom_nums(ii,3) = nz_preom
        eom_nums(ii,4) = tz_preom
        eom_results(ii,2) = (all_orbit%kx(i0)**2 + all_orbit%ky(i0)**2 + all_orbit%kz(i0)**2)
        eom_results(ii,3) = all_orbit%e(i0)
        eom_results(ii,4) = real(fock_mtx(i0,i0))
        eom_results(ii,5) = -1.d0 * real(preom_eigs(1))
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
        eom_results(ii,5) = real(paeom_eigs(1))
        CALL deallocate_paeom
     end IF
  end DO
     
  ! Print Results
  IF ( iam == 0 ) write(6,*) ' -- EOM Energy Band --'
  ii = 0
  DO ind = -4, 4
     IF ( ind == 0 ) cycle
     ii = ii + 1
     IF ( iam == 0 ) then
        write(6,'(i4,4x,f4.2,4x,5f16.8)') below_ef, rho, eom_results(ii,1), eom_results(ii,2), &
             eom_results(ii,3), eom_results(ii,4), eom_results(ii,5)        
        ! write(6,'(5(i4,a2),5f16.8)') ind, ':', eom_nums(ii,1), ',', eom_nums(ii,2), ',', &
        !      eom_nums(ii,3), ',', eom_nums(ii,4), ':', eom_results(ii,1), eom_results(ii,2), &
        !      eom_results(ii,3), eom_results(ii,4), eom_results(ii,5)
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
  INTEGER :: ii, i1,i2
  INTEGER(i8) :: ndim, ndim2, ndim3

  CALL get_eom_ind
  IF ( iam == 0 ) write(6,*)
  IF ( iam == 0 ) write(6,*) '...Computing PA-EOM states...'
  IF ( iam == 0 ) write(6,'(A15,28x,A5,I3,3x,A5,I3,3x,A5,I3,3x,A5,I3,3x,A5,I3)') &
       'PA-EOM state: ','nx =',nx_paeom,'ny =',ny_paeom,'nz =',nz_paeom,'tz =',tz_paeom
  IF ( iam == 0 ) write(6,*)

  ! PA-EOM parameters
  eom_states = 1
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

  ALLOCATE( paeom_eigs(eom_states) )
  paeom_eigs = 0.d0
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
  CALL compute_right_paeom_states
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  ! Normalize and print states
  IF ( iam == 0 ) write(6,*) ' -- PA-EOM States --'
  i1 = r1_paeom_ind(1)
  i2 = r1_paeom_ind(2)
  DO ii = 1, eom_states
     IF ( iam == 0 ) then
        write(6,'(5(i4,a2),2f16.8,a2,6f16.8)') ii, ':', nx_paeom, ',', ny_paeom, ',', nz_paeom, ',', tz_paeom, ':', &
             kf**2, (all_orbit%kx(i1)**2 + all_orbit%ky(i1)**2 + all_orbit%kz(i1)**2), ':', &
             all_orbit%e(i1), real(fock_mtx(i1,i1)), real(fock_mtx(i1,i2)), real(hbar1b_I2(i1,i1)), &
             real(hbar1b_I2(i1,i2)), real(paeom_eigs(ii))
     end IF
  end DO
  IF ( iam == 0 ) write(6,*) ' -------------------'
  
end SUBROUTINE compute_paeom_states


SUBROUTINE compute_right_paeom_states
  USE parallel
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: iterations
  INTEGER(i8) :: ii
  REAL(dp) :: startwtime, endwtime
  REAL(dp) :: sigma, tolerance, min_val
  COMPLEX(dpc) :: norm, sum1, e_temp, e1
  COMPLEX(dpc), ALLOCATABLE :: start_vector(:)
  COMPLEX(dpc), ALLOCATABLE :: h_new(:,:)
  COMPLEX(dpc), ALLOCATABLE :: lanc_r(:,:)
  COMPLEX(dpc), ALLOCATABLE :: eigs_temp(:), h_temp(:,:), vec_temp(:,:)
  INTEGER, ALLOCATABLE :: vec_no(:)
  LOGICAL :: selfconsistency
  INTEGER :: k, p,q, i,j
  REAL(dp) :: rand1
  INTEGER, ALLOCATABLE, DIMENSION(:) :: seed
  
  ! EOM parameters
  tolerance = 1.e-6
  
  paeom_eigs = 0.d0
  
  ALLOCATE( vec_no(2*eom_states) )
  ALLOCATE( h_new(eom_iterations, eom_iterations) )
  ALLOCATE( start_vector(eom%my_start:eom%my_stop) )
  h_new  = 0.d0
  vec_no = 0 
  start_vector = 0.d0

  CALL random_seed(size=k)
  ALLOCATE(seed(k))
  CALL random_seed(get=seed)
  ! Initialize PA-EOM vec
  DO ii = eom%my_start, eom%my_stop
     CALL random_number(rand1)
     start_vector(ii) = ( (-0.5d0 + rand1) )/real(ii)**2
     ! start_vector(ii) = 1.d0/real(ii)**2
     IF ( ii > eom_ndim2 ) start_vector(ii) = 0.d0 ! initialize r3 to zero
  end DO

  ! Normalize vector
  norm = sum( start_vector(:)*start_vector(:) )
  CALL mpi_allreduce(mpi_in_place,norm,1,mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'compute_right_paeom_states', 'allreduce')
  start_vector(:) = start_vector(:)/sqrt(norm)
  
  ALLOCATE( lanc_r(eom%my_start:eom%my_stop, 0:eom_iterations) )  
  lanc_r = 0.d0
  lanc_r(:,0) = start_vector(:)

  sigma = 1000.d0
  e1 = 100.d0
  min_val = 1.e-4
  selfconsistency = .FALSE.
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  iterations = 1
  DO WHILE ( .not.selfconsistency .and. iterations <= eom_iterations-1 )
     p = iterations-1
     
     norm = sum( start_vector(:)*start_vector(:) )
     CALL mpi_allreduce(mpi_in_place,norm,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
     CALL check_mpi(ierror, 'compute_right_paeom_states', 'allreduce')
     ! Calculate new lanczos vector by H|lanc> -> |lanc> 
     CALL mpi_barrier(mpi_comm_world,ierror)
     startwtime = MPI_WTIME()
     CALL populate_paeom(eom%my_start, eom%my_stop, start_vector)
     CALL build_hbar3b_paeom
     CALL paeom_1p
     CALL paeom_2p1h
     IF ( eom_approx > 0 ) then
        CALL paeom_1p_3p2h
        CALL paeom_2p1h_3p2h
        CALL paeom_3p2h
     end IF
     CALL populate_paeom_vec(eom%my_start, eom%my_stop, start_vector)
     CALL mpi_barrier(mpi_comm_world,ierror)
     endwtime = MPI_WTIME()
     IF ( iam == master ) write(6,*) 'Total execution for right Hbar.Veom', endwtime - startwtime

     ! Orthogonalize the work vector and setup the Hessenberg matrix
     DO q = 0, p
        sum1 = sum(conjg(lanc_r(:,q))*start_vector(:))
        CALL mpi_allreduce(sum1, h_new(q+1,p+1), 1, mpi_complex16, mpi_sum, mpi_comm_world, ierror)
        CALL check_mpi(ierror, 'compute_right_paeom_states', 'allreduce')
        start_vector(:) = start_vector(:) - h_new(q+1,p+1)*lanc_r(:,q)
     end DO

     norm = sum( conjg(start_vector(:))*start_vector(:) )
     CALL mpi_allreduce(mpi_in_place, norm, 1, mpi_complex16, mpi_sum, mpi_comm_world, ierror)
     CALL check_mpi(ierror, 'compute_right_paeom_states', 'allreduce')
     
     h_new(p+2,p+1) = sqrt(norm)
     
     start_vector(:) = start_vector(:)/sqrt(norm)
     lanc_r(:,p+1) = start_vector(:)

     ALLOCATE( vec_temp(iterations,iterations) ) 
     ALLOCATE( eigs_temp(iterations) ) 
     ALLOCATE( h_temp(iterations,iterations) )
     vec_temp = 0.d0; eigs_temp = 0.d0

     DO i= 1, iterations
        DO j = 1, iterations
           h_temp(i,j) = dcmplx(h_new(i,j))
        end DO
     end DO

     CALL lapack_diag_cmplx(h_temp, vec_temp, eigs_temp, iterations)
     
     ! Get energies and vectors
     ii = 0
     e_temp = 0.d0
     DO i = 1, iterations
        IF ( abs(eigs_temp(i)) < min_val ) cycle
        ii = ii + 1
        IF ( ii <= eom_states ) then
           vec_no(ii) = i
           e_temp = e_temp + eigs_temp(i)                                                  
        end IF
     end DO
     sigma = abs( e1 - e_temp )
     e1 = e_temp

     ii = 0
     DO i = 1, iterations
        IF ( abs(eigs_temp(i)) < min_val ) cycle
        ii = ii + 1
        IF ( ii <= eom_states ) then                                                                                                    
           IF ( iam == 0 ) then
              write(6,'(i4,2x,2f16.8)') ii, eigs_temp(i)
           end IF
        end IF
     end DO
     IF ( iam == 0 ) write(6,'(a,i5,e16.8)') ' Iterations, Delta_E: ', iterations, sigma
     IF ( iam == 0 ) write(6,*)
     
     iterations = iterations + 1
     IF ( ( iterations == eom_iterations .or. (abs(sigma) < tolerance ) ) .and. &
          iterations > min(eom_iterations-1,12) )  then
        selfconsistency = .TRUE. 
     else
        selfconsistency = .FALSE.
     end IF

     IF ( selfconsistency .or. iterations-1 == eom_ndim ) then
        ii = 0
        DO i = 1, iterations-1
           IF ( abs(eigs_temp(i)) < min_val ) cycle
           ii = ii + 1
           IF ( ii <= eom_states ) then
              paeom_eigs(ii) = eigs_temp(vec_no(ii))
           end IF
        end DO

     end IF
     DEALLOCATE( vec_temp, eigs_temp, h_temp )
  end DO

  DEALLOCATE( vec_no ) 
  DEALLOCATE( h_new ) 
  DEALLOCATE( start_vector,lanc_r )

end SUBROUTINE compute_right_paeom_states


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
  INTEGER :: ii, i1,i2
  INTEGER(i8) :: ndim, ndim2, ndim3

  CALL get_eom_ind
  IF ( iam == 0 ) write(6,*)
  IF ( iam == 0 ) write(6,*) '...Computing PR-EOM states...'
  IF ( iam == 0 ) write(6,'(A15,28x,A5,I3,3x,A5,I3,3x,A5,I3,3x,A5,I3,3x,A5,I3)') &
       'PR-EOM state: ','nx =',nx_preom,'ny =',ny_preom,'nz =',nz_preom,'tz =',tz_preom
  IF ( iam == 0 ) write(6,*)

  ! PA-EOM parameters
  eom_states = 1
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

  ALLOCATE( preom_eigs(eom_states) )
  preom_eigs = 0.d0
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
  CALL compute_right_preom_states
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  ! Normalize and print states
  IF ( iam == 0 ) write(6,*) ' -- PR-EOM States --'
  i1 = r1_preom_ind(1)
  i2 = r1_preom_ind(2)
  DO ii = 1, eom_states
     IF ( iam == 0 ) then
        write(6,'(5(i4,a2),2f16.8,a2,6f16.8)') ii, ':', nx_preom, ',', ny_preom, ',', nz_preom, ',', tz_preom, ':', &
             kf**2, (all_orbit%kx(i1)**2 + all_orbit%ky(i1)**2 + all_orbit%kz(i1)**2), ':', &
             all_orbit%e(i1), real(fock_mtx(i1,i1)), real(fock_mtx(i1,i2)), real(hbar1b_I3(i1,i1)), &
             real(hbar1b_I3(i1,i2)), -1.d0 * real(preom_eigs(ii))
     end IF
  end DO
  IF ( iam == 0 ) write(6,*) ' -------------------'
  
end SUBROUTINE compute_preom_states


SUBROUTINE compute_right_preom_states
  USE parallel
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: iterations
  INTEGER(i8) :: ii
  REAL(dp) :: startwtime, endwtime
  REAL(dp) :: sigma, tolerance, min_val
  COMPLEX(dpc) :: norm, sum1, e_temp, e1
  COMPLEX(dpc), ALLOCATABLE :: start_vector(:)
  COMPLEX(dpc), ALLOCATABLE :: h_new(:,:)
  COMPLEX(dpc), ALLOCATABLE :: lanc_r(:,:)
  COMPLEX(dpc), ALLOCATABLE :: eigs_temp(:), h_temp(:,:), vec_temp(:,:)
  INTEGER, ALLOCATABLE :: vec_no(:)
  LOGICAL :: selfconsistency
  INTEGER :: k, p,q, i,j
  REAL(dp) :: rand1
  INTEGER, ALLOCATABLE, DIMENSION(:) :: seed
  
  ! EOM parameters
  tolerance = 1.e-6
  
  preom_eigs = 0.d0
  
  ALLOCATE( vec_no(2*eom_states) )
  ALLOCATE( h_new(eom_iterations, eom_iterations) )
  ALLOCATE( start_vector(eom%my_start:eom%my_stop) )
  h_new  = 0.d0
  vec_no = 0 
  start_vector = 0.d0

  CALL random_seed(size=k)
  ALLOCATE(seed(k))
  CALL random_seed(get=seed)
  ! Initialize PR-EOM vec
  DO ii = eom%my_start, eom%my_stop
     CALL random_number(rand1)
     start_vector(ii) = ( (-0.5d0 + rand1) )/real(ii)**2
     ! start_vector(ii) = 1.d0/real(ii)**2
     IF ( ii > eom_ndim2 ) start_vector(ii) = 0.d0 ! initialize r3 to zero
  end DO

  ! Normalize vector
  norm = sum( start_vector(:)*start_vector(:) )
  CALL mpi_allreduce(mpi_in_place,norm,1,mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'compute_right_preom_states', 'allreduce')
  start_vector(:) = start_vector(:)/sqrt(norm)
  
  ALLOCATE( lanc_r(eom%my_start:eom%my_stop, 0:eom_iterations) )  
  lanc_r = 0.d0
  lanc_r(:,0) = start_vector(:)

  sigma = 1000.d0
  e1 = 100.d0
  min_val = 1.e-4
  selfconsistency = .FALSE.
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  iterations = 1
  DO WHILE ( .not.selfconsistency .and. iterations <= eom_iterations-1 )
     p = iterations-1
     
     norm = sum( start_vector(:)*start_vector(:) )
     CALL mpi_allreduce(mpi_in_place,norm,1,mpi_real8,mpi_sum,mpi_comm_world,ierror)
     CALL check_mpi(ierror, 'compute_right_preom_states', 'allreduce')
     ! Calculate new lanczos vector by H|lanc> -> |lanc> 
     CALL mpi_barrier(mpi_comm_world,ierror)
     startwtime = MPI_WTIME()
     CALL populate_preom(eom%my_start, eom%my_stop, start_vector)
     CALL build_hbar3b_preom
     CALL preom_1h
     CALL preom_1p2h
     IF ( eom_approx > 0 ) then
        CALL preom_1h_2p3h
        CALL preom_1p2h_2p3h
        CALL preom_2p3h
     end IF
     CALL populate_preom_vec(eom%my_start, eom%my_stop, start_vector)
     CALL mpi_barrier(mpi_comm_world,ierror)
     endwtime = MPI_WTIME()
     IF ( iam == master ) write(6,*) 'Total execution for right Hbar.Veom', endwtime - startwtime

     ! Orthogonalize the work vector and setup the Hessenberg matrix
     DO q = 0, p
        sum1 = sum(conjg(lanc_r(:,q))*start_vector(:))
        CALL mpi_allreduce(sum1, h_new(q+1,p+1), 1, mpi_complex16, mpi_sum, mpi_comm_world, ierror)
        CALL check_mpi(ierror, 'compute_right_preom_states', 'allreduce')
        start_vector(:) = start_vector(:) - h_new(q+1,p+1)*lanc_r(:,q)
     end DO

     norm = sum( conjg(start_vector(:))*start_vector(:) )
     CALL mpi_allreduce(mpi_in_place, norm, 1, mpi_complex16, mpi_sum, mpi_comm_world, ierror)
     CALL check_mpi(ierror, 'compute_right_preom_states', 'allreduce')
     
     h_new(p+2,p+1) = sqrt(norm)
     
     start_vector(:) = start_vector(:)/sqrt(norm)
     lanc_r(:,p+1) = start_vector(:)

     ALLOCATE( vec_temp(iterations,iterations) ) 
     ALLOCATE( eigs_temp(iterations) ) 
     ALLOCATE( h_temp(iterations,iterations) )
     vec_temp = 0.d0; eigs_temp = 0.d0

     DO i= 1, iterations
        DO j = 1, iterations
           h_temp(i,j) = dcmplx(h_new(i,j))
        end DO
     end DO

     CALL lapack_diag_cmplx(h_temp, vec_temp, eigs_temp, iterations)
     
     ! Get energies and vectors
     ii = 0
     e_temp = 0.d0
     DO i = 1, iterations
        IF ( abs(eigs_temp(i)) < min_val ) cycle
        ii = ii + 1
        IF ( ii <= eom_states ) then
           vec_no(ii) = i
           e_temp = e_temp + eigs_temp(i)                                                  
        end IF
     end DO
     sigma = abs( e1 - e_temp )
     e1 = e_temp

     ii = 0
     DO i = 1, iterations
        IF ( abs(eigs_temp(i)) < min_val ) cycle
        ii = ii + 1
        IF ( ii <= eom_states ) then                                                                                                    
           IF ( iam == 0 ) then
              write(6,'(i4,2x,2f16.8)') ii, eigs_temp(i)
           end IF
        end IF
     end DO
     IF ( iam == 0 ) write(6,'(a,i5,e16.8)') ' Iterations, Delta_E: ', iterations, sigma
     IF ( iam == 0 ) write(6,*)
     
     iterations = iterations + 1
     IF ( ( iterations == eom_iterations .or. (abs(sigma) < tolerance ) ) .and. &
          iterations > min(eom_iterations-1,12) )  then
        selfconsistency = .TRUE. 
     else
        selfconsistency = .FALSE.
     end IF

     IF ( selfconsistency .or. iterations-1 == eom_ndim ) then
        ii = 0
        DO i = 1, iterations-1
           IF ( abs(eigs_temp(i)) < min_val ) cycle
           ii = ii + 1
           IF ( ii <= eom_states ) then
              preom_eigs(ii) = eigs_temp(vec_no(ii))
           end IF
        end DO

     end IF
     DEALLOCATE( vec_temp, eigs_temp, h_temp )
  end DO

  DEALLOCATE( vec_no ) 
  DEALLOCATE( h_new ) 
  DEALLOCATE( start_vector,lanc_r )

end SUBROUTINE compute_right_preom_states


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
