
SUBROUTINE cc_energy(energy, out)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE contracts
  USE mpi_check
  
  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: energy
  LOGICAL, INTENT(IN) :: out
  INTEGER :: ch, bra_confs,ket_confs
  INTEGER :: bra,ket, bra_min,bra_max
  COMPLEX(dpc) :: sum, t2, v2b

  CALL assert_built('t2',           'cc_energy')
  CALL assert_built('interactions', 'cc_energy')

  sum = 0.d0
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)
     
     IF ( check_my_channel_v2b_pphh(ch) == 0 ) cycle
     bra_min = mapping_v2b_pphh(iam+1,ch,2)
     bra_max = mapping_v2b_pphh(iam+1,ch,3)
     
     DO bra = bra_min, bra_max
        DO ket = 1, ket_confs
           t2   = t2_ccm(ch)%cval(bra,ket)
           v2b  = conjg(v2b_pphh(ch)%cval(bra,ket))
           sum = sum + t2*v2b
        end DO
     end DO
  end DO
  CALL mpi_allreduce(mpi_in_place,sum,1,mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'cc_energy', 'allreduce')
  
  ecorr2 = sum
  ecorr = ecorr2
  energy = e0 + ecorr
  IF ( iam == 0 .and. out ) write(6,*) 'CC energy', real(energy), real(e0), real(ecorr)

  IF ( test == 2 ) CALL build_energy_test
  
END SUBROUTINE cc_energy


SUBROUTINE cc_energy_3b(energy, out)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE mpi_check
  
  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: energy
  LOGICAL, INTENT(IN) :: out
  INTEGER :: ch3, ch1,ch2, ket_confs
  INTEGER :: bra,ket, bra0,ket0, bra_min,bra_max
  INTEGER :: a,b,c,i,j,k, cind1,kind1
  COMPLEX(dpc) :: t3, v3b
  COMPLEX(dpc) :: sum

  sum = 0.d0
  DO ch3 = ch3_min, ch3_max
     DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
        c       = clist_t3(ch3)%ival2(cind1,1)
        ch1     = clist_t3(ch3)%ival2(cind1,2)
        bra_min = mapping_t3(ch3)%ival2(cind1,1)
        bra_max = mapping_t3(ch3)%ival2(cind1,2)
        DO kind1 = 1, klimit_t3(ch3)
           k         = klist_t3(ch3)%ival2(kind1,1)
           ch2       = klist_t3(ch3)%ival2(kind1,2)
           ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
  
           DO bra  = bra_min, bra_max
              bra0 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
              a    = lookup_2b_configs(3,ch1)%ival2(1,bra0)
              b    = lookup_2b_configs(3,ch1)%ival2(2,bra0)
              IF ( c <= b ) cycle
              DO ket  = 1, ket_confs
                 ket0 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
                 i    = lookup_2b_configs(1,ch2)%ival2(1,ket0)
                 j    = lookup_2b_configs(1,ch2)%ival2(2,ket0)
                 IF ( k <= j ) cycle
                 ! t3.W3
                 v3b = t3_ccm0(ch3)%val2(cind1,kind1)%cval(bra,ket)
                 t3  = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket)
                 sum = sum + t3*conjg(v3b)
              end DO
           end DO
        end DO
     end DO
  end DO
  CALL mpi_allreduce(mpi_in_place,sum,1,mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'cc_energy_3b', 'allreduce')
  
  ecorr3 = sum
  ecorr = ecorr2 + ecorr3
  energy = e0 + ecorr
  IF ( iam == 0 .and. out ) write(6,*) '3b correlation energy', real(ecorr3)
  IF ( iam == 0 .and. out ) write(6,*) 'CC energy', real(energy), real(e0), real(ecorr2)
  
  IF ( test == 2 ) CALL build_energy_test
  
END SUBROUTINE cc_energy_3b


SUBROUTINE cc_energy_lambda(energy, out)
  USE parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  USE operator_storage
  USE chiral_potentials
  USE mpi_check
  
  IMPLICIT NONE
  COMPLEX(dpc), INTENT(INOUT) :: energy
  LOGICAL, INTENT(IN) :: out
  INTEGER :: ch3, ch1,ch2, ket_confs
  INTEGER :: bra,ket, bra0,ket0, bra_min,bra_max
  INTEGER :: a,b,c,i,j,k, cind1,kind1
  COMPLEX(dpc) :: t3, v3b, denom
  COMPLEX(dpc) :: sum

  IF ( iam == 0 ) WRITE(6,*) "...Computing T3(Lam)..."

  ! <cab|t|kij> <-- P(c/ab|k/ij).[- <cd|t|ij>.<ab|v|kd> + <ab|t|lk>.<lc|v|ij>]
  CALL t3_eqn

  ! ! <cab|t|kij> <-- + P(c/ab|k/ij).<cd|t|kl>.<lab|w|dij>
  ! CALL t3_w3_eqn1
  ! ! <cab|t|kij> <-- + P(k/ij).(1/2).<de|t|ij>.<cab|w|kde>
  ! CALL t3_w3_eqn2
  ! ! <cab|t|kij> <-- + P(c/ab).(1/2).<ab|t|lm>.<clm|w|kij>
  ! CALL t3_w3_eqn3

  ! P(k/ij).<cab|t|kij>
  CALL antisymmetrize_t3

  DO ch3 = ch3_min, ch3_max
     DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
        c       = clist_t3(ch3)%ival2(cind1,1)
        ch1     = clist_t3(ch3)%ival2(cind1,2)
        bra_min = mapping_t3(ch3)%ival2(cind1,1)
        bra_max = mapping_t3(ch3)%ival2(cind1,2)
        DO kind1 = 1, klimit_t3(ch3)
           k         = klist_t3(ch3)%ival2(kind1,1)
           ch2       = klist_t3(ch3)%ival2(kind1,2)
           ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
           
           !$omp parallel default(shared) private(bra,ket,bra0,ket0, a,b,i,j, denom, v3b)
           !$omp do schedule(dynamic)
           DO bra = bra_min, bra_max
              DO ket = 1, ket_confs
                 bra0 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
                 a    = lookup_2b_configs(3,ch1)%ival2(1,bra0)
                 b    = lookup_2b_configs(3,ch1)%ival2(2,bra0)
                 IF ( a == c .or. b == c ) cycle
                 ket0 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
                 i    = lookup_2b_configs(1,ch2)%ival2(1,ket0)
                 j    = lookup_2b_configs(1,ch2)%ival2(2,ket0)
                 IF ( i == k .or. j == k ) cycle
                 ! v3b = v3int(a,b,c,i,j,k)
                 denom = ( fock_mtx(i,i) + fock_mtx(j,j) + fock_mtx(k,k) &
                      - fock_mtx(a,a) - fock_mtx(b,b) - fock_mtx(c,c) )
                 ! ! <abc|t|ijk> <-- <abc|w|ijk>
                 ! t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket) + v3b
                 ! <abc|t|ijk>/E^{abc}_{ijk}
                 t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket)/denom
              end DO
           end DO
           !$omp end do
           !$omp end parallel
        end DO
     end DO
  end DO
  CALL mpi_barrier(mpi_comm_world,ierror)  
  IF ( test == 3 ) CALL build_t3_eqn_test
  
  sum = 0.d0
  DO ch3 = ch3_min, ch3_max
     DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
        c       = clist_t3(ch3)%ival2(cind1,1)
        ch1     = clist_t3(ch3)%ival2(cind1,2)
        bra_min = mapping_t3(ch3)%ival2(cind1,1)
        bra_max = mapping_t3(ch3)%ival2(cind1,2)
        DO kind1 = 1, klimit_t3(ch3)
           k         = klist_t3(ch3)%ival2(kind1,1)
           ch2       = klist_t3(ch3)%ival2(kind1,2)
           ket_confs = number_2b_t3(ch3)%ival2(2,ch2)

           DO bra  = bra_min, bra_max
              bra0 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
              a    = lookup_2b_configs(3,ch1)%ival2(1,bra0)
              b    = lookup_2b_configs(3,ch1)%ival2(2,bra0)
              IF ( c <= b ) cycle
              DO ket  = 1, ket_confs
                 ket0 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
                 i    = lookup_2b_configs(1,ch2)%ival2(1,ket0)
                 j    = lookup_2b_configs(1,ch2)%ival2(2,ket0)
                 IF ( k <= j ) cycle
                 t3  = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket)
                 denom = ( fock_mtx(i,i) + fock_mtx(j,j) + fock_mtx(k,k) &
                      - fock_mtx(a,a) - fock_mtx(b,b) - fock_mtx(c,c) )
                 sum = sum + t3*conjg(t3)*denom
              end DO
           end DO
        end DO
     end DO
  end DO
  CALL mpi_allreduce(mpi_in_place,sum,1,mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'cc_energy_lambda', 'allreduce')
  
  ecorr3 = sum
  ecorr = ecorr2 + ecorr3
  energy = e0 + ecorr
  IF ( iam == 0 .and. out ) write(6,*) '3b correlation energy', real(ecorr3)
  IF ( iam == 0 .and. out ) write(6,*) 'CC energy', real(energy), real(e0), real(ecorr2)
  
  IF ( test == 2 ) CALL build_energy_test
  
END SUBROUTINE cc_energy_lambda

! <cab|t|kij> <-- P(c/ab).<cd|t|kl>.<lab|w|dij>
SUBROUTINE t3_w3_eqn1
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE chiral_potentials
  USE parallel
  USE omp_lib
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch1,ch2,ch3, ch0
  INTEGER :: bra,ket, bra0,bra1, bra_min,bra_max
  INTEGER :: bra_confs,ket_confs, bra_confs0
  INTEGER :: a,b,c, i,j,k, d,l, cind1,kind1, cind10
  INTEGER :: bra_ab, bra_cb, bra_ac
  INTEGER :: phase_ab, phase_cb, phase_ac
  INTEGER :: dim1,dim2,dim3
  REAL(dp) :: startwtime, endwtime
  COMPLEX(dpc) :: v3b
  INTEGER, ALLOCATABLE :: clist_inv(:)
  COMPLEX(dpc), ALLOCATABLE :: t3(:), t2(:)
  TYPE (superblock_storage), ALLOCATABLE :: w3_temp(:)
  
  IF ( iam == 0 ) WRITE(6,'(A29)',advance='no') " ...Computing T3 <- W3_1...  "
  startwtime = MPI_WTIME()
  CALL mpi_barrier(mpi_comm_world,ierror)
  ALLOCATE( clist_inv(below_ef+1:tot_orbs) )
     
  DO ch3 = ch3_min, ch3_max
     ! Build inverse c map
     clist_inv = 0
     DO cind1 = 1, climit_t3(ch3)
        c = clist_t3(ch3)%ival2(cind1,1)
        clist_inv(c) = cind1
     end DO
     
     DO kind1 = 1, klimit_t3(ch3)
        k         = klist_t3(ch3)%ival2(kind1,1)
        ch2       = klist_t3(ch3)%ival2(kind1,2)
        ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
        
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ALLOCATE( w3_temp(climit_t3(ch3)) )
        DO cind1 = 1, climit_t3(ch3)
           c = clist_t3(ch3)%ival2(cind1,1)
           ch1 = clist_t3(ch3)%ival2(cind1,2)
           bra_confs = number_2b_t3(ch3)%ival2(1,ch1)
           
           ALLOCATE( w3_temp(cind1)%val1(bra_confs) )
           ch0 = ph_channel_2bcross%ival2(c,k)
           bra_confs0 = number_2bcross(2,ch0) ! l-d           
           DO bra = 1, bra_confs
              ALLOCATE( w3_temp(cind1)%val1(bra)%cval(bra_confs0,ket_confs) )
              w3_temp(cind1)%val1(bra)%cval = 0.d0
           end DO
        end DO
        
        DO cind1 = 1, climit_t3(ch3)
           c = clist_t3(ch3)%ival2(cind1,1)
           ch1 = clist_t3(ch3)%ival2(cind1,2)
           bra_confs = number_2b_t3(ch3)%ival2(1,ch1)
           ch0 = ph_channel_2bcross%ival2(c,k)
           bra_confs0 = number_2bcross(2,ch0) ! l-d
           !$omp parallel default(shared) private(bra,bra0,ket, l,d,a,b,i,j, v3b)
           !$omp do schedule(dynamic)
           DO bra0 = 1, bra_confs0
              DO bra = 1, bra_confs
                 DO ket = 1, ket_confs
                    l   = lookup_2bcross_configs(2,ch0)%ival2(1,bra0)
                    d   = lookup_2bcross_configs(2,ch0)%ival2(2,bra0)
                    bra1 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)                 
                    a   = lookup_2b_configs(3,ch1)%ival2(1,bra1)
                    b   = lookup_2b_configs(3,ch1)%ival2(2,bra1)
                    IF ( a == c .or. b == c ) cycle
                    i   = lookup_2b_configs(1,ch2)%ival2(1,ket)
                    j   = lookup_2b_configs(1,ch2)%ival2(2,ket)
                    IF ( i == k .or. j == k ) cycle
                    v3b = v3int(l,a,b,d,i,j)
                    w3_temp(cind1)%val1(bra)%cval(bra0,ket) = v3b
                 end DO
              end DO
           end DO
           !$omp end do
           !$omp end parallel
        end DO
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)           
           c       = clist_t3(ch3)%ival2(cind1,1)
           ch1     = clist_t3(ch3)%ival2(cind1,2)
           bra_min = mapping_t3(ch3)%ival2(cind1,1)
           bra_max = mapping_t3(ch3)%ival2(cind1,2)

           DO bra  = bra_min, bra_max
              bra1 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
              a    = lookup_2b_configs(3,ch1)%ival2(1,bra1)
              b    = lookup_2b_configs(3,ch1)%ival2(2,bra1)
              IF ( a == c .or. b == c ) cycle
              phase_ab = -1
              bra_ab = bra
              cind10 = clist_inv(c)
              ch0 = ph_channel_2bcross%ival2(c,k)
              bra0 = ph_config_2bcross%ival2(c,k)

              dim1 = 1
              dim2 = ket_confs
              dim3 = number_2bcross(2,ch0)
              IF ( dim3 == 0 ) cycle
              ALLOCATE( t2(dim3) )
              t2(:) = t2_ccm_cross(ch0)%cval(bra0,:)
              ALLOCATE( t3(dim2) )
              t3(:) = 0.d0
              call ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(real(phase_ab),0.d0), t2, dim1, &
                   w3_temp(cind10)%val1(bra_ab)%cval, dim3, dcmplx(1.d0,0.d0), t3, dim1 )
              t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,:) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,:) + t3(:)
              DEALLOCATE( t2,t3 )
           end DO
           
           DO bra  = bra_min, bra_max
              bra1 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
              a    = lookup_2b_configs(3,ch1)%ival2(1,bra1)
              b    = lookup_2b_configs(3,ch1)%ival2(2,bra1)
              IF ( a == c .or. b == c ) cycle
              phase_cb = 1
              bra_cb = pp_config_2b%ival2(c,b)
              IF ( b < c ) phase_cb = -phase_cb
              cind10 = clist_inv(a)
              ch0 = ph_channel_2bcross%ival2(a,k)
              bra0 = ph_config_2bcross%ival2(a,k)

              dim1 = 1
              dim2 = ket_confs
              dim3 = number_2bcross(2,ch0)
              IF ( dim3 == 0 ) cycle
              ALLOCATE( t2(dim3) )
              t2(:) = t2_ccm_cross(ch0)%cval(bra0,:)
              ALLOCATE( t3(dim2) )
              t3(:) = 0.d0
              call ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(real(phase_cb),0.d0), t2, dim1, &
                   w3_temp(cind10)%val1(bra_cb)%cval, dim3, dcmplx(1.d0,0.d0), t3, dim1 )
              t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,:) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,:) + t3(:)
              DEALLOCATE( t2,t3 )
           end DO

           DO bra  = bra_min, bra_max
              bra1 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
              a    = lookup_2b_configs(3,ch1)%ival2(1,bra1)
              b    = lookup_2b_configs(3,ch1)%ival2(2,bra1)
              IF ( a == c .or. b == c ) cycle
              phase_ac = 1
              bra_ac = pp_config_2b%ival2(a,c)
              IF ( c < a ) phase_ac = -phase_ac
              cind10 = clist_inv(b)
              ch0 = ph_channel_2bcross%ival2(b,k)
              bra0 = ph_config_2bcross%ival2(b,k)

              dim1 = 1
              dim2 = ket_confs
              dim3 = number_2bcross(2,ch0)
              IF ( dim3 == 0 ) cycle
              ALLOCATE( t2(dim3) )
              t2(:) = t2_ccm_cross(ch0)%cval(bra0,:)
              ALLOCATE( t3(dim2) )
              t3(:) = 0.d0
              call ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(real(phase_ac),0.d0), t2, dim1, &
                   w3_temp(cind10)%val1(bra_ac)%cval, dim3, dcmplx(1.d0,0.d0), t3, dim1 )
              t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,:) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,:) + t3(:)
              DEALLOCATE( t2,t3 )
           end DO

        end DO

        DO cind1 = 1, climit_t3(ch3)
           ch1 = clist_t3(ch3)%ival2(cind1,2)
           bra_confs = number_2b_t3(ch3)%ival2(1,ch1)
           DO bra = 1, bra_confs
              DEALLOCATE( w3_temp(cind1)%val1(bra)%cval )
           end DO
           DEALLOCATE( w3_temp(cind1)%val1 )
        end DO
        DEALLOCATE( w3_temp )

     end DO
  end DO
  DEALLOCATE( clist_inv )
  
  CALL mpi_barrier(mpi_comm_world,ierror)  
  endwtime = MPI_WTIME()
  IF ( iam == 0 ) WRITE(6,'(A29,A29,f12.6,A5)') bsp29, " ...Computing T3 <- W3_1...  ", (endwtime - startwtime), ' sec.'

end SUBROUTINE t3_w3_eqn1

! <cab|t|kij> <-- (1/2).<de|t|ij>.<cab|w|kde>
SUBROUTINE t3_w3_eqn2
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE chiral_potentials
  USE parallel
  USE omp_lib
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch1,ch2,ch3
  INTEGER :: bra,bra1, ket0, bra_min,bra_max
  INTEGER :: ket_confs, ket_confs0
  INTEGER :: a,b,c, k, d,e, cind1,kind1
  INTEGER :: dim1,dim2,dim3
  REAL(dp) :: startwtime, endwtime
  COMPLEX(dpc) :: v3b
  COMPLEX(dpc), ALLOCATABLE :: w3_temp(:,:)
  
  IF ( iam == 0 ) WRITE(6,'(A29)',advance='no') " ...Computing T3 <- W3_2...  "
  startwtime = MPI_WTIME()
  CALL mpi_barrier(mpi_comm_world,ierror)
  
  DO ch3 = ch3_min, ch3_max

     DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
        c       = clist_t3(ch3)%ival2(cind1,1)
        ch1     = clist_t3(ch3)%ival2(cind1,2)
        bra_min = mapping_t3(ch3)%ival2(cind1,1)
        bra_max = mapping_t3(ch3)%ival2(cind1,2)
        DO kind1 = 1, klimit_t3(ch3)
           k = klist_t3(ch3)%ival2(kind1,1)
           ch2 = klist_t3(ch3)%ival2(kind1,2)
           ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
           ket_confs0 = number_2b(3,ch2)

           ALLOCATE( w3_temp(bra_min:bra_max, ket_confs0) )
           w3_temp = 0.d0
           
           !$omp parallel default(shared) private(bra,ket0, a,b,d,e, v3b)
           !$omp do schedule(dynamic)
           DO bra = bra_min, bra_max
              DO ket0 = 1, ket_confs0
                 bra1 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
                 a    = lookup_2b_configs(3,ch1)%ival2(1,bra1)
                 b    = lookup_2b_configs(3,ch1)%ival2(2,bra1)
                 IF ( a == c .or. b == c ) cycle
                 d    = lookup_2b_configs(3,ch2)%ival2(1,ket0)
                 e    = lookup_2b_configs(3,ch2)%ival2(2,ket0)
                 v3b = v3int(c,a,b,k,d,e)
                 w3_temp(bra,ket0) = v3b
              end DO
           end DO
           !$omp end do
           !$omp end parallel

           dim1 = bra_max - bra_min + 1
           dim2 = ket_confs
           dim3 = ket_confs0
           CALL ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), w3_temp(bra_min:bra_max,:), dim1, &
                t2_ccm(ch2)%cval, dim3, dcmplx(1.d0,0.d0), t3_ccm(ch3)%val2(cind1,kind1)%cval(bra_min:bra_max,:), dim1 )
           
           DEALLOCATE( w3_temp )           
        end DO
     end DO
     
  end DO

  CALL mpi_barrier(mpi_comm_world,ierror)  
  endwtime = MPI_WTIME()
  IF ( iam == 0 ) WRITE(6,'(A29,A29,f12.6,A5)') bsp29, " ...Computing T3 <- W3_2...  ", (endwtime - startwtime), ' sec.'

end SUBROUTINE t3_w3_eqn2

! <cab|t|kij> <-- (1/2).<ab|t|lm>.<clm|w|kij>
SUBROUTINE t3_w3_eqn3
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE chiral_potentials
  USE parallel
  USE omp_lib
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch1,ch2,ch3
  INTEGER :: bra,ket, bra0,bra1,ket1, bra_min,bra_max
  INTEGER :: ket_confs, bra_confs0
  INTEGER :: a,b,c, i,j,k, l,m, cind1,kind1, cind10
  INTEGER :: ch_ab, ch_cb, ch_ac, bra_ab, bra_cb, bra_ac
  INTEGER :: phase_ab, phase_cb, phase_ac
  INTEGER :: dim1,dim2,dim3
  REAL(dp) :: startwtime, endwtime
  COMPLEX(dpc) :: v3b
  INTEGER, ALLOCATABLE :: clist_inv(:)
  COMPLEX(dpc), ALLOCATABLE :: t3(:), t2(:)
  TYPE (superblock_storage) :: w3_temp
  
  IF ( iam == 0 ) WRITE(6,'(A29)',advance='no') " ...Computing T3 <- W3_3...  "
  startwtime = MPI_WTIME()
  CALL mpi_barrier(mpi_comm_world,ierror)
  ALLOCATE( clist_inv(below_ef+1:tot_orbs) )
     
  DO ch3 = ch3_min, ch3_max
     ! Build inverse c map
     clist_inv = 0
     DO cind1 = 1, climit_t3(ch3)
        c = clist_t3(ch3)%ival2(cind1,1)
        clist_inv(c) = cind1
     end DO

     DO kind1 = 1, klimit_t3(ch3)
        k         = klist_t3(ch3)%ival2(kind1,1)
        ch2       = klist_t3(ch3)%ival2(kind1,2)
        ket_confs = number_2b_t3(ch3)%ival2(2,ch2)

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ALLOCATE( w3_temp%val1(climit_t3(ch3)) )
        DO cind1 = 1, climit_t3(ch3)
           c = clist_t3(ch3)%ival2(cind1,1)
           ch1 = clist_t3(ch3)%ival2(cind1,2)
           bra_confs0 = number_2b(1,ch1)
           ALLOCATE( w3_temp%val1(cind1)%cval(bra_confs0,ket_confs) )
           w3_temp%val1(cind1)%cval = 0.d0
        end DO

        DO cind1 = 1, climit_t3(ch3)
           c = clist_t3(ch3)%ival2(cind1,1)
           ch1 = clist_t3(ch3)%ival2(cind1,2)
           bra_confs0 = number_2b(1,ch1)
           !$omp parallel default(shared) private(bra0,ket, l,m,i,j, v3b)
           !$omp do schedule(static) collapse(2)
           DO bra0 = 1, bra_confs0
              DO ket  = 1, ket_confs
                 ket1 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
                 i    = lookup_2b_configs(1,ch2)%ival2(1,ket1)
                 j    = lookup_2b_configs(1,ch2)%ival2(2,ket1)
                 IF ( i == k .or. j == k ) cycle
                 l   = lookup_2b_configs(1,ch1)%ival2(1,bra0)
                 m   = lookup_2b_configs(1,ch1)%ival2(2,bra0)
                 v3b = v3int(c,l,m,k,i,j)
                 w3_temp%val1(cind1)%cval(bra0,ket) = v3b/3.d0
              end DO
           end DO
           !$omp end do
           !$omp end parallel
        end DO
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
           c       = clist_t3(ch3)%ival2(cind1,1)
           ch1     = clist_t3(ch3)%ival2(cind1,2)
           bra_min = mapping_t3(ch3)%ival2(cind1,1)
           bra_max = mapping_t3(ch3)%ival2(cind1,2)

           DO bra  = bra_min, bra_max
              bra1 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
              a    = lookup_2b_configs(3,ch1)%ival2(1,bra1)
              b    = lookup_2b_configs(3,ch1)%ival2(2,bra1)
              IF ( a == c .or. b == c ) cycle
              cind10 = clist_inv(c)
              ch_ab = ch1
              phase_ab = 1
              bra_ab = bra
              
              dim1 = 1
              dim2 = ket_confs
              dim3 = number_2b(1,ch_ab)
              IF ( dim3 == 0 ) cycle
              ALLOCATE( t2(dim3) )
              t2(:) = t2_ccm(ch_ab)%cval(bra_ab,:)
              ALLOCATE( t3(dim2) )
              t3(:) = 0.d0
              CALL ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(real(phase_ab),0.d0), t2, dim1, &
                   w3_temp%val1(cind10)%cval, dim3, dcmplx(1.d0,0.d0), t3, dim1 )
              t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,:) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,:) + t3(:)
              DEALLOCATE( t2,t3 )
           end DO

           DO bra  = bra_min, bra_max
              bra1 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
              a    = lookup_2b_configs(3,ch1)%ival2(1,bra1)
              b    = lookup_2b_configs(3,ch1)%ival2(2,bra1)
              IF ( a == c .or. b == c ) cycle
              cind10 = clist_inv(a)
              ch_cb = pp_channel_2b%ival2(c,b)
              phase_cb = -1
              bra_cb = pp_config_2b%ival2(c,b)
              IF ( b < c ) phase_cb = -phase_cb
              
              dim1 = 1
              dim2 = ket_confs
              dim3 = number_2b(1,ch_cb)
              IF ( dim3 == 0 ) cycle
              ALLOCATE( t2(dim3) )
              t2(:) = t2_ccm(ch_cb)%cval(bra_cb,:)
              ALLOCATE( t3(dim2) )
              t3(:) = 0.d0
              CALL ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(real(phase_cb),0.d0), t2, dim1, &
                   w3_temp%val1(cind10)%cval, dim3, dcmplx(1.d0,0.d0), t3, dim1 )
              t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,:) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,:) + t3(:)
              DEALLOCATE( t2,t3 )
           end DO

           DO bra  = bra_min, bra_max
              bra1 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
              a    = lookup_2b_configs(3,ch1)%ival2(1,bra1)
              b    = lookup_2b_configs(3,ch1)%ival2(2,bra1)
              IF ( a == c .or. b == c ) cycle
              cind10 = clist_inv(b)
              ch_ac = pp_channel_2b%ival2(a,c)
              phase_ac = -1
              bra_ac = pp_config_2b%ival2(a,c)
              IF ( c < a ) phase_ac = -phase_ac
              
              dim1 = 1
              dim2 = ket_confs
              dim3 = number_2b(1,ch_ac)
              IF ( dim3 == 0 ) cycle
              ALLOCATE( t2(dim3) )
              t2(:) = t2_ccm(ch_ac)%cval(bra_ac,:)
              ALLOCATE( t3(dim2) )
              t3(:) = 0.d0
              CALL ZGEMM ( 'n', 'n', dim1, dim2, dim3, dcmplx(real(phase_ac),0.d0), t2, dim1, &
                   w3_temp%val1(cind10)%cval, dim3, dcmplx(1.d0,0.d0), t3, dim1 )
              t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,:) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,:) + t3(:)
              DEALLOCATE( t2,t3 )
           end DO
        end DO

        DO cind1 = 1, climit_t3(ch3)
           DEALLOCATE( w3_temp%val1(cind1)%cval )
        end DO
        DEALLOCATE( w3_temp%val1 )

     end DO
  end DO
  DEALLOCATE( clist_inv )
  
  CALL mpi_barrier(mpi_comm_world,ierror)  
  endwtime = MPI_WTIME()
  IF ( iam == 0 ) WRITE(6,'(A29,A29,f12.6,A5)') bsp29, " ...Computing T3 <- W3_3...  ", (endwtime - startwtime), ' sec.'

end SUBROUTINE t3_w3_eqn3
