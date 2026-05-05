
SUBROUTINE allocate_hbar_t2_iter
  !
  CALL allocate_hbar1b_I2
  CALL allocate_hbar1b_I3
  CALL allocate_hbar2b_I4
  CALL allocate_hbar2b_I5e_cross
  !
end SUBROUTINE allocate_hbar_t2_iter

SUBROUTINE deallocate_hbar_t2_iter
  USE constants
  !
  !
  CALL deallocate_hbar1b_I2
  CALL deallocate_hbar1b_I3
  CALL deallocate_hbar2b_I4
  CALL deallocate_hbar2b_I5e_cross
  !
end SUBROUTINE deallocate_hbar_t2_iter

SUBROUTINE build_hbar_t2_iter
  USE build_status
  USE contracts
  !
  CALL assert_built('t2',           'build_hbar_t2_iter')
  CALL assert_built('interactions', 'build_hbar_t2_iter')
  !
  CALL build_hbar1b_I2
  CALL build_hbar1b_I3
  CALL build_hbar2b_I4
  CALL build_hbar2b_I5e_cross
  !
  hbar_built = .TRUE.
  !
end SUBROUTINE build_hbar_t2_iter


SUBROUTINE allocate_hbar_paeom
  !
  CALL allocate_hbar1b_I2
  CALL allocate_hbar1b_I3
  CALL allocate_hbar2b_I3pa
  CALL allocate_hbar2b_I5pa_cross
  CALL allocate_hbar2b_I6pa
  CALL allocate_hbar3b_paeom
  !
end SUBROUTINE allocate_hbar_paeom

SUBROUTINE build_hbar_paeom
  USE constants
  !
  IF ( test0 == 4 ) then
     test = 4
     CALL build_tamp_test
  end IF
  CALL build_hbar1b_I2
  CALL build_hbar1b_I3
  CALL build_hbar2b_I3pa
  CALL build_hbar2b_I5pa_cross
  CALL build_hbar2b_I6pa
  !
end SUBROUTINE build_hbar_paeom


SUBROUTINE allocate_hbar_preom
  !
  CALL allocate_hbar1b_I2
  CALL allocate_hbar1b_I3
  CALL allocate_hbar2b_I4pr
  CALL allocate_hbar2b_I5pr_cross
  CALL allocate_hbar2b_I7pr
  CALL allocate_hbar3b_preom
  !
end SUBROUTINE allocate_hbar_preom

SUBROUTINE build_hbar_preom
  USE constants
  !
  IF ( test0 == 4 ) then
     test = 4
     CALL build_tamp_test
  end IF
  CALL build_hbar1b_I2
  CALL build_hbar1b_I3
  CALL build_hbar2b_I4pr
  CALL build_hbar2b_I5pr_cross
  CALL build_hbar2b_I7pr
  !
end SUBROUTINE build_hbar_preom



!
!     \                 \              \a
!     \a                \a             \
!     \                 \              \ ============
!     xxxxxxxxx   <--   ========   +   \/ \        /\
!     \                 \              /\ |k     c| |l
!     \b                \b           b/ \/        \/
!     \                 \            /  ------------
!
! <a|x|b> <-- <a|v|b> - (1/2).<kl|v|bc>.<ac|t|kl>
!
SUBROUTINE build_hbar1b_I2
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, bra,ket
  INTEGER :: a,b, c1,c2
  INTEGER :: dim1,dim2,dim3
  COMPLEX(dpc) :: x2
  COMPLEX(dpc), allocatable :: temp_mtx(:,:)

  IF ( iam == 0 ) write(6,*) "Computing Hbar intermediate pp (I2)..."
  
  ! <a|x|b>
  hbar1b_I2 = 0.d0
  ! <a|x|b> <-- <a|v|b>
  hbar1b_I2(:,:) = hbar1b_I2(:,:) + fock_mtx(below_ef+1:tot_orbs, below_ef+1:tot_orbs)

  ! <a|x|b> <-- -(1/2).<kl|v|bc>.<ac|t|kl> (<bc|v|kl>)
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( number_2b(3,ch) == 0 ) cycle
     ! <ac|t|kl>.<kl|v|bc>
     dim1 = number_2b(3,ch)
     dim2 = dim1
     dim3 = number_2b(1,ch)
     ALLOCATE ( temp_mtx(dim1,dim2) )
     temp_mtx = 0.d0
     CALL ZGEMM ( 'n', 't', dim1, dim2, dim3, cmplx(1.d0,0.d0,8), t2_ccm(ch)%cval, dim1, &
          conjg(v2b_pphh(ch)%cval), dim2, cmplx(1.d0,0.d0,8), temp_mtx, dim1 )
     
     DO bra = 1, dim1
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        c1  = lookup_2b_configs(3,ch)%ival2(2,bra)
        DO ket = 1, dim2
           b   = lookup_2b_configs(3,ch)%ival2(1,ket)
           c2  = lookup_2b_configs(3,ch)%ival2(2,ket)
           x2  = temp_mtx(bra,ket)
           IF ( c1 == c2 ) hbar1b_I2(a,b)   = hbar1b_I2(a,b)   - x2
           IF ( c1 == b  ) hbar1b_I2(a,c2)  = hbar1b_I2(a,c2)  + x2
           IF ( a  == c2 ) hbar1b_I2(c1,b)  = hbar1b_I2(c1,b)  + x2
           IF ( a  == b  ) hbar1b_I2(c1,c2) = hbar1b_I2(c1,c2) - x2
        end DO
     end DO
     DEALLOCATE ( temp_mtx )     
  end DO
  
  IF ( test > 1 ) CALL build_hbar1b_I2_test
  
end SUBROUTINE build_hbar1b_I2

SUBROUTINE allocate_hbar1b_I2
  USE constants
  USE operator_storage
  USE mem_tracker
  IF ( ALLOCATED(hbar1b_I2) ) return
  ALLOCATE( hbar1b_I2(below_ef+1:tot_orbs, below_ef+1:tot_orbs) )
  hbar1b_I2 = 0.d0
  CALL mem_register('hbar', REAL((tot_orbs - below_ef)**2 * 16.d0, dp))
end SUBROUTINE allocate_hbar1b_I2

SUBROUTINE deallocate_hbar1b_I2
  USE constants
  USE operator_storage
  IF ( .not. ALLOCATED(hbar1b_I2) ) return
  DEALLOCATE( hbar1b_I2 )
end SUBROUTINE deallocate_hbar1b_I2


!
!     \                 \              \j
!     \j                \j             \
!     \                 \              \ ============
!     xxxxxxxxx   <--   ========   +   \/ \        /\
!     \                 \              /\ |c     k| |d
!     \i                \i           i/ \/        \/
!     \                 \            /  ------------
!
! <i|x|j> <-- <i|v|j> + (1/2).<ik|v|cd>.<cd|t|jk>
!
SUBROUTINE build_hbar1b_I3
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, bra,ket
  INTEGER :: i,j, k1,k2
  INTEGER :: dim1,dim2,dim3
  COMPLEX(dpc) :: x2
  COMPLEX(dpc), allocatable :: temp_mtx(:,:)

  IF ( iam == 0 ) write(6,*) "Computing Hbar intermediate hh (I3)..."

  ! <i|x|j>
  hbar1b_I3 = 0.d0
  ! <i|x|j> <-- <i|v|j>
  hbar1b_I3(:,:) = hbar1b_I3(:,:) + fock_mtx(1:below_ef, 1:below_ef)

  ! <i|x|j> <-- + (1/2).<ik|v|cd>.<cd|t|jk> (<cd|v|ik>)
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( number_2b(3,ch) == 0 ) cycle
     ! <ik|v|cd>.<cd|t|jk>
     dim1 = number_2b(1,ch)
     dim2 = dim1
     dim3 = number_2b(3,ch)
     ALLOCATE ( temp_mtx(dim1,dim2) )
     temp_mtx = 0.d0
     CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), conjg(v2b_pphh(ch)%cval), dim3, &
          t2_ccm(ch)%cval, dim3, dcmplx(1.d0,0.d0), temp_mtx, dim1 )
     
     DO bra = 1, dim1
        i   = lookup_2b_configs(1,ch)%ival2(1,bra)
        k1  = lookup_2b_configs(1,ch)%ival2(2,bra)
        DO ket = 1, dim2
           j   = lookup_2b_configs(1,ch)%ival2(1,ket)
           k2  = lookup_2b_configs(1,ch)%ival2(2,ket)
           x2  = temp_mtx(bra,ket)
           IF ( k1 == k2 ) hbar1b_I3(i,j)   = hbar1b_I3(i,j)   + x2
           IF ( k1 == j  ) hbar1b_I3(i,k2)  = hbar1b_I3(i,k2)  - x2
           IF ( i  == k2 ) hbar1b_I3(k1,j)  = hbar1b_I3(k1,j)  - x2
           IF ( i  == j  ) hbar1b_I3(k1,k2) = hbar1b_I3(k1,k2) + x2
        end DO
     end DO
     DEALLOCATE ( temp_mtx )
  end DO

  IF ( test > 1 ) CALL build_hbar1b_I3_test
  
end SUBROUTINE build_hbar1b_I3

SUBROUTINE allocate_hbar1b_I3
  USE constants
  USE operator_storage
  USE mem_tracker
  IF ( ALLOCATED(hbar1b_I3) ) return
  ALLOCATE( hbar1b_I3(1:below_ef, 1:below_ef) )
  hbar1b_I3 = 0.d0
  CALL mem_register('hbar', REAL(below_ef**2 * 16.d0, dp))
end SUBROUTINE allocate_hbar1b_I3

SUBROUTINE deallocate_hbar1b_I3
  USE constants
  USE operator_storage
  IF ( .not. ALLOCATED(hbar1b_I3) ) return
  DEALLOCATE( hbar1b_I3 )
end SUBROUTINE deallocate_hbar1b_I3


!
!     \         \           \         \       \                  /
!    k\         \l         k\         \l      \                 /
!     \         \           \         \      k\                /l
!     xxxxxxxxxxx    <--    ===========   +   \  ===========  /
!     \         \           \         \       \c/\        /\d/
!    i\         \j         i\         \j      \/ \i     j/ \/
!     \         \           \         \      ----\------/----
!
!
! <ij|x|kl> <-- <ij|v|kl> + (1/2).<ij|v|cd>.<cd|t|kl>
! All Procs
SUBROUTINE build_hbar2b_I4
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, ket_confs, cd_confs
  INTEGER :: dim1,dim2,dim3
  
  IF ( iam == 0 ) write(6,*) "Computing Hbar intermediate hhhh (I4)..."

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     ket_confs = number_2b(1,ch)

     ! <ij|x|kl>
     hbar2b_I4(ch)%cval = 0.d0     
     ! <ij|x|kl> <-- <ij|v|kl>
     hbar2b_I4(ch)%cval = hbar2b_I4(ch)%cval + v2b_hhhh(ch)%cval
     
     ! <ij|x|kl> <-- + (1/2).<ij|v|cd>.<cd|t|kl>
     cd_confs = number_2b(3,ch)
     IF ( cd_confs > 0 ) THEN
        dim1 = ket_confs
        dim2 = ket_confs
        dim3 = cd_confs
        CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, cmplx(1.d0,0.d0,8), conjg(v2b_pphh(ch)%cval), dim3, &
             t2_ccm(ch)%cval, dim3, cmplx(1.d0,0.d0,8), hbar2b_I4(ch)%cval, dim1 )
     end IF     
  end DO
  
  IF ( test > 1 ) CALL build_hbar2b_I4_test
  
end SUBROUTINE build_hbar2b_I4

SUBROUTINE allocate_hbar2b_I4
  USE constants
  USE operator_storage
  USE configurations
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, bra_confs, ket_confs
  IF ( ALLOCATED(hbar2b_I4) ) return
  ALLOCATE( hbar2b_I4(channels_2b%number_confs) )
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     bra_confs = number_2b(1,ch)
     ket_confs = bra_confs
     ALLOCATE( hbar2b_I4(ch)%cval(bra_confs,ket_confs) )
     hbar2b_I4(ch)%cval = 0.d0
     CALL mem_register('hbar', REAL(bra_confs*ket_confs * 16.d0, dp))
  end DO
end SUBROUTINE allocate_hbar2b_I4

SUBROUTINE deallocate_hbar2b_I4
  USE constants
  USE operator_storage
  USE configurations
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch
  IF ( .not. ALLOCATED(hbar2b_I4) ) return
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     DEALLOCATE( hbar2b_I4(ch)%cval )
  end DO
  DEALLOCATE( hbar2b_I4 )
end SUBROUTINE deallocate_hbar2b_I4


!
!     \           \           \        \            \                  /
!    j\           \a         j\        \a           \                 /
!     \           \           \        \           j\                /a
!     x:::::x:::::x    <--    ==========   + (1/2). \  ===========  /
!     \           \           \        \            \c/\        /\k/
!    i\           \b         i\        \b           \/ \i     b/ \/
!     \           \           \        \            ---\------/---
!
! <ia|x*|jb>   <-- <ia|v|jb>   - (1/2).<ik|v|cb>.<ca|t|jk>
! <i-b|x*|j-a> <-- <i-b|v|j-a> - (1/2).<i-b|v|c-k>.<c-k|t|j-a>
! Only Procs with phhp_cross
SUBROUTINE build_hbar2b_I5e_cross
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, ck_confs
  INTEGER :: dim1,dim2,dim3
  
  IF ( iam == 0 ) write(6,*) "Computing Hbar intermediate hphp (I5* cross)..."

  DO ch = 1, channels_2bcross%number_confs
     IF ( check_my_channel_v2b_phhp_cross(ch) == 0 ) cycle
     bra_confs = number_2bcross(2,ch)
     ket_confs = bra_confs

     ! <ia|x*|jb> = <i-b|x*|j-a>
     hbar2b_I5e_cross(ch)%cval = 0.d0
     ! <ia|x*|jb> <-- <ia|v|jb>
     hbar2b_I5e_cross(ch)%cval = hbar2b_I5e_cross(ch)%cval + v2b_hphp_cross(ch)%cval

     ! <i-b|x*|j-a> <-- -(1/2).<ik|v|cb>.<ca|t|jk>
     !                  -(1/2).<i-b|v|c-k>.<c-k|t|j-a> (<cb|v|ik>.<ca|t|jk> -> <c-k|v|i-b>.<c-k|t|j-a>)
     ck_confs = number_2bcross(3,ch)
     IF ( ck_confs > 0 ) THEN
        dim1 = bra_confs
        dim2 = ket_confs
        dim3 = ck_confs
        CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, cmplx(-0.5d0,0.d0,8), conjg(v2b_phhp_cross(ch)%cval), dim3, &
             t2_ccm_cross(ch)%cval, dim3, cmplx(1.d0,0.d0,8), hbar2b_I5e_cross(ch)%cval, dim1 )
     end IF
  end DO

  IF ( test > 1 ) CALL build_hbar2b_I5e_test
  
end SUBROUTINE build_hbar2b_I5e_cross

SUBROUTINE allocate_hbar2b_I5e_cross
  USE constants
  USE operator_storage
  USE configurations
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, bra_confs, ket_confs
  IF ( ALLOCATED(hbar2b_I5e_cross) ) return
  ALLOCATE( hbar2b_I5e_cross(channels_2bcross%number_confs) )
  DO ch = 1, channels_2bcross%number_confs
     IF ( check_my_channel_v2b_phhp_cross(ch) == 0 ) cycle
     bra_confs = number_2bcross(2,ch)
     ket_confs = bra_confs
     ALLOCATE( hbar2b_I5e_cross(ch)%cval(bra_confs,ket_confs) )
     hbar2b_I5e_cross(ch)%cval = 0.d0
     CALL mem_register('hbar', REAL(bra_confs*ket_confs * 16.d0, dp))
  end DO
end SUBROUTINE allocate_hbar2b_I5e_cross

SUBROUTINE deallocate_hbar2b_I5e_cross
  USE constants
  USE operator_storage
  USE configurations
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch
  IF ( .not. ALLOCATED(hbar2b_I5e_cross) ) return
  DO ch = 1, channels_2bcross%number_confs
     IF ( check_my_channel_v2b_phhp_cross(ch) == 0 ) cycle
     DEALLOCATE( hbar2b_I5e_cross(ch)%cval )
  end DO
  DEALLOCATE( hbar2b_I5e_cross )
end SUBROUTINE deallocate_hbar2b_I5e_cross


!
!     \         \           \         \       \                  /
!    a\         \b         a\         \b      \                 /
!     \         \           \         \      a\                /b
!     xxxxxxxxxxx    <--    ===========   +   \  ===========  /
!     \         \           \         \       \k/\        /\l/
!    c\         \d         c\         \d      \/ \c     d/ \/
!     \         \           \         \       ---\------/---
!
! <ab|x|cd> <-- <ab|v|cd> + (1/2).<ab|t|kl>.<kl|v|cd>
!
SUBROUTINE build_hbar2b_I3pa
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  USE chiral_potentials
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, bra_confs, kl_confs
  INTEGER :: bra,ket,ket0, bra_min,bra_max
  INTEGER :: a,b,c,d, h
  INTEGER :: dim1,dim2,dim3, size1
  COMPLEX(dpc) :: v2b, v3b
  
  IF ( iam == 0 ) write(6,*) "Computing Hbar intermediate pppp (I3)..."

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)

     IF ( check_my_channel_paeom_pph(ch) == 0 ) cycle
     bra_min = mapping_paeom_pph(iam+1,ch,2)
     bra_max = mapping_paeom_pph(iam+1,ch,3)
     
     ! <ab|x|cd>
     hbar2b_I3(ch)%cval = 0.d0
     
     ! <ab|x|cd> <-- <ab|v|cd>
     !$omp parallel default(shared) private(bra,ket,ket0, a,b,c,d,h, v2b,v3b)
     !$omp do schedule(static)
     DO bra = bra_min, bra_max
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        b   = lookup_2b_configs(3,ch)%ival2(2,bra)
        size1 = bra + bra_confs-1
        DO ket0 = bra, size1
           ket = mod(ket0-1,bra_confs) + 1
           IF ( bra_min <= ket .and. ket < bra ) cycle
           c   = lookup_2b_configs(3,ch)%ival2(1,ket)
           d   = lookup_2b_configs(3,ch)%ival2(2,ket)
           v2b = v2int(a,b,c,d)
           v3b = 0.d0
           IF ( tnf_approx > 0 ) then
              DO h = 1, below_ef
                 v3b = v3b + v3int(a,b,h,c,d,h)
              end DO
           end IF
           hbar2b_I3(ch)%cval(bra,ket) = hbar2b_I3(ch)%cval(bra,ket) + v2b + v3b
        end DO
     end DO
     !$omp end do
     !$omp end parallel
     DO bra = bra_min, bra_max
        DO ket = bra_min, bra-1
           hbar2b_I3(ch)%cval(bra,ket) = conjg(hbar2b_I3(ch)%cval(ket,bra))
        end DO
     end DO
     
     ! <ab|x|cd> <-- + (1/2).<ab|t|kl>.<kl|v|cd>
     kl_confs = number_2b(1,ch)
     IF ( kl_confs > 0 ) THEN
        dim1 = bra_max-bra_min+1
        dim2 = bra_confs
        dim3 = kl_confs
        CALL ZGEMM ( 'n', 't', dim1, dim2, dim3, cmplx(1.d0,0.d0,8), t2_ccm(ch)%cval(bra_min:bra_max,:), dim1, &
             conjg(v2b_pphh(ch)%cval), dim2, cmplx(1.d0,0.d0,8), hbar2b_I3(ch)%cval(bra_min:bra_max,:), dim1 )
     end IF
  end DO
  
  IF ( test > 1 ) CALL build_hbar2b_I3pa_test
  
end SUBROUTINE build_hbar2b_I3pa

SUBROUTINE allocate_hbar2b_I3pa
  USE constants
  USE operator_storage
  USE configurations
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, bra_confs, bra_min,bra_max
  IF ( ALLOCATED(hbar2b_I3) ) return
  ALLOCATE( hbar2b_I3(channels_2b%number_confs) )
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     IF ( check_my_channel_paeom_pph(ch) == 0 ) cycle
     bra_min = mapping_paeom_pph(iam+1,ch,2)
     bra_max = mapping_paeom_pph(iam+1,ch,3)
     ALLOCATE( hbar2b_I3(ch)%cval(bra_min:bra_max,bra_confs) )     
     hbar2b_I3(ch)%cval = 0.d0
     CALL mem_register('hbar', REAL((bra_max-bra_min+1)*bra_confs * 16.d0, dp))
  end DO
end SUBROUTINE allocate_hbar2b_I3pa

SUBROUTINE deallocate_hbar2b_I3pa
  USE constants
  USE operator_storage
  USE configurations
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch
  IF ( .not. ALLOCATED(hbar2b_I3) ) return
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     IF ( check_my_channel_paeom_pph(ch) == 0 ) cycle
     DEALLOCATE( hbar2b_I3(ch)%cval )     
  end DO
  DEALLOCATE( hbar2b_I3 )
end SUBROUTINE deallocate_hbar2b_I3pa


!
!     \         \           \        \          \                  /
!    j\         \a         j\        \a         \                 /
!     \         \           \        \         j\                /a
!     xxxxxxxxxxx    <--    ==========    +     \  ===========  /
!     \         \           \        \          \c/\        /\k/
!    i\         \b         i\        \b         \/ \i     b/ \/
!     \         \           \        \          ---\------/---
!
! <ia|x|jb>   <-- <ia|v|jb>   - <ik|v|cb>.<ca|t|jk>
! <i-b|x|j-a> <-- <i-b|v|j-a> - <i-b|v|c-k>.<c-k|t|j-a>
!
SUBROUTINE build_hbar2b_I5pa_cross
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch2, ket_confs, ck_confs
  INTEGER :: ket_min,ket_max
  INTEGER :: dim1,dim2,dim3
  
  IF ( iam == 0 ) write(6,*) "Computing Hbar intermediate hphp (I5 cross)..."

  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_paeom_ind_cross(2*ch2) == 0 ) cycle
     ket_confs = number_2bcross(2,ch2)
     
     IF ( check_my_channel_paeom_php_cross(ch2) == 0 ) cycle
     ket_min = mapping_paeom_php_cross(iam+1,ch2,2)
     ket_max = mapping_paeom_php_cross(iam+1,ch2,3)
     
     ! <ia|x|jb> = <i-b|x|j-a>
     hbar2b_I5_cross(ch2)%cval = 0.d0     
     ! <i-b|x|j-a> <-- <i-b|v|j-a>
     hbar2b_I5_cross(ch2)%cval(:,ket_min:ket_max) = hbar2b_I5_cross(ch2)%cval(:,ket_min:ket_max) &
          + v2b_hphp_cross(ch2)%cval(:,ket_min:ket_max)
     
     ! <i-b|x|j-a> <-- -<ik|v|cb>.<ca|t|jk>
     !                 -<i-b|v|c-k>.<c-k|t|j-a> (<cb|v|ik>.<ca|t|jk> -> <c-k|v|i-b>.<c-k|t|j-a>)
     ck_confs = number_2bcross(3,ch2)
     IF ( ck_confs > 0 ) THEN
        dim1 = ket_confs
        dim2 = ket_max-ket_min+1
        dim3 = ck_confs
        CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, cmplx(-1.d0,0.d0,8), conjg(v2b_phhp_cross(ch2)%cval), &
             dim3, t2_ccm_cross(ch2)%cval(:,ket_min:ket_max), dim3, cmplx(1.d0,0.d0,8), &
             hbar2b_I5_cross(ch2)%cval(:,ket_min:ket_max), dim1 )
     end IF
  end DO
  
  IF ( test > 1 ) CALL build_hbar2b_I5pa_test
  
end SUBROUTINE build_hbar2b_I5pa_cross

SUBROUTINE allocate_hbar2b_I5pa_cross
  USE constants
  USE operator_storage
  USE configurations
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, ket_confs, ket_min,ket_max
  IF ( ALLOCATED(hbar2b_I5_cross) ) return
  ALLOCATE( hbar2b_I5_cross(channels_2bcross%number_confs) )
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch) == 0 ) cycle
     IF ( r2_paeom_ind_cross(2*ch) == 0 ) cycle
     ket_confs = number_2bcross(2,ch)
     IF ( check_my_channel_paeom_php_cross(ch) == 0 ) cycle
     ket_min = mapping_paeom_php_cross(iam+1,ch,2)
     ket_max = mapping_paeom_php_cross(iam+1,ch,3)
     ALLOCATE( hbar2b_I5_cross(ch)%cval(ket_confs,ket_min:ket_max) )
     hbar2b_I5_cross(ch)%cval = 0.d0
     CALL mem_register('hbar', REAL(ket_confs*(ket_max-ket_min+1) * 16.d0, dp))
  end DO
end SUBROUTINE allocate_hbar2b_I5pa_cross

SUBROUTINE deallocate_hbar2b_I5pa_cross
  USE constants
  USE operator_storage
  USE configurations
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch
  IF ( .not. ALLOCATED(hbar2b_I5_cross) ) return
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch) == 0 ) cycle
     IF ( r2_paeom_ind_cross(2*ch) == 0 ) cycle
     IF ( check_my_channel_paeom_php_cross(ch) == 0 ) cycle
     DEALLOCATE( hbar2b_I5_cross(ch)%cval )
  end DO
  DEALLOCATE( hbar2b_I5_cross )
end SUBROUTINE deallocate_hbar2b_I5pa_cross


!
!     \       \  /       \       \  /      \         \     /          /  \     /
!    b\      a\ /i      b\      a\ /i     b\        i\    /a        b/   \i   /a
!     \       \/         \       \/        =======   \   /          /    \   /
!     xxxxxxxxx    <--   =========    +    \    /\   \  /   +   ==========  /
!     \                  \                 \  d| |k  \ /       /\ /     k\ /
!    c\                 c\                c\   \/    \/      c/l\/       \/
!     \                  \                 \  ---------      /  -----------
!
! <ab|x|ic> <-- <ab|v|ic> + P(ab).<ad|t|ik>.<kb|v|dc> + (1/2).<ab|t|kl>.<kl|v|ic>
!
SUBROUTINE build_hbar2b_I6pa
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, ch0,ch2, bra_confs,ket_confs, bra_min,bra_max
  INTEGER :: bra,ket, bra0,ket0, bra2, kd, kl_confs,kd_confs
  INTEGER :: i,a,b,c, k,d
  INTEGER :: dim1,dim2,dim3, phase
  COMPLEX(dpc) :: v2, t2
  
  IF ( iam == 0 ) write(6,*) "Computing Hbar intermediate pphp (I6)..."

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(2,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(2,ch)
     
     IF ( check_my_channel_paeom_pph(ch) == 0 ) cycle
     bra_min = mapping_paeom_pph(iam+1,ch,2)
     bra_max = mapping_paeom_pph(iam+1,ch,3)

     ! <ab|x|ic>
     hbar2b_I6(ch)%cval = 0.d0     
     ! <ab|x|ic> <-- <ab|v|ic>
     hbar2b_I6(ch)%cval(bra_min:bra_max,:) = hbar2b_I6(ch)%cval(bra_min:bra_max,:) + v2b_pphp(ch)%cval(bra_min:bra_max,:)
     
     ! <ab|x|ic> <-- +(1/2).<ab|t|kl>.<kl|v|ic>
     kl_confs = number_2b(1,ch)
     IF ( kl_confs > 0 ) THEN
        dim1 = bra_max - bra_min + 1
        dim2 = ket_confs
        dim3 = kl_confs
        CALL ZGEMM ( 'n', 't', dim1, dim2, dim3, cmplx(1.d0,0.d0,8), t2_ccm(ch)%cval(bra_min:bra_max,:), dim1, &
             conjg(v2b_hphh(ch)%cval), dim2, cmplx(1.d0,0.d0,8), hbar2b_I6(ch)%cval(bra_min:bra_max,:), dim1 )
     end IF

     ! <ab|x|ic> <-- +P(ab).<ad|t|ik>.<kb|v|dc> = +P(ab).<ad|t|ki>.<kb|v|cd>
     ! <ab|x|ic> = <a-i|x|c-b> <-- P(ab).<a-i|t|k-d>.<k-d|v|c-b>
     DO bra = bra_min, bra_max
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        b   = lookup_2b_configs(3,ch)%ival2(2,bra)
        DO ket = 1, ket_confs
           i   = lookup_2b_configs(2,ch)%ival2(1,ket)
           c   = lookup_2b_configs(2,ch)%ival2(2,ket)

           ch2 = ph_channel_2bcross%ival2(a,i)
           IF ( ch2 == 0 ) cycle
           bra2 = ph_config_2bcross%ival2(a,i)
           kd_confs = number_2bcross(2,ch2)
           DO kd = 1, kd_confs
              k  = lookup_2bcross_configs(2,ch2)%ival2(1,kd)
              d  = lookup_2bcross_configs(2,ch2)%ival2(2,kd)
              t2 = t2_ccm_cross(ch2)%cval(bra2,kd)
              
              ch0 = hp_channel_2b%ival2(k,b)
              IF ( ch0 == 0 ) cycle
              IF ( ch0 /= pp_channel_2b%ival2(c,d) ) cycle
              
              phase = 1
              bra0 = hp_config_2b%ival2(k,b)
              ket0 = pp_config_2b%ival2(c,d)
              IF ( d < c ) phase = -phase
              v2 = phase * conjg(v2b_pphp(ch0)%cval(ket0,bra0))
              hbar2b_I6(ch)%cval(bra,ket) = hbar2b_I6(ch)%cval(bra,ket) + t2*v2
           end DO
        end DO
        ! P(ab)
        DO ket = 1, ket_confs
           i   = lookup_2b_configs(2,ch)%ival2(1,ket)
           c   = lookup_2b_configs(2,ch)%ival2(2,ket)
           
           ch2 = ph_channel_2bcross%ival2(b,i)
           IF ( ch2 == 0 ) cycle
           bra2 = ph_config_2bcross%ival2(b,i)
           kd_confs = number_2bcross(2,ch2)
           DO kd = 1, kd_confs
              k  = lookup_2bcross_configs(2,ch2)%ival2(1,kd)
              d  = lookup_2bcross_configs(2,ch2)%ival2(2,kd)
              t2 = t2_ccm_cross(ch2)%cval(bra2,kd)
              
              ch0 = hp_channel_2b%ival2(k,a)
              IF ( ch0 == 0 ) cycle
              IF ( ch0 /= pp_channel_2b%ival2(c,d) ) cycle
              
              phase = 1
              bra0 = hp_config_2b%ival2(k,a)
              ket0 = pp_config_2b%ival2(c,d)
              IF ( d < c ) phase = -phase
              v2 = phase * conjg(v2b_pphp(ch0)%cval(ket0,bra0))
              hbar2b_I6(ch)%cval(bra,ket) = hbar2b_I6(ch)%cval(bra,ket) - t2*v2
           end DO
        end DO
     end DO
  end DO

  IF ( test > 1 ) CALL build_hbar2b_I6pa_test
  
end SUBROUTINE build_hbar2b_I6pa

SUBROUTINE allocate_hbar2b_I6pa
  USE constants
  USE operator_storage
  USE configurations
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, ket_confs, bra_min,bra_max
  IF ( ALLOCATED(hbar2b_I6) ) return
  ALLOCATE( hbar2b_I6(channels_2bcross%number_confs) )
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(2,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     ket_confs = number_2b(2,ch)     
     IF ( check_my_channel_paeom_pph(ch) == 0 ) cycle
     bra_min = mapping_paeom_pph(iam+1,ch,2)
     bra_max = mapping_paeom_pph(iam+1,ch,3)
     ALLOCATE( hbar2b_I6(ch)%cval(bra_min:bra_max,ket_confs) )
     hbar2b_I6(ch)%cval = 0.d0
     CALL mem_register('hbar', REAL((bra_max-bra_min+1)*ket_confs * 16.d0, dp))
  end DO
end SUBROUTINE allocate_hbar2b_I6pa

SUBROUTINE deallocate_hbar2b_I6pa
  USE constants
  USE operator_storage
  USE configurations
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch
  IF ( .not. ALLOCATED(hbar2b_I6) ) return
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(2,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     IF ( check_my_channel_paeom_pph(ch) == 0 ) cycle
     DEALLOCATE( hbar2b_I6(ch)%cval )
  end DO
  DEALLOCATE( hbar2b_I6 )
end SUBROUTINE deallocate_hbar2b_I6pa


!
!                         ============
!     xxxxxxxxx   <--    / \        /\
!     \                 /  |c     k| |d
!     \i              i/  /        \/
!     \               /  ############
!
! <i|z3b|> <-- (1/2).<ik|v|cd>.<cd|r|k>
!
SUBROUTINE build_hbar3b_paeom
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, bra,ket, bra_min,bra_max
  INTEGER :: k,k1,k2, dim1,dim2,dim3
  COMPLEX(dpc) :: x2
  COMPLEX(dpc), allocatable :: temp_mtx(:,:)

  IF ( iam == 0 ) write(6,*) "Computing PA-EOM Hbar 3b intermediate (I3)..."
  
  ! <i|z3b|>
  hbar3b_paeom_I3 = 0.d0
  ! <i|z3b|> <-- + (1/2).<ik|v|cd>.<cd|r|k> (<cd|v|ik>)
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
  
     IF ( check_my_channel_paeom_pph(ch) == 0 ) cycle
     bra_min = mapping_paeom_pph(iam+1,ch,2)
     bra_max = mapping_paeom_pph(iam+1,ch,3)
     
     ! <ik|v|cd>.<cd|r|k>
     dim1 = number_2b(1,ch)
     dim2 = 2
     dim3 = bra_max-bra_min + 1
     ALLOCATE ( temp_mtx(dim1,2) )
     temp_mtx = 0.d0
     CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, cmplx(1.d0,0.d0,8), conjg(v2b_pphh(ch)%cval(bra_min:bra_max,:)), dim3, &
          r2_paeom(ch)%cval(bra_min:bra_max,:), dim3, cmplx(1.d0,0.d0,8), temp_mtx, dim1 )
     
     DO bra = 1, dim1
        k1  = lookup_2b_configs(1,ch)%ival2(1,bra)
        k2  = lookup_2b_configs(1,ch)%ival2(2,bra)
        DO ket = 1, 2
           k   = r2_paeom_ind(2*(ch-1)+ket)
           x2  = temp_mtx(bra,ket)
           IF ( k2 == k ) hbar3b_paeom_I3(k1) = hbar3b_paeom_I3(k1) + x2
           IF ( k1 == k ) hbar3b_paeom_I3(k2) = hbar3b_paeom_I3(k2) - x2
        end DO
     end DO
     DEALLOCATE ( temp_mtx )
  end DO
  CALL mpi_allreduce(mpi_in_place,hbar3b_paeom_I3,size(hbar3b_paeom_I3),mpi_COMPLEX16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'build_hbar3b_paeom', 'allreduce')
  
  IF ( test > 1 ) CALL build_hbar3b_paeom_test
  
end SUBROUTINE build_hbar3b_paeom

SUBROUTINE allocate_hbar3b_paeom
  USE constants
  USE operator_storage
  USE mem_tracker
  IF ( ALLOCATED(hbar3b_paeom_I3) ) return
  ALLOCATE( hbar3b_paeom_I3(1:below_ef) )
  hbar3b_paeom_I3 = 0.d0
  CALL mem_register('hbar', REAL(below_ef * 16.d0, dp))
end SUBROUTINE allocate_hbar3b_paeom

SUBROUTINE deallocate_hbar3b_paeom
  USE constants
  USE operator_storage
  IF ( .not. ALLOCATED(hbar3b_paeom_I3) ) return
  DEALLOCATE( hbar3b_paeom_I3 )
end SUBROUTINE deallocate_hbar3b_paeom



!
!     \         \           \         \       \                  /
!    k\         \l         k\         \l      \                 /
!     \         \           \         \      k\                /l
!     xxxxxxxxxxx    <--    ===========   +   \  ===========  /
!     \         \           \         \       \c/\        /\d/
!    i\         \j         i\         \j      \/ \i     j/ \/
!     \         \           \         \      ----\------/----
!
!
! <ij|x|kl> <-- <ij|v|kl> + (1/2).<ij|v|cd>.<cd|t|kl>
!
SUBROUTINE build_hbar2b_I4pr
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, bra_confs, cd_confs
  INTEGER :: ket_min,ket_max
  INTEGER :: dim1,dim2,dim3
  
  IF ( iam == 0 ) write(6,*) "Computing Hbar intermediate hhhh (I4)..."

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     bra_confs = number_2b(1,ch)

     IF ( check_my_channel_preom_phh(ch) == 0 ) cycle
     ket_min = mapping_preom_phh(iam+1,ch,2)
     ket_max = mapping_preom_phh(iam+1,ch,3)
     
     ! <ij|x|kl>
     hbar2b_I4(ch)%cval = 0.d0
     ! <ij|x|kl> <-- <ij|v|kl>
     hbar2b_I4(ch)%cval(:,ket_min:ket_max) = hbar2b_I4(ch)%cval(:,ket_min:ket_max) + v2b_hhhh(ch)%cval(:,ket_min:ket_max)
     
     ! <ij|x|kl> <-- + (1/2).<ij|v|cd>.<cd|t|kl>
     cd_confs = number_2b(3,ch)
     IF ( cd_confs > 0 ) THEN
        dim1 = bra_confs
        dim2 = ket_max - ket_min + 1
        dim3 = cd_confs
        CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), conjg(v2b_pphh(ch)%cval), dim3, &
             t2_ccm(ch)%cval(:,ket_min:ket_max), dim3, dcmplx(1.d0,0.d0), hbar2b_I4(ch)%cval(:,ket_min:ket_max), dim1 )
     end IF
  end DO
  
  IF ( test > 1 ) CALL build_hbar2b_I4_test
  
end SUBROUTINE build_hbar2b_I4pr

SUBROUTINE allocate_hbar2b_I4pr
  USE constants
  USE operator_storage
  USE configurations
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, bra_confs, ket_min,ket_max
  IF ( ALLOCATED(hbar2b_I4) ) return
  ALLOCATE( hbar2b_I4(channels_2b%number_confs) )
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     bra_confs = number_2b(1,ch)
     IF ( check_my_channel_preom_phh(ch) == 0 ) cycle
     ket_min = mapping_preom_phh(iam+1,ch,2)
     ket_max = mapping_preom_phh(iam+1,ch,3)
     ALLOCATE( hbar2b_I4(ch)%cval(bra_confs,ket_min:ket_max) )
     hbar2b_I4(ch)%cval = 0.d0
     CALL mem_register('hbar', REAL(bra_confs*(ket_max-ket_min+1) * 16.d0, dp))
  end DO
end SUBROUTINE allocate_hbar2b_I4pr

SUBROUTINE deallocate_hbar2b_I4pr
  USE constants
  USE operator_storage
  USE configurations
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch
  IF ( .not. ALLOCATED(hbar2b_I4) ) return
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     IF ( check_my_channel_preom_phh(ch) == 0 ) cycle
     DEALLOCATE( hbar2b_I4(ch)%cval )
  end DO
  DEALLOCATE( hbar2b_I4 )
end SUBROUTINE deallocate_hbar2b_I4pr


!
!     \         \           \        \          \                  /
!    j\         \a         j\        \a         \                 /
!     \         \           \        \         j\                /a
!     xxxxxxxxxxx    <--    ==========    +     \  ===========  /
!     \         \           \        \          \c/\        /\k/
!    i\         \b         i\        \b         \/ \i     b/ \/
!     \         \           \        \          ---\------/---
!
! <ia|x|jb>   <-- <ia|v|jb>   - <ik|v|cb>.<ca|t|jk>
! <i-b|x|j-a> <-- <i-b|v|j-a> - <i-b|v|c-k>.<c-k|t|j-a>
!
SUBROUTINE build_hbar2b_I5pr_cross
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch2, ket_confs, ck_confs
  INTEGER :: ket_min,ket_max
  INTEGER :: dim1,dim2,dim3
  
  IF ( iam == 0 ) write(6,*) "Computing Hbar intermediate hphp (I5 cross)..."

  DO ch2 = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch2) == 0 ) cycle
     IF ( r2_preom_ind_cross(2*ch2) == 0 ) cycle
     ket_confs = number_2bcross(2,ch2)
     
     IF ( check_my_channel_preom_hhp_cross(ch2) == 0 ) cycle
     ket_min = mapping_preom_hhp_cross(iam+1,ch2,2)
     ket_max = mapping_preom_hhp_cross(iam+1,ch2,3)
     
     ! <ia|x|jb> = <i-b|x|j-a>
     hbar2b_I5_cross(ch2)%cval = 0.d0     
     ! <i-b|x|j-a> <-- <i-b|v|j-a>
     hbar2b_I5_cross(ch2)%cval(:,ket_min:ket_max) = hbar2b_I5_cross(ch2)%cval(:,ket_min:ket_max) &
          + v2b_hphp_cross(ch2)%cval(:,ket_min:ket_max)
     
     ! <i-b|x|j-a> <-- -<ik|v|cb>.<ca|t|jk>
     !                 -<i-b|v|c-k>.<c-k|t|j-a> (<cb|v|ik>.<ca|t|jk> -> <c-k|v|i-b>.<c-k|t|j-a>)
     ck_confs = number_2bcross(3,ch2)
     IF ( ck_confs > 0 ) THEN
        dim1 = ket_confs
        dim2 = ket_max-ket_min+1
        dim3 = ck_confs
        CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, dcmplx(-1.d0,0.d0), conjg(v2b_phhp_cross(ch2)%cval), &
             dim3, t2_ccm_cross(ch2)%cval(:,ket_min:ket_max), dim3, dcmplx(1.d0,0.d0), &
             hbar2b_I5_cross(ch2)%cval(:,ket_min:ket_max), dim1 )
     end IF
  end DO
  
  IF ( test > 1 ) CALL build_hbar2b_I5pr_test
  
end SUBROUTINE build_hbar2b_I5pr_cross

SUBROUTINE allocate_hbar2b_I5pr_cross
  USE constants
  USE operator_storage
  USE configurations
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, ket_confs, ket_min,ket_max
  IF ( ALLOCATED(hbar2b_I5_cross) ) return
  ALLOCATE( hbar2b_I5_cross(channels_2bcross%number_confs) )
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch) == 0 ) cycle
     IF ( r2_preom_ind_cross(2*ch) == 0 ) cycle
     ket_confs = number_2bcross(2,ch)
     IF ( check_my_channel_preom_hhp_cross(ch) == 0 ) cycle
     ket_min = mapping_preom_hhp_cross(iam+1,ch,2)
     ket_max = mapping_preom_hhp_cross(iam+1,ch,3)
     ALLOCATE( hbar2b_I5_cross(ch)%cval(ket_confs,ket_min:ket_max) )
     hbar2b_I5_cross(ch)%cval = 0.d0
     CALL mem_register('hbar', REAL(ket_confs*(ket_max-ket_min+1) * 16.d0, dp))
  end DO
end SUBROUTINE allocate_hbar2b_I5pr_cross

SUBROUTINE deallocate_hbar2b_I5pr_cross
  USE constants
  USE operator_storage
  USE configurations
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch
  IF ( .not. ALLOCATED(hbar2b_I5_cross) ) return
  DO ch = 1, channels_2bcross%number_confs
     IF ( number_2bcross(2,ch) == 0 ) cycle
     IF ( r2_preom_ind_cross(2*ch) == 0 ) cycle
     IF ( check_my_channel_preom_hhp_cross(ch) == 0 ) cycle
     DEALLOCATE( hbar2b_I5_cross(ch)%cval )
  end DO
  DEALLOCATE( hbar2b_I5_cross )
end SUBROUTINE deallocate_hbar2b_I5pr_cross


!
!     \        \  /       \       \  /     \         \     /          /  \     /
!    j\       k\ /a      j\      k\ /a    j\        k\    /a        j/   \a   /k
!     \        \/         \       \/       xxxxxxx   \   /          /    \   /
!     xxxxxxxxxx    <--   =========    +   \    /\   \  /   +   xxxxxxxxxx  /
!     \                   \                \  l| |c  \ /       /\ /     d\ /
!    i\                  i\               i\   \/    \/      i/c\/       \/
!     \                   \                \  ---------      /  -----------
!
! <ia|x|jk> <-- <ia|v|jk> + P(jk).<ca|t|lk>.<il|v|jc> + (1/2).<ia|v|cd>.<cd|t|jk>
!
SUBROUTINE build_hbar2b_I7pr
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, ch0,ch2, bra_confs,ket_confs, ket_min,ket_max
  INTEGER :: bra,ket, bra0,ket0, ket2, cl, cd_confs,cl_confs
  INTEGER :: i,a,j,k, c,l
  INTEGER :: dim1,dim2,dim3, phase
  COMPLEX(dpc) :: v2, t2
  
  IF ( iam == 0 ) write(6,*) "Computing Hbar intermediate hphh (I7)..."

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(2,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     bra_confs = number_2b(2,ch)
     ket_confs = number_2b(1,ch)
     
     IF ( check_my_channel_preom_phh(ch) == 0 ) cycle
     ket_min = mapping_preom_phh(iam+1,ch,2)
     ket_max = mapping_preom_phh(iam+1,ch,3)

     ! <ia|x|jk>
     hbar2b_I7(ch)%cval = 0.d0
     ! <ia|x|jk> <-- <ia|v|jk>
     hbar2b_I7(ch)%cval(:,ket_min:ket_max) = hbar2b_I7(ch)%cval(:,ket_min:ket_max) + v2b_hphh(ch)%cval(:,ket_min:ket_max)
     
     ! <ia|x|jk> <-- +(1/2).<ia|v|cd>.<cd|t|jk>
     cd_confs = number_2b(3,ch)
     IF ( cd_confs > 0 ) THEN
        dim1 = bra_confs
        dim2 = ket_max - ket_min + 1
        dim3 = cd_confs
        CALL ZGEMM ( 't', 'n', dim1, dim2, dim3, dcmplx(1.d0,0.d0), conjg(v2b_pphp(ch)%cval), dim3, &
             t2_ccm(ch)%cval(:,ket_min:ket_max), dim3, dcmplx(1.d0,0.d0), hbar2b_I7(ch)%cval(:,ket_min:ket_max), dim1 )
     end IF

     ! <ia|x|jk> <-- +P(jk).<ca|t|lk>.<il|x|jc> = -P(jk).<ca|t|kl>.<il|v|jc>
     ! <ia|x|jk> = <i-j|x|k-a> <-- -P(jk).<c-l|t|k-a>.<i-j|v|c-l>
     DO ket = ket_min, ket_max
        j   = lookup_2b_configs(1,ch)%ival2(1,ket)
        k   = lookup_2b_configs(1,ch)%ival2(2,ket)
        DO bra = 1, bra_confs
           i   = lookup_2b_configs(2,ch)%ival2(1,bra)
           a   = lookup_2b_configs(2,ch)%ival2(2,bra)
           
           ch2 = hp_channel_2bcross%ival2(k,a)
           IF ( ch2 == 0 ) cycle
           ket2 = hp_config_2bcross%ival2(k,a)
           cl_confs = number_2bcross(3,ch2)
           DO cl = 1, cl_confs
              c  = lookup_2bcross_configs(3,ch2)%ival2(1,cl)
              l  = lookup_2bcross_configs(3,ch2)%ival2(2,cl)
              t2 = t2_ccm_cross(ch2)%cval(cl,ket2)

              ch0 = hp_channel_2b%ival2(j,c)
              IF ( ch0 == 0 ) cycle
              IF ( ch0 /= hh_channel_2b%ival2(i,l) ) cycle

              phase = 1
              ket0 = hp_config_2b%ival2(j,c)
              bra0 = hh_config_2b%ival2(i,l)
              IF ( l < i ) phase = -phase
              v2 = phase * conjg(v2b_hphh(ch0)%cval(ket0,bra0))
              hbar2b_I7(ch)%cval(bra,ket) = hbar2b_I7(ch)%cval(bra,ket) - t2*v2
           end DO
        end DO
        ! P(jk)
        DO bra = 1, bra_confs
           i   = lookup_2b_configs(2,ch)%ival2(1,bra)
           a   = lookup_2b_configs(2,ch)%ival2(2,bra)

           ch2 = hp_channel_2bcross%ival2(j,a)
           IF ( ch2 == 0 ) cycle
           ket2 = hp_config_2bcross%ival2(j,a)
           cl_confs = number_2bcross(3,ch2)
           DO cl = 1, cl_confs
              c  = lookup_2bcross_configs(3,ch2)%ival2(1,cl)
              l  = lookup_2bcross_configs(3,ch2)%ival2(2,cl)
              t2 = t2_ccm_cross(ch2)%cval(cl,ket2)

              ch0 = hp_channel_2b%ival2(k,c)
              IF ( ch0 == 0 ) cycle
              IF ( ch0 /= hh_channel_2b%ival2(i,l) ) cycle

              phase = 1
              ket0 = hp_config_2b%ival2(k,c)
              bra0 = hh_config_2b%ival2(i,l)
              IF ( l < i ) phase = -phase
              v2 = phase * conjg(v2b_hphh(ch0)%cval(ket0,bra0))
              hbar2b_I7(ch)%cval(bra,ket) = hbar2b_I7(ch)%cval(bra,ket) + t2*v2
           end DO
        end DO
     end DO
  end DO

  IF ( test > 1 ) CALL build_hbar2b_I7pr_test
  
end SUBROUTINE build_hbar2b_I7pr

SUBROUTINE allocate_hbar2b_I7pr
  USE constants
  USE operator_storage
  USE configurations
  USE mem_tracker
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, bra_confs, ket_min,ket_max
  IF ( ALLOCATED(hbar2b_I7) ) return
  ALLOCATE( hbar2b_I7(channels_2bcross%number_confs) )
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(2,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     bra_confs = number_2b(2,ch)
     IF ( check_my_channel_preom_phh(ch) == 0 ) cycle
     ket_min = mapping_preom_phh(iam+1,ch,2)
     ket_max = mapping_preom_phh(iam+1,ch,3)
     ALLOCATE( hbar2b_I7(ch)%cval(bra_confs,ket_min:ket_max) )
     hbar2b_I7(ch)%cval = 0.d0
     CALL mem_register('hbar', REAL(bra_confs*(ket_max-ket_min+1) * 16.d0, dp))
  end DO
end SUBROUTINE allocate_hbar2b_I7pr

SUBROUTINE deallocate_hbar2b_I7pr
  USE constants
  USE operator_storage
  USE configurations
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch
  IF ( .not. ALLOCATED(hbar2b_I7) ) return
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(2,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     IF ( check_my_channel_preom_phh(ch) == 0 ) cycle
     DEALLOCATE( hbar2b_I7(ch)%cval )
  end DO
  DEALLOCATE( hbar2b_I7 )
end SUBROUTINE deallocate_hbar2b_I7pr


!
!                         ============
!     xxxxxxxxx   <--    / \        /\
!     \                 /  |k     c| |l
!     \a              a/  /        \/
!     \               /  ############
!
! <|z3b|a> <-- -(1/2).<kl|v|ac>.<c|r|kl>
!
SUBROUTINE build_hbar3b_preom
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  USE mpi_check
  
  IMPLICIT NONE
  INTEGER :: ch, bra,ket, ket_min,ket_max
  INTEGER :: c,c1,c2, dim1,dim2,dim3
  COMPLEX(dpc) :: x2
  COMPLEX(dpc), allocatable :: temp_mtx(:,:)

  IF ( iam == 0 ) write(6,*) "Computing PR-EOM Hbar 3b intermediate (I2)..."
  
  ! <|z3b|a>
  hbar3b_preom_I2 = 0.d0
  ! <|z3b|a> <-- -(1/2).<kl|v|ac>.<c|r|kl>
  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
  
     IF ( check_my_channel_preom_phh(ch) == 0 ) cycle
     ket_min = mapping_preom_phh(iam+1,ch,2)
     ket_max = mapping_preom_phh(iam+1,ch,3)
     
     ! <c|r|kl>.<kl|v|ac>
     dim1 = 2
     dim2 = number_2b(3,ch)
     dim3 = ket_max-ket_min + 1
     ALLOCATE ( temp_mtx(2,dim2) )
     temp_mtx = 0.d0
     CALL ZGEMM ( 'n', 't', dim1, dim2, dim3, dcmplx(1.d0,0.d0), r2_preom(ch)%cval(:,ket_min:ket_max), dim1, &
          conjg(v2b_pphh(ch)%cval(:,ket_min:ket_max)), dim2, dcmplx(1.d0,0.d0), temp_mtx, dim1 )

     DO ket = 1, dim2
        c1  = lookup_2b_configs(3,ch)%ival2(1,ket)
        c2  = lookup_2b_configs(3,ch)%ival2(2,ket)
        DO bra = 1, 2
           c   = r2_preom_ind(2*(ch-1)+bra)
           x2  = temp_mtx(bra,ket)
           IF ( c2 == c ) hbar3b_preom_I2(c1) = hbar3b_preom_I2(c1) - x2
           IF ( c1 == c ) hbar3b_preom_I2(c2) = hbar3b_preom_I2(c2) + x2
        end DO
     end DO     
     DEALLOCATE ( temp_mtx )
  end DO
  CALL mpi_allreduce(mpi_in_place,hbar3b_preom_I2,size(hbar3b_preom_I2),mpi_complex16,mpi_sum,mpi_comm_world,ierror)
  CALL check_mpi(ierror, 'build_hbar3b_preom', 'allreduce')
  
  IF ( test > 1 ) CALL build_hbar3b_preom_test
  
end SUBROUTINE build_hbar3b_preom

SUBROUTINE allocate_hbar3b_preom
  USE constants
  USE operator_storage
  USE mem_tracker
  IF ( ALLOCATED(hbar3b_preom_I2) ) return
  ALLOCATE( hbar3b_preom_I2(below_ef+1:tot_orbs) )
  hbar3b_preom_I2 = 0.d0
  CALL mem_register('hbar', REAL((tot_orbs-below_ef) * 16.d0, dp))
end SUBROUTINE allocate_hbar3b_preom

SUBROUTINE deallocate_hbar3b_preom
  USE constants
  USE operator_storage
  IF ( .not. ALLOCATED(hbar3b_preom_I2) ) return
  DEALLOCATE( hbar3b_preom_I2 )
end SUBROUTINE deallocate_hbar3b_preom
