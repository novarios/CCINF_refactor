
SUBROUTINE setup_tests
  USE parallel
  USE single_particle_orbits
  USE constants  
  USE operator_storage
  USE configurations
  USE chiral_potentials

  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, bra_min,bra_max, bra,ket
  INTEGER :: p,q,r,s, h,h1,h2,h3, a,b,c,d,i,j,k,l
  COMPLEX(dpc) :: v2b, v3b, e0_test
  
  IF ( iam == 0 ) write(6,*)
  IF ( iam == 0 ) write(6,*) '...Setting up Test Interaction...'

  ! Vacuum Energy
  e0_test = 0.d0
  
  DO h = 1, below_ef
     e0_test = e0_test + all_orbit%e(h)
  end DO

  DO h1 = 1, below_ef-1
     DO h2 = h1+1, below_ef
        v2b = v2int(h1,h2,h1,h2)
        e0_test = e0_test + v2b
     end DO
  end DO

  IF ( tnf_approx > 0 ) then
     DO h1 = 1, below_ef-2
        DO h2 = h1+1, below_ef-1
           DO h3 = h2+1, below_ef
              v3b = v3int(h1,h2,h3,h1,h2,h3)
              e0_test = e0_test + v3b
           end DO
        end DO
     end DO
  end IF
  
  IF ( test == 1 .and. iam == 0 ) then
     IF ( abs( e0_test - e0 ) > 1.e-5 ) then
        WRITE(6,*) "vacuum_energy: ", e0_test, e0
        stop
     end IF
  end IF
  
  
  ! Fock Matrix
  ALLOCATE( fock_test(1:tot_orbs, 1:tot_orbs) )
  fock_test = 0.d0
  
  DO p = 1, tot_orbs
     fock_test(p,p) = fock_test(p,p) + all_orbit%e(p)
  end DO

  DO p = 1, tot_orbs
     DO q = 1, tot_orbs
        DO h = 1, below_ef
           v2b = v2int(p,h,q,h)
           fock_test(p,q) = fock_test(p,q) + v2b
           ! IF ( abs(v2b) > 1.e-6 ) write(6,*) 'fock_test: ', p,q,h, v2b
        end DO
     end DO
  end DO

  IF ( tnf_approx > 0 ) then
     DO p = 1, tot_orbs
        DO q = 1, tot_orbs
           DO h1 = 1, below_ef-1
              DO h2 = h1+1, below_ef
                 v3b = v3int(p,h1,h2,q,h1,h2)
                 fock_test(p,q) = fock_test(p,q) + v3b
              end DO
           end DO
        end DO
     end DO
  end IF

  IF ( test == 1 .and. iam == 0 ) then
     DO p = 1, tot_orbs
        DO q = 1, tot_orbs
           IF ( abs( fock_test(p,q) - fock_mtx(p,q) ) > 1.e-5 ) then
              WRITE(6,*) "fock_mtx: ", p, q, fock_test(p,q), fock_mtx(p,q)
              stop
           end IF
        end DO
     end DO
  end IF


  ! 2N Interaction
  ALLOCATE( v2b_test(1:tot_orbs, 1:tot_orbs, 1:tot_orbs, 1:tot_orbs) )
  v2b_test = 0.d0
  
  DO p = 1, tot_orbs-1
     DO q = p+1, tot_orbs
        DO r = 1, tot_orbs-1
           DO s = r+1, tot_orbs
              v2b = v2int(p,q,r,s)
              v2b_test(p,q,r,s) = v2b_test(p,q,r,s) + v2b
              v2b_test(q,p,r,s) = v2b_test(q,p,r,s) - v2b
              v2b_test(p,q,s,r) = v2b_test(p,q,s,r) - v2b
              v2b_test(q,p,s,r) = v2b_test(q,p,s,r) + v2b
              IF ( tnf_approx > 0 ) then
                 DO h = 1, below_ef
                    IF ( h == p .or. h == q ) cycle
                    IF ( h == r .or. h == s ) cycle
                    v3b = v3int(p,q,h,r,s,h)
                    v2b_test(p,q,r,s) = v2b_test(p,q,r,s) + v3b
                    v2b_test(q,p,r,s) = v2b_test(q,p,r,s) - v3b
                    v2b_test(p,q,s,r) = v2b_test(p,q,s,r) - v3b
                    v2b_test(q,p,s,r) = v2b_test(q,p,s,r) + v3b
                 end DO
              end IF
           end DO
        end DO
     end DO
  end DO

  IF ( test == 1 .and. iam == 0 ) then
     ! hphp_cross
     DO ch = 1, channels_2bcross%number_confs
        bra_confs = number_2bcross(2,ch)
        IF ( bra_confs <= 0 ) cycle
        DO bra = 1, bra_confs
           i   = lookup_2bcross_configs(2,ch)%ival2(1,bra)
           b   = lookup_2bcross_configs(2,ch)%ival2(2,bra)
           DO ket = 1, bra_confs
              j   = lookup_2bcross_configs(2,ch)%ival2(1,ket)
              a   = lookup_2bcross_configs(2,ch)%ival2(2,ket)
              IF ( abs( v2b_test(i,a,j,b) - v2b_hphp_cross(ch)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "v2b_hphp: ", i, a, j, b, v2b_test(i,a,j,b), v2b_hphp_cross(ch)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO

     ! hhhh
     DO ch = 1, channels_2b%number_confs
        bra_confs = number_2b(1,ch)
        IF ( bra_confs <= 0 ) cycle
        DO bra = 1, bra_confs
           i   = lookup_2b_configs(1,ch)%ival2(1,bra)
           j   = lookup_2b_configs(1,ch)%ival2(2,bra)
           DO ket = 1, bra_confs
              k   = lookup_2b_configs(1,ch)%ival2(1,ket)
              l   = lookup_2b_configs(1,ch)%ival2(2,ket)
              IF ( abs( v2b_test(i,j,k,l) - v2b_hhhh(ch)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "v2b_hhhh: ", i, j, k, l, v2b_test(i,j,k,l), v2b_hhhh(ch)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO

     ! hphh
     DO ch = 1, channels_2b%number_confs
        bra_confs = number_2b(2,ch)
        ket_confs = number_2b(1,ch)
        IF ( bra_confs * ket_confs <= 0 ) cycle
        DO bra = 1, bra_confs
           i   = lookup_2b_configs(2,ch)%ival2(1,bra)
           a   = lookup_2b_configs(2,ch)%ival2(2,bra)
           DO ket = 1, ket_confs
              j   = lookup_2b_configs(1,ch)%ival2(1,ket)
              k   = lookup_2b_configs(1,ch)%ival2(2,ket)
              IF ( abs( v2b_test(i,a,j,k) - v2b_hphh(ch)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "v2b_hphh: ", i, a, j, k, v2b_test(i,a,j,k), v2b_hphh(ch)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO

     ! pphp
     DO ch = 1, channels_2b%number_confs
        bra_confs = number_2b(3,ch)
        ket_confs = number_2b(2,ch)
        IF ( bra_confs * ket_confs <= 0 ) cycle
        DO bra = 1, bra_confs
           a   = lookup_2b_configs(3,ch)%ival2(1,bra)
           b   = lookup_2b_configs(3,ch)%ival2(2,bra)
           DO ket = 1, ket_confs
              i   = lookup_2b_configs(2,ch)%ival2(1,ket)
              c   = lookup_2b_configs(2,ch)%ival2(2,ket)
              IF ( abs( v2b_test(a,b,i,c) - v2b_pphp(ch)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "v2b_pphp: ", a, b, i, c, v2b_test(a,b,i,c), v2b_pphp(ch)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO

     ! pphh
     DO ch = 1, channels_2b%number_confs
        bra_confs = number_2b(3,ch)
        ket_confs = number_2b(1,ch)
        IF ( bra_confs * ket_confs <= 0 ) cycle
        DO bra = 1, bra_confs
           a   = lookup_2b_configs(3,ch)%ival2(1,bra)
           b   = lookup_2b_configs(3,ch)%ival2(2,bra)
           DO ket = 1, ket_confs
              i   = lookup_2b_configs(1,ch)%ival2(1,ket)
              j   = lookup_2b_configs(1,ch)%ival2(2,ket)
              IF ( abs( v2b_test(a,b,i,j) - v2b_pphh(ch)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "v2b_pphh: ", a, b, i, j, v2b_test(a,b,i,j), v2b_pphh(ch)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO

  end IF

  IF ( test == 1 ) then
     ! pppp
     DO ch = 1, channels_2b%number_confs
        IF ( check_my_channel_v2b_pphh(ch) == 0 ) cycle        
        bra_min = mapping_v2b_pphh(iam+1,ch,2)
        bra_max = mapping_v2b_pphh(iam+1,ch,3) ! only connect to t2_pphh
        bra_confs = number_2b(3,ch)
        DO bra = bra_min, bra_max
           a   = lookup_2b_configs(3,ch)%ival2(1,bra)
           b   = lookup_2b_configs(3,ch)%ival2(2,bra)
           DO ket = 1, bra_confs
              c   = lookup_2b_configs(3,ch)%ival2(1,ket)
              d   = lookup_2b_configs(3,ch)%ival2(2,ket)
              IF ( abs( v2b_test(a,b,c,d) - v2b_pppp(ch)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "v2b_pppp: ", a, b, c, d, v2b_test(a,b,c,d), v2b_pppp(ch)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO
     
  end IF
           
end SUBROUTINE setup_tests


SUBROUTINE build_tamp_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: ch,ch1,ch2,ch3
  INTEGER :: bra_confs,ket_confs
  INTEGER :: bra,ket, bra0,ket0
  INTEGER :: a,b,c,i,j,k, cind1,kind1

  IF ( .not. ALLOCATED(t2_ccm_test) ) then
     ALLOCATE( t2_ccm_test(below_ef+1:tot_orbs, below_ef+1:tot_orbs, 1:below_ef, 1:below_ef) )
  end IF
  t2_ccm_test = 0.d0
  
  DO ch = 1, channels_2b%number_confs
     bra_confs = number_2b(3,ch)
     ket_confs = number_2b(1,ch)
     IF ( bra_confs * ket_confs <= 0 ) cycle
     DO bra = 1, bra_confs
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        b   = lookup_2b_configs(3,ch)%ival2(2,bra)
        DO ket = 1, ket_confs
           i   = lookup_2b_configs(1,ch)%ival2(1,ket)
           j   = lookup_2b_configs(1,ch)%ival2(2,ket)
           t2_ccm_test(a,b,i,j) = t2_ccm(ch)%cval(bra,ket)
           t2_ccm_test(b,a,i,j) = -1.d0 * t2_ccm(ch)%cval(bra,ket)
           t2_ccm_test(a,b,j,i) = -1.d0 * t2_ccm(ch)%cval(bra,ket)
           t2_ccm_test(b,a,j,i) = t2_ccm(ch)%cval(bra,ket)
        end DO
     end DO
  end DO

  IF ( t3_switch ) then
     IF ( .not. ALLOCATED(t3_ccm_test) ) then
        ALLOCATE( t3_ccm_test(below_ef+1:tot_orbs, below_ef+1:tot_orbs, below_ef+1:tot_orbs, 1:below_ef, 1:below_ef, 1:below_ef) )
     end IF
     t3_ccm_test = 0.d0
     
     DO ch3 = ch3_min, ch3_max
        DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
           c         = clist_t3(ch3)%ival2(cind1,1)
           ch1       = clist_t3(ch3)%ival2(cind1,2)
           bra_confs = number_2b_t3(ch3)%ival2(1,ch1)
           DO kind1 = 1, klimit_t3(ch3)
              k         = klist_t3(ch3)%ival2(kind1,1)
              ch2       = klist_t3(ch3)%ival2(kind1,2)
              ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
              DO bra  = 1, bra_confs
                 bra0 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
                 a    = lookup_2b_configs(3,ch1)%ival2(1,bra0)
                 b    = lookup_2b_configs(3,ch1)%ival2(2,bra0)
                 DO ket  = 1, ket_confs
                    ket0 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
                    i    = lookup_2b_configs(1,ch2)%ival2(1,ket0)
                    j    = lookup_2b_configs(1,ch2)%ival2(2,ket0)
                    t3_ccm_test(a,b,c,i,j,k) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket)
                    t3_ccm_test(b,a,c,i,j,k) = -1.d0 * t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket)
                    t3_ccm_test(a,b,c,j,i,k) = -1.d0 * t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket)
                    t3_ccm_test(b,a,c,j,i,k) = t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket)
                 end DO
              end DO
           end DO
        end DO
     end DO     
  end IF
  
end SUBROUTINE build_tamp_test


! <a|x|b> <-- <a|v|b> - (1/2).<kl|v|bc>.<ac|t|kl>
SUBROUTINE build_hbar1b_I2_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: a,b, k,l,c
  COMPLEX(dpc) :: v2, t2

  IF ( iam == 0 ) write(6,*) "Testing I2_1b..."
  
  ! <a|x|b>
  IF ( .not. ALLOCATED(hbar1b_I2_test) ) then
     ALLOCATE( hbar1b_I2_test(below_ef+1:tot_orbs, below_ef+1:tot_orbs) )
  end IF
  hbar1b_I2_test = 0.d0
  
  DO a = below_ef+1, tot_orbs
     DO b = below_ef+1, tot_orbs
        ! <a|x|b> <-- <a|v|b>
        hbar1b_I2_test(a,b) = hbar1b_I2_test(a,b) + fock_test(a,b)        
        ! <a|x|b> <-- -(1/2).<kl|v|bc>.<ac|t|kl>
        DO c = below_ef+1, tot_orbs
           IF ( c == a .or. c == b ) cycle
           DO k = 1, below_ef-1
              DO l = k+1, below_ef ! x 1/2
                 t2 = t2_ccm_test(a,c,k,l)
                 v2 = v2b_test(k,l,b,c)
                 hbar1b_I2_test(a,b) = hbar1b_I2_test(a,b) - t2*v2
              end DO
           end DO
        end DO
     end DO
  end DO  
  
  IF ( test == 2 ) then
     DO a = below_ef+1, tot_orbs
        DO b = below_ef+1, tot_orbs
           IF ( abs( hbar1b_I2_test(a,b) - hbar1b_I2(a,b) ) > 1.e-5 ) then
              WRITE(6,*) "hbar1b I2: ", a, b, hbar1b_I2_test(a,b), hbar1b_I2(a,b)
              stop
           end IF
        end DO
     end DO
  end IF
  
end SUBROUTINE build_hbar1b_I2_test


! <i|x|j> <-- <i|v|j> + (1/2).<ik|v|cd>.<cd|t|jk>
SUBROUTINE build_hbar1b_I3_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: i,j, c,d,k
  COMPLEX(dpc) :: v2, t2

  IF ( iam == 0 ) write(6,*) "Testing I3_1b..."
  
  ! <i|x|j>
  IF ( .not. ALLOCATED(hbar1b_I3_test) ) then
     ALLOCATE( hbar1b_I3_test(1:below_ef, 1:below_ef) )
  end IF
  hbar1b_I3_test = 0.d0
  
  DO i = 1, below_ef
     DO j = 1, below_ef
        ! <i|x|j> <-- <i|v|j>
        hbar1b_I3_test(i,j) = hbar1b_I3_test(i,j) + fock_test(i,j)
        ! <i|x|j> <-- +(1/2).<ik|v|cd>.<cd|t|jk>
        DO k = 1, below_ef
           IF ( k == i .or. k == j ) cycle
           DO c = below_ef+1, tot_orbs-1
              DO d = c+1, tot_orbs ! x 1/2
                 t2 = t2_ccm_test(c,d,j,k)
                 v2 = v2b_test(i,k,c,d)
                 hbar1b_I3_test(i,j) = hbar1b_I3_test(i,j) + t2*v2
              end DO
           end DO
        end DO
     end DO
  end DO  
  
  IF ( test == 2 ) then
     DO i = 1, below_ef
        DO j = 1, below_ef
           IF ( abs( hbar1b_I3_test(i,j) - hbar1b_I3(i,j) ) > 1.e-5 ) then
              WRITE(6,*) "hbar1b I3: ", i, j, hbar1b_I3_test(i,j), hbar1b_I3(i,j)
              stop
           end IF
        end DO
     end DO
  end IF
  
end SUBROUTINE build_hbar1b_I3_test


! <ij|x|kl> <-- <ij|v|kl> + (1/2).<ij|v|cd>.<cd|t|kl>
SUBROUTINE build_hbar2b_I4_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch, bra_confs, bra,ket
  INTEGER :: i,j,k,l, c,d
  COMPLEX(dpc) :: v2, t2

  IF ( iam == 0 ) write(6,*) "Testing I4_2b..."
  
  ! <ij|x|kl>
  IF ( .not. ALLOCATED(hbar2b_I4_test) ) then
     ALLOCATE( hbar2b_I4_test(1:below_ef, 1:below_ef, 1:below_ef, 1:below_ef) )
  end IF
  hbar2b_I4_test = 0.d0
  
  DO i = 1, below_ef-1
     DO j = i+1, below_ef ! i < j
        DO k = 1, below_ef-1
           DO l = k+1, below_ef ! k < l
              ! <ij|x|kl> <-- <ij|v|kl>
              hbar2b_I4_test(i,j,k,l) = hbar2b_I4_test(i,j,k,l) + v2b_test(i,j,k,l)
              ! <ij|x|kl> <-- +(1/2).<ij|v|cd>.<cd|t|kl>
              DO c = below_ef+1, tot_orbs-1
                 DO d = c+1, tot_orbs ! x 1/2
                    t2 = t2_ccm_test(c,d,k,l)
                    v2 = v2b_test(i,j,c,d)
                    hbar2b_I4_test(i,j,k,l) = hbar2b_I4_test(i,j,k,l) + t2*v2
                 end DO
              end DO
              hbar2b_I4_test(j,i,l,k) = hbar2b_I4_test(i,j,k,l)
              hbar2b_I4_test(j,i,k,l) = -1.d0 * hbar2b_I4_test(i,j,k,l)
              hbar2b_I4_test(i,j,l,k) = -1.d0 * hbar2b_I4_test(i,j,k,l)
           end DO
        end DO
     end DO
  end DO
  
  IF ( test == 2 ) then
     DO ch = 1, channels_2b%number_confs
        bra_confs = number_2b(1,ch)
        IF ( bra_confs <= 0 ) cycle
        DO bra = 1, bra_confs
           i   = lookup_2b_configs(1,ch)%ival2(1,bra)
           j   = lookup_2b_configs(1,ch)%ival2(2,bra)
           DO ket = 1, bra_confs
              k   = lookup_2b_configs(1,ch)%ival2(1,ket)
              l   = lookup_2b_configs(1,ch)%ival2(2,ket)
              IF ( abs( hbar2b_I4_test(i,j,k,l) - hbar2b_I4(ch)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "hbar2b I4: ", i, j, k, l, hbar2b_I4_test(i,j,k,l), hbar2b_I4(ch)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO
  end IF
  
end SUBROUTINE build_hbar2b_I4_test


! <ia|x*|jb> <-- <ia|v|jb> - (1/2).<ik|v|cb>.<ca|t|jk>
SUBROUTINE build_hbar2b_I5e_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch, bra_confs, bra,ket
  INTEGER :: i,a,j,b, k,c
  INTEGER :: Nx1,Ny1,Nz1,Tz1
  COMPLEX(dpc) :: v2, t2

  IF ( iam == 0 ) write(6,*) "Testing I5e_2b..."
  
  ! <ia|x*|jb>
  IF ( .not. ALLOCATED(hbar2b_I5e_test) ) then
     ALLOCATE( hbar2b_I5e_test(1:below_ef, below_ef+1:tot_orbs, 1:below_ef, below_ef+1:tot_orbs) )
  end IF
  hbar2b_I5e_test = 0.d0
  
  DO i = 1, below_ef
     DO a = below_ef+1, tot_orbs
        Nx1 = all_orbit%nx(i) + all_orbit%nx(a)
        Ny1 = all_orbit%ny(i) + all_orbit%ny(a)
        Nz1 = all_orbit%nz(i) + all_orbit%nz(a)
        Tz1 = all_orbit%tz(i) + all_orbit%tz(a)
        DO j = 1, below_ef
           DO b = below_ef+1, tot_orbs
              IF ( all_orbit%nx(j) + all_orbit%nx(b) /= Nx1 ) cycle
              IF ( all_orbit%ny(j) + all_orbit%ny(b) /= Ny1 ) cycle
              IF ( all_orbit%nz(j) + all_orbit%nz(b) /= Nz1 ) cycle
              IF ( all_orbit%tz(j) + all_orbit%tz(b) /= Tz1 ) cycle
              
              ! <ia|x*|jb> <-- <ia|v|jb>
              hbar2b_I5e_test(i,a,j,b) = hbar2b_I5e_test(i,a,j,b) + v2b_test(i,a,j,b)
              ! <ia|x*|jb> <-- -(1/2).<ik|v|cb>.<ca|t|jk>
              DO c = below_ef+1, tot_orbs
                 IF ( c == a .or. c == b ) cycle
                 DO k = 1, below_ef
                    IF ( k == i .or. k == j ) cycle
                    t2 = t2_ccm_test(c,a,j,k)
                    v2 = v2b_test(i,k,c,b)
                    hbar2b_I5e_test(i,a,j,b) = hbar2b_I5e_test(i,a,j,b) - 0.5d0 * t2*v2
                 end DO
              end DO
           end DO
        end DO
     end DO
  end DO
  
  IF ( test == 2 ) then
     DO ch = 1, channels_2bcross%number_confs
        IF ( check_my_channel_v2b_phhp_cross(ch) == 0 ) cycle
        bra_confs = number_2bcross(2,ch)
        DO bra = 1, bra_confs
           i   = lookup_2bcross_configs(2,ch)%ival2(1,bra)
           b   = lookup_2bcross_configs(2,ch)%ival2(2,bra)
           DO ket = 1, bra_confs
              j   = lookup_2bcross_configs(2,ch)%ival2(1,ket)
              a   = lookup_2bcross_configs(2,ch)%ival2(2,ket)
              IF ( abs( hbar2b_I5e_test(i,a,j,b) - hbar2b_I5e_cross(ch)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "hbar2b I5e_cross: ", i, a, j, b, hbar2b_I5e_test(i,a,j,b), hbar2b_I5e_cross(ch)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO
  end IF
  
end SUBROUTINE build_hbar2b_I5e_test


SUBROUTINE build_energy_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  USE chiral_potentials
  
  IMPLICIT NONE
  INTEGER :: a,b,i,j, k,c
  COMPLEX(dpc) :: v2,v3, t2,t3, sum1

  IF ( iam == 0 ) write(6,*) "Testing Ecorr..."

  sum1 = 0.d0
  
  DO a = below_ef+1, tot_orbs
     DO b = below_ef+1, tot_orbs
        IF ( b == a ) cycle
        DO i = 1, below_ef
           DO j = 1, below_ef
              IF ( j == i ) cycle
              v2 = v2b_test(i,j,a,b)
              t2 = t2_ccm_test(a,b,i,j)
              sum1 = sum1 + 0.25d0 * t2*v2
           end DO
        end DO
     end DO
  end DO
     
  IF ( tnf_approx > 1 .and. t3_switch ) then
     DO a = below_ef+1, tot_orbs
        DO b = below_ef+1, tot_orbs
           IF ( b == a ) cycle
           DO c = below_ef+1, tot_orbs
              IF ( c == a .or. c == b ) cycle
              DO i = 1, below_ef
                 DO j = 1, below_ef
                    IF ( j == i ) cycle
                    DO k = 1, below_ef
                       IF ( k == i .or. k == j ) cycle
                       v3 = v3int(i,j,k,a,b,c)
                       t3 = t3_ccm_test(a,b,c,i,j,k)
                       sum1 = sum1 + t3*v3/36.d0
                    end DO
                 end DO
              end DO
           end DO
        end DO
     end DO
  end IF

  IF ( test == 2 ) then
     WRITE(6,*) "Ecorr_test: ", sum1, ecorr
  end IF
  
end SUBROUTINE build_energy_test


! <ab|t|ij> <-- <ab|v|ij> + P(ab).<cb|t|ij>.<a|x|c> - P(ij).<ab|t|ik>.<k|x|j> + (1/2).<ab|v|cd>.<cd|t|ij> + (1/2).<ab|t|kl>.<kl|x|ij> - P(ab|ij).<ac|t|kj>.<kb|x*|ic>
! <ab|t|ij> <-- + (1/2).P(ab).<cda|t|kji>.<kb|v|cd> - (1/2).P(ij).<bca|t|kli>.<kl|v|jc> + <abc|t|ijk>.<k|v|c>
SUBROUTINE build_t2_eqn_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, bra,ket
  INTEGER :: a,b,i,j, k,l,c,d
  COMPLEX(dpc) :: denom

  IF ( iam == 0 ) write(6,*) "Testing t2..."
  
  ! <ab|t|ij>
  IF ( .not. ALLOCATED(t2_ccm_eqn_test) ) then
     ALLOCATE( t2_ccm_eqn_test(below_ef+1:tot_orbs, below_ef+1:tot_orbs, 1:below_ef, 1:below_ef) )
  end IF
  t2_ccm_eqn_test = 0.d0
  
  DO a = below_ef+1, tot_orbs-1
     DO b = a+1, tot_orbs ! a < b
        DO i = 1, below_ef-1
           DO j = i+1, below_ef ! i < j
              ! <ab|t|ij> <-- <ab|v|ij>
              t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) + v2b_test(a,b,i,j)
              ! <ab|t|ij> <-- +P(ab).<cb|t|ij>.<a|x|c>
              DO c = below_ef+1, tot_orbs
                 IF ( c == a .or. c == b ) cycle
                 t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) + hbar1b_I2_test(a,c)*t2_ccm_test(c,b,i,j)
                 t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) - hbar1b_I2_test(b,c)*t2_ccm_test(c,a,i,j)
              end DO
              ! <ab|t|ij> <-- -P(ij).<ab|t|ik>.<k|x|j>
              DO k = 1, below_ef
                 IF ( k == i .or. k == j ) cycle
                 t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) - hbar1b_I3_test(k,j)*t2_ccm_test(a,b,i,k)
                 t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) + hbar1b_I3_test(k,i)*t2_ccm_test(a,b,j,k)
              end DO
              ! <ab|t|ij> <-- +(1/2).<ab|v|cd>.<cd|t|ij>
              DO c = below_ef+1, tot_orbs-1
                 DO d = c+1, tot_orbs ! x 1/2
                    t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) + v2b_test(a,b,c,d)*t2_ccm_test(c,d,i,j)
                 end DO
              end DO
              ! <ab|t|ij> <-- +(1/2).<ab|t|kl>.<kl|x|ij>
              DO k = 1, below_ef-1
                 DO l = k+1, below_ef ! x 1/2
                    t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) + hbar2b_I4_test(k,l,i,j)*t2_ccm_test(a,b,k,l)
                 end DO
              end DO
              ! <ab|t|ij> <-- -P(ab|ij).<ac|t|kj>.<kb|x*|ic>
              DO c = below_ef+1, tot_orbs
                 DO k = 1, below_ef
                    t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) - hbar2b_I5e_test(k,b,i,c)*t2_ccm_test(a,c,k,j)
                    t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) + hbar2b_I5e_test(k,a,i,c)*t2_ccm_test(b,c,k,j)
                    t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) + hbar2b_I5e_test(k,b,j,c)*t2_ccm_test(a,c,k,i)
                    t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) - hbar2b_I5e_test(k,a,j,c)*t2_ccm_test(b,c,k,i)
                 end DO
              end DO
              IF ( t3_switch ) then
                 ! <ab|t|ij> <-- +(1/2).P(ab).<cda|t|kji>.<kb|v|cd>
                 DO k = 1, below_ef
                    IF ( k == i .or. k == j ) cycle
                    DO c = below_ef+1, tot_orbs-1
                       DO d = c+1, tot_orbs ! x 1/2
                          t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) + v2b_test(k,b,c,d)*t3_ccm_test(c,d,a,k,j,i)
                          t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) - v2b_test(k,a,c,d)*t3_ccm_test(c,d,b,k,j,i)
                       end DO
                    end DO
                 end DO
                 ! <ab|t|ij> <-- -(1/2).P(ij).<bca|t|kli>.<kl|v|jc>
                 DO c = below_ef+1, tot_orbs
                    IF ( c == a .or. c == b ) cycle
                    DO k = 1, below_ef-1
                       DO l = k+1, below_ef ! x 1/2
                          t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) - v2b_test(k,l,j,c)*t3_ccm_test(b,c,a,k,l,i)
                          t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j) + v2b_test(k,l,i,c)*t3_ccm_test(b,c,a,k,l,j)
                       end DO
                    end DO
                 end DO
              end IF
              denom = hbar1b_I3(i,i) + hbar1b_I3(j,j) - hbar1b_I2(a,a) - hbar1b_I2(b,b)
              t2_ccm_eqn_test(a,b,i,j) = t2_ccm_eqn_test(a,b,i,j)/denom
              t2_ccm_eqn_test(b,a,j,i) = t2_ccm_eqn_test(a,b,i,j)
              t2_ccm_eqn_test(b,a,i,j) = -1.d0 * t2_ccm_eqn_test(a,b,i,j)
              t2_ccm_eqn_test(a,b,j,i) = -1.d0 * t2_ccm_eqn_test(a,b,i,j)
           end DO
        end DO
     end DO
  end DO
  
  IF ( test == 3 ) then
     DO ch = 1, channels_2b%number_confs
        bra_confs = number_2b(3,ch)
        ket_confs = number_2b(1,ch)
        IF ( bra_confs <= 0 ) cycle
        IF ( ket_confs <= 0 ) cycle
        DO bra = 1, bra_confs
           a   = lookup_2b_configs(3,ch)%ival2(1,bra)
           b   = lookup_2b_configs(3,ch)%ival2(2,bra)
           DO ket = 1, ket_confs
              i   = lookup_2b_configs(1,ch)%ival2(1,ket)
              j   = lookup_2b_configs(1,ch)%ival2(2,ket)
              IF ( abs( t2_ccm_eqn_test(a,b,i,j) - t2_ccm_eqn(ch)%cval(bra,ket) ) > 1.e-10 ) then
                 WRITE(6,*) "t2_eqn: ", a, b, i, j, t2_ccm_eqn_test(a,b,i,j), t2_ccm_eqn(ch)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO
  end IF
  
end SUBROUTINE build_t2_eqn_test


! <abc|t|ijk> <-- + P(c/ab).<c|x|d>.<abd|t|ijk> - P(k/ij).<abc|t|ijl>.<l|x|k>
! <abc|t|ijk> <-- - P(c/ab|k/ij).<cd|t|ij>.<ab|v|kd> + P(c/ab|k/ij).<ab|t|lk>.<lc|v|ij>
SUBROUTINE build_t3_eqn_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE ang_mom_functions
  USE chiral_potentials
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch1,ch2,ch3
  INTEGER :: bra_confs,ket_confs
  INTEGER :: bra,ket, bra0,ket0
  INTEGER :: a,b,i,j, c,d,k,l, e,m 
  INTEGER :: Nx1,Ny1,Nz1,Tz1, cind1,kind1
  COMPLEX(dpc) :: v3, denom, sum1,sum2,sum3
  sum1 = 0.d0
  sum2 = 0.d0
  sum3 = 0.d0

  IF ( iam == 0 ) write(6,*) "Testing t3..."
  
  ! <ab|t|ij>
  IF ( .not. ALLOCATED(t3_ccm_test) ) then
     ALLOCATE( t3_ccm_test(below_ef+1:tot_orbs, below_ef+1:tot_orbs, below_ef+1:tot_orbs, 1:below_ef, 1:below_ef, 1:below_ef) )
  end IF
  t3_ccm_test = 0.d0
  
  DO a = below_ef+1, tot_orbs-2
     DO b = a+1, tot_orbs-1 ! a < b
        DO c = b+1, tot_orbs ! b < c
           IF ( cut3b(a,b,c) ) cycle
           Nx1 = all_orbit%nx(a) + all_orbit%nx(b) + all_orbit%nx(c)
           Ny1 = all_orbit%ny(a) + all_orbit%ny(b) + all_orbit%ny(c)
           Nz1 = all_orbit%nz(a) + all_orbit%nz(b) + all_orbit%nz(c)
           Tz1 = all_orbit%tz(a) + all_orbit%tz(b) + all_orbit%tz(c)
           DO i = 1, below_ef-2
              DO j = i+1, below_ef-1 ! i < j
                 DO k = j+1, below_ef ! j < k
                    IF ( cut3b(i,j,k) ) cycle
                    IF ( all_orbit%nx(i) + all_orbit%nx(j) + all_orbit%nx(k) /= Nx1 ) cycle
                    IF ( all_orbit%ny(i) + all_orbit%ny(j) + all_orbit%ny(k) /= Ny1 ) cycle
                    IF ( all_orbit%nz(i) + all_orbit%nz(j) + all_orbit%nz(k) /= Nz1 ) cycle
                    IF ( all_orbit%tz(i) + all_orbit%tz(j) + all_orbit%tz(k) /= Tz1 ) cycle

                    ! <abc|t|ijk> <-- <abd|v|ijk>
                    IF ( tnf_approx > 1 ) then
                       v3 = v3int(a,b,c,i,j,k)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + v3
                    end IF
                    ! <abc|t|ijk> <-- -P(c/ab|k/ij).<cd|t|ij>.<ab|v|kd>
                    DO d = below_ef+1, tot_orbs
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - v2b_test(a,b,k,d)*t2_ccm_test(c,d,i,j)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + v2b_test(c,b,k,d)*t2_ccm_test(a,d,i,j)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + v2b_test(a,c,k,d)*t2_ccm_test(b,d,i,j)
                       
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + v2b_test(a,b,i,d)*t2_ccm_test(c,d,k,j)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - v2b_test(c,b,i,d)*t2_ccm_test(a,d,k,j)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - v2b_test(a,c,i,d)*t2_ccm_test(b,d,k,j)
                       
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + v2b_test(a,b,j,d)*t2_ccm_test(c,d,i,k)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - v2b_test(c,b,j,d)*t2_ccm_test(a,d,i,k)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - v2b_test(a,c,j,d)*t2_ccm_test(b,d,i,k)
                    end DO
                    ! <abc|t|ijk> <-- +P(c/ab|k/ij).<ab|t|lk>.<lc|v|ij>
                    DO l = 1, below_ef
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + v2b_test(l,c,i,j)*t2_ccm_test(a,b,l,k)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - v2b_test(l,a,i,j)*t2_ccm_test(c,b,l,k)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - v2b_test(l,b,i,j)*t2_ccm_test(a,c,l,k)

                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - v2b_test(l,c,k,j)*t2_ccm_test(a,b,l,i)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + v2b_test(l,a,k,j)*t2_ccm_test(c,b,l,i)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + v2b_test(l,b,k,j)*t2_ccm_test(a,c,l,i)
                       
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - v2b_test(l,c,i,k)*t2_ccm_test(a,b,l,j)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + v2b_test(l,a,i,k)*t2_ccm_test(c,b,l,j)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + v2b_test(l,b,i,k)*t2_ccm_test(a,c,l,j)
                    end DO

                    IF ( cc_approx == 1 ) then
                       ! <abc|t|ijk> <-- <abd|v|ijk>
                       v3 = v3int(a,b,c,i,j,k)
                       t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + v3
                       ! <abc|t|ijk> <-- + P(c/ab|k/ij).<cd|t|kl>.<lab|w|dij>
                       DO l = 1, below_ef
                          DO d = below_ef+1, tot_orbs
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + t2_ccm_test(c,d,k,l)*v3int(l,a,b,d,i,j)
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - t2_ccm_test(a,d,k,l)*v3int(l,c,b,d,i,j)
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - t2_ccm_test(b,d,k,l)*v3int(l,a,c,d,i,j)

                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - t2_ccm_test(c,d,i,l)*v3int(l,a,b,d,k,j)
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + t2_ccm_test(a,d,i,l)*v3int(l,c,b,d,k,j)
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + t2_ccm_test(b,d,i,l)*v3int(l,a,c,d,k,j)
                             
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - t2_ccm_test(c,d,j,l)*v3int(l,a,b,d,i,k)
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + t2_ccm_test(a,d,j,l)*v3int(l,c,b,d,i,k)
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + t2_ccm_test(b,d,j,l)*v3int(l,a,c,d,i,k)
                          end DO
                       end DO
                       ! <abc|t|ijk> <-- + P(k/ij).<de|t|ij>.<cab|w|kde>
                       DO d = below_ef+1, tot_orbs-1
                          DO e = d+1, tot_orbs ! x 1/2
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + t2_ccm_test(d,e,i,j)*v3int(c,a,b,k,d,e)
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - t2_ccm_test(d,e,k,j)*v3int(c,a,b,i,d,e)
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - t2_ccm_test(d,e,i,k)*v3int(c,a,b,j,d,e)
                          end DO
                       end DO
                       ! <abc|t|ijk> <-- + P(c/ab).<ab|t|lm>.<clm|w|kij>
                       DO l = 1, below_ef-1
                          DO m = l+1, below_ef ! x 1/2
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) + t2_ccm_test(a,b,l,m)*v3int(c,l,m,k,i,j)
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - t2_ccm_test(c,b,l,m)*v3int(a,l,m,k,i,j)
                             t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k) - t2_ccm_test(a,c,l,m)*v3int(b,l,m,k,i,j)
                          end DO
                       end DO
                    end IF
                    
                    IF ( tnf_approx > 1 ) then
                       denom = hbar1b_I3(i,i) + hbar1b_I3(j,j) + hbar1b_I3(k,k) - hbar1b_I2(a,a) - hbar1b_I2(b,b) - hbar1b_I2(c,c)
                    else
                       denom = fock_test(i,i) + fock_test(j,j) + fock_test(k,k) - fock_test(a,a) - fock_test(b,b) - fock_test(c,c)
                    end IF
                    
                    t3_ccm_test(a,b,c,i,j,k) = t3_ccm_test(a,b,c,i,j,k)/denom
                    t3_ccm_test(b,c,a,i,j,k) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(c,a,b,i,j,k) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(b,a,c,i,j,k) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(a,c,b,i,j,k) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(c,b,a,i,j,k) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)

                    t3_ccm_test(a,b,c,j,k,i) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(b,c,a,j,k,i) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(c,a,b,j,k,i) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(b,a,c,j,k,i) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(a,c,b,j,k,i) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(c,b,a,j,k,i) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)

                    t3_ccm_test(a,b,c,k,i,j) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(b,c,a,k,i,j) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(c,a,b,k,i,j) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(b,a,c,k,i,j) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(a,c,b,k,i,j) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(c,b,a,k,i,j) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)

                    t3_ccm_test(a,b,c,j,i,k) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(b,c,a,j,i,k) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(c,a,b,j,i,k) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(b,a,c,j,i,k) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(a,c,b,j,i,k) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(c,b,a,j,i,k) = t3_ccm_test(a,b,c,i,j,k)

                    t3_ccm_test(a,b,c,i,k,j) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(b,c,a,i,k,j) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(c,a,b,i,k,j) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(b,a,c,i,k,j) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(a,c,b,i,k,j) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(c,b,a,i,k,j) = t3_ccm_test(a,b,c,i,j,k)

                    t3_ccm_test(a,b,c,k,j,i) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(b,c,a,k,j,i) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(c,a,b,k,j,i) = -1.d0 * t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(b,a,c,k,j,i) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(a,c,b,k,j,i) = t3_ccm_test(a,b,c,i,j,k)
                    t3_ccm_test(c,b,a,k,j,i) = t3_ccm_test(a,b,c,i,j,k)
                 end DO
              end DO
           end DO
        end DO
     end DO
  end DO

  IF ( test == 3 ) then
     DO ch3 = ch3_min, ch3_max
        DO cind1 = climits_t3(ch3,1), climits_t3(ch3,2)
           c         = clist_t3(ch3)%ival2(cind1,1)
           ch1       = clist_t3(ch3)%ival2(cind1,2)
           bra_confs = number_2b_t3(ch3)%ival2(1,ch1)
           DO kind1 = 1, klimit_t3(ch3)
              k         = klist_t3(ch3)%ival2(kind1,1)
              ch2       = klist_t3(ch3)%ival2(kind1,2)
              ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
              DO bra  = 1, bra_confs
                 bra0 = pp_config_t3(ch3)%ival1(ch1)%ival1(bra)
                 a    = lookup_2b_configs(3,ch1)%ival2(1,bra0)
                 b    = lookup_2b_configs(3,ch1)%ival2(2,bra0)
                 IF ( c <= b ) cycle
                 DO ket  = 1, ket_confs
                    ket0 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
                    i    = lookup_2b_configs(1,ch2)%ival2(1,ket0)
                    j    = lookup_2b_configs(1,ch2)%ival2(2,ket0)
                    IF ( k <= j ) cycle
                    IF ( abs( t3_ccm_test(a,b,c,i,j,k) - t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket) ) > 1.e-5 ) then
                       WRITE(6,*) "t3_eqn: ", a, b, c, i, j, k, t3_ccm_test(a,b,c,i,j,k), &
                            t3_ccm(ch3)%val2(cind1,kind1)%cval(bra,ket)
                       stop
                    end IF
                 end DO
              end DO
           end DO
        end DO
     end DO
  end IF
  
end SUBROUTINE build_t3_eqn_test







!
! PAEOM
!
SUBROUTINE build_paeom_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch,ch1,ch2,ch3, bra_confs,ket_confs
  INTEGER :: bra_min,bra_max, bra,ket
  INTEGER :: a,b,c,j,k, aind
  
  ! <a|r|>
  IF ( .not. ALLOCATED(r1_paeom_test) ) then
     ALLOCATE( r1_paeom_test(below_ef+1:tot_orbs) )
  end IF
  r1_paeom_test = 0.d0

  DO aind = 1, 2
     a    = r1_paeom_ind(aind)
     r1_paeom_test(a) = r1_paeom(aind)
  end DO

  ! <ab|r|j>
  IF ( .not. ALLOCATED(r2_paeom_test) ) then
     ALLOCATE( r2_paeom_test(below_ef+1:tot_orbs,below_ef+1:tot_orbs,1:below_ef) )
  end IF
  r2_paeom_test = 0.d0

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(3,ch) == 0 ) cycle
     IF ( r2_paeom_ind(2*ch) == 0 ) cycle
     bra_confs = number_2b(3,ch)
     DO bra = 1, bra_confs
        a   = lookup_2b_configs(3,ch)%ival2(1,bra)
        b   = lookup_2b_configs(3,ch)%ival2(2,bra)
        DO ket = 1, 2
           j   = r2_paeom_ind(2*(ch-1)+ket)
           r2_paeom_test(a,b,j) = r2_paeom(ch)%cval(bra,ket)
           r2_paeom_test(b,a,j) = -1.d0 * r2_paeom(ch)%cval(bra,ket)
        end DO
     end DO
  end DO

  IF ( eom_approx > 0 ) then
     
     ! <abc|r|jk>
     IF ( .not. ALLOCATED( r3_paeom_test ) ) ALLOCATE( r3_paeom_test(below_ef+1:tot_orbs, below_ef+1:tot_orbs, &
          below_ef+1:tot_orbs, 1:below_ef, 1:below_ef) )
     r3_paeom_test = 0.d0
     
     DO ch3 = ch3_paeom_min, ch3_paeom_max
        DO aind = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
           a         = clist_paeom3(ch3)%ival2(aind,1)
           ch1       = clist_paeom3(ch3)%ival2(aind,2)
           bra_min   = mapping_paeom3(ch3)%ival2(aind,1)
           bra_max   = mapping_paeom3(ch3)%ival2(aind,2)
           ch2       = ch2_paeom3(ch3)
           ket_confs = number_2b(1,ch2)
           DO bra = bra_min, bra_max
              b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
              c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
              DO ket = 1, ket_confs
                 j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
                 k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
                 r3_paeom_test(a,b,c,j,k) = r3_paeom(ch3)%val1(aind)%cval(bra,ket)
                 r3_paeom_test(a,c,b,j,k) = -1.d0 * r3_paeom(ch3)%val1(aind)%cval(bra,ket)
                 r3_paeom_test(a,b,c,k,j) = -1.d0 * r3_paeom(ch3)%val1(aind)%cval(bra,ket)
                 r3_paeom_test(a,c,b,k,j) = r3_paeom(ch3)%val1(aind)%cval(bra,ket)
              end DO
           end DO
        end DO
     end DO
     
  end IF
  
end SUBROUTINE build_paeom_test

! <a|r|> <-- <a|x|c>.<c|r|> - (1/2).<la|v|cd>.<cd|r|l> + (1/4).<kl|v|cd>.<acd|r|kl>
SUBROUTINE build_r1_paeom_eqn_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: a, aind, c,d,k,l
  COMPLEX(dpc) :: v1b, v2b, r1, r2, r3
  
  ! <a|r|>
  IF ( .not. ALLOCATED( r1_paeom_eqn_test ) ) ALLOCATE( r1_paeom_eqn_test(below_ef+1:tot_orbs) )
  r1_paeom_eqn_test = 0.d0

  DO a = below_ef+1, tot_orbs
     ! <a|r|> <-- + <a|x|c>.<c|r|>
     DO c = below_ef+1, tot_orbs
        v1b = hbar1b_I2_test(a,c)
        r1 = r1_paeom_test(c)
        r1_paeom_eqn_test(a) = r1_paeom_eqn_test(a) + v1b*r1
     end DO
     
     ! <a|r|> <-- - (1/2).<la|v|cd>.<cd|r|l>
     DO l = 1, below_ef
        DO c = below_ef+1, tot_orbs-1
           DO d = c+1, tot_orbs ! x 1/2
              v2b = v2b_test(l,a,c,d)
              r2 = r2_paeom_test(c,d,l)
              r1_paeom_eqn_test(a) = r1_paeom_eqn_test(a) - v2b*r2
           end DO
        end DO
     end DO
     
     ! <a|r|> <-- + (1/4).<kl|v|cd>.<acd|r|kl>
     IF ( eom_approx > 0 ) then
        DO k = 1, below_ef-1
           DO l = k+1, below_ef ! x 1/2
              DO c = below_ef+1, tot_orbs-1
                 DO d = c+1, tot_orbs ! x 1/2
                    v2b = v2b_test(k,l,c,d)
                    r3 = r3_paeom_test(a,c,d,k,l)
                    r1_paeom_eqn_test(a) = r1_paeom_eqn_test(a) + v2b*r3
                 end DO
              end DO
           end DO
        end DO
     end IF
  end DO
  
  
  IF ( test == 4 ) then
     
     IF ( iam == 0 ) write(6,*) ' Testing PA-EOM right 1p...'
     DO aind = 1, 2
        a    = r1_paeom_ind(aind)
        r1_paeom_test(a) = r1_paeom(aind)
        IF ( abs( r1_paeom_eqn(aind) - r1_paeom_eqn_test(a) ) > 1.e-5 ) then
           WRITE(6,*) "r1_eqn ", a, r1_paeom_eqn(aind), r1_paeom_eqn_test(a)
           stop
        end IF
     end DO
     
  end IF
  
end SUBROUTINE build_r1_paeom_eqn_test

! <ab|r|j> <-- - <ab|x|jc>.<c|r|> + P(ab).<a|x|c>.<cb|r|j> - <k|x|j>.<ab|r|k> + (1/2).<ab|x|cd>.<cd|r|j> - P(ab).<ka|x|jc>.<cb|r|k> - <k|z3b|>.<ab|t|kj>
!              + (1/2).P(ab).<kb|v|cd>.<acd|r|kj> - (1/2).<kl|v|jc>.<abc|r|kl>
SUBROUTINE build_r2_paeom_eqn_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch,bra_confs,bra,ket
  INTEGER :: a,b,j, c,d,k,l
  COMPLEX(dpc) :: v1b, v2b, r1, r2, r3, t2
  
  ! <ab|r|j>
  IF ( .not. ALLOCATED( r2_paeom_eqn_test ) ) ALLOCATE( r2_paeom_eqn_test(below_ef+1:tot_orbs,below_ef+1:tot_orbs,1:below_ef) )
  r2_paeom_eqn_test = 0.d0

  DO a = below_ef+1, tot_orbs-1
     DO b = a+1, tot_orbs ! a <-> b
        DO j = 1, below_ef

           ! <ab|r|j> <-- - <ab|x|jc>.<c|r|>
           DO c = below_ef+1, tot_orbs
              v2b = hbar2b_I6_test(a,b,j,c)
              r1 = r1_paeom_test(c)
              r2_paeom_eqn_test(a,b,j) = r2_paeom_eqn_test(a,b,j) - v2b*r1
           end DO
           
           ! <ab|r|j> <-- + P(ab).<a|x|c>.<cb|r|j>
           DO c = below_ef+1, tot_orbs
              v1b = hbar1b_I2_test(a,c)
              r2 = r2_paeom_test(c,b,j)
              r2_paeom_eqn_test(a,b,j) = r2_paeom_eqn_test(a,b,j) + v1b*r2
              v1b = hbar1b_I2_test(b,c)
              r2 = r2_paeom_test(c,a,j)
              r2_paeom_eqn_test(a,b,j) = r2_paeom_eqn_test(a,b,j) - v1b*r2
           end DO
           
           ! <ab|r|j> <-- - <k|x|j>.<ab|r|k>
           DO k = 1, below_ef
              v1b = hbar1b_I3_test(k,j)
              r2 = r2_paeom_test(a,b,k)
              r2_paeom_eqn_test(a,b,j) = r2_paeom_eqn_test(a,b,j) - v1b*r2
           end DO

           ! <ab|r|j> <-- + (1/2).<ab|x|cd>.<cd|r|j>
           DO c = below_ef+1, tot_orbs-1
              DO d = c+1, tot_orbs ! x 1/2
                 v2b = hbar2b_I3_test(a,b,c,d)
                 r2 = r2_paeom_test(c,d,j)
                 r2_paeom_eqn_test(a,b,j) = r2_paeom_eqn_test(a,b,j) + v2b*r2
              end DO
           end DO
           
           ! <ab|r|j> <-- - P(ab).<ka|x|jc>.<cb|r|k>
           DO c = below_ef+1, tot_orbs
              DO k = 1, below_ef
                 v2b = hbar2b_I5_test(k,a,j,c)
                 r2 = r2_paeom_test(c,b,k)
                 r2_paeom_eqn_test(a,b,j) = r2_paeom_eqn_test(a,b,j) - v2b*r2
                 v2b = hbar2b_I5_test(k,b,j,c)
                 r2 = r2_paeom_test(c,a,k)
                 r2_paeom_eqn_test(a,b,j) = r2_paeom_eqn_test(a,b,j) + v2b*r2
              end DO
           end DO
           
           ! <ab|r|j> <-- - <k|z3b|>.<ab|t|kj>
           DO k = 1, below_ef
              v1b = hbar3b_paeom_I3_test(k)
              t2 = t2_ccm_test(a,b,k,j)
              r2_paeom_eqn_test(a,b,j) = r2_paeom_eqn_test(a,b,j) - v1b*t2
           end DO
                      
           IF ( eom_approx > 0 ) then
              ! <ab|r|j> <-- + (1/2).P(ab).<kb|v|cd>.<acd|r|kj>
              DO k = 1, below_ef
                 DO c = below_ef+1, tot_orbs-1
                    DO d = c+1, tot_orbs ! x 1/2
                       v2b = v2b_test(k,b,c,d)
                       r3 = r3_paeom_test(a,c,d,k,j)
                       r2_paeom_eqn_test(a,b,j) = r2_paeom_eqn_test(a,b,j) + v2b*r3
                       v2b = v2b_test(k,a,c,d)
                       r3 = r3_paeom_test(b,c,d,k,j)
                       r2_paeom_eqn_test(a,b,j) = r2_paeom_eqn_test(a,b,j) - v2b*r3
                    end DO
                 end DO
              end DO
              
              ! <ab|r|j> <-- - (1/2).<kl|v|jc>.<abc|r|kl>
              DO k = 1, below_ef-1
                 DO l = k+1, below_ef ! x 1/2
                    DO c = below_ef+1, tot_orbs
                       v2b = v2b_test(k,l,j,c)
                       r3 = r3_paeom_test(a,b,c,k,l)
                       r2_paeom_eqn_test(a,b,j) = r2_paeom_eqn_test(a,b,j) - v2b*r3
                    end DO
                 end DO
              end DO
           end IF
           
           r2_paeom_eqn_test(b,a,j) = -1.d0 * r2_paeom_eqn_test(a,b,j)
           
        end DO
     end DO
  end DO
  

  IF ( test == 4 ) then

     IF ( iam == 0 ) write(6,*) ' Testing PA-EOM right 2p1h...'
     DO ch = 1, channels_2b%number_confs
        IF ( number_2b(3,ch) == 0 ) cycle
        IF ( r2_paeom_ind(2*ch) == 0 ) cycle
        bra_confs = number_2b(3,ch)
        DO bra = 1, bra_confs
           a   = lookup_2b_configs(3,ch)%ival2(1,bra)
           b   = lookup_2b_configs(3,ch)%ival2(2,bra)
           DO ket = 1, 2
              j   = r2_paeom_ind(2*(ch-1)+ket)
              IF ( abs( r2_paeom_eqn(ch)%cval(bra,ket) - r2_paeom_eqn_test(a,b,j) ) > 1.e-5 ) then
                 WRITE(6,*) "r2_eqn ", a, b, j, r2_paeom_eqn(ch)%cval(bra,ket), r2_paeom_eqn_test(a,b,j)
                 stop
              end IF
           end DO
        end DO
     end DO
     
  end IF

end SUBROUTINE build_r2_paeom_eqn_test

! <abc|r|jk> <-- + P(bc/a).<a|x|d>.<dbc|r|jk> - P(jk).<l|x|k>.<abc|r|jl> - P(bc/a|jk).<bc|v|kd>.<ad|r|j> - P(bc/a).<la|x|jk>.<bc|r|l>
SUBROUTINE build_r3_paeom_eqn_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch3,ch1,ch2, ket_confs
  INTEGER :: bra_min,bra_max, bra,ket, aind
  INTEGER :: Nx1,Ny1,Nz1,Tz1
  INTEGER :: a,b,c,j,k, d,l
  COMPLEX(dpc) :: v2b, r2, r3
  
  ! <abc|r|jk>
  DO a = below_ef+1, tot_orbs-2
     DO b = a+1, tot_orbs-1 ! a < b
        DO c = b+1, tot_orbs ! b < c
           Nx1 = all_orbit%nx(a) + all_orbit%nx(b) + all_orbit%nx(c)
           Ny1 = all_orbit%ny(a) + all_orbit%ny(b) + all_orbit%ny(c)
           Nz1 = all_orbit%nz(a) + all_orbit%nz(b) + all_orbit%nz(c)
           Tz1 = all_orbit%tz(a) + all_orbit%tz(b) + all_orbit%tz(c)
           DO j = 1, below_ef-1
              DO k = j+1, below_ef ! j < k
                 IF ( Nx1 - all_orbit%nx(j) - all_orbit%nx(k) /= nx_paeom ) cycle
                 IF ( Ny1 - all_orbit%ny(j) - all_orbit%ny(k) /= ny_paeom ) cycle
                 IF ( Nz1 - all_orbit%nz(j) - all_orbit%nz(k) /= nz_paeom ) cycle
                 IF ( Tz1 - all_orbit%tz(j) - all_orbit%tz(k) /= tz_paeom ) cycle

                 ! <abc|r|jk> <-- + P(bc/a).<a|x|d>.<dbc|r|jk> - P(jk).<l|x|k>.<abc|r|jl>
                 r3 = r3_paeom_test(a,b,c,j,k)
                 r3_paeom_test(a,b,c,j,k) = 0.d0
                 r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) + hbar1b_I2_test(a,a) * r3
                 r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) + hbar1b_I2_test(b,b) * r3
                 r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) + hbar1b_I2_test(c,c) * r3
                 r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) - hbar1b_I3_test(j,j) * r3
                 r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) - hbar1b_I3_test(k,k) * r3

                 ! <abc|r|jk> <-- - P(bc/a|jk).<bc|v|kd>.<ad|r|j>
                 DO d = below_ef+1, tot_orbs
                    v2b = v2b_test(b,c,k,d)
                    r2 = r2_paeom_test(a,d,j)
                    r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) - v2b*r2
                    v2b = v2b_test(a,c,k,d)
                    r2 = r2_paeom_test(b,d,j)
                    r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) + v2b*r2
                    v2b = v2b_test(b,a,k,d)
                    r2 = r2_paeom_test(c,d,j)
                    r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) + v2b*r2
                    v2b = v2b_test(b,c,j,d)
                    r2 = r2_paeom_test(a,d,k)
                    r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) + v2b*r2
                    v2b = v2b_test(a,c,j,d)
                    r2 = r2_paeom_test(b,d,k)
                    r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) - v2b*r2
                    v2b = v2b_test(b,a,j,d)
                    r2 = r2_paeom_test(c,d,k)
                    r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) - v2b*r2
                 end DO

                 ! <abc|r|jk> <-- - P(bc/a).<la|x|jk>.<bc|r|l>
                 DO l = 1, below_ef
                    v2b = v2b_test(l,a,j,k)
                    r2 = r2_paeom_test(b,c,l)
                    r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) - v2b*r2
                    v2b = v2b_test(l,b,j,k)
                    r2 = r2_paeom_test(a,c,l)
                    r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) + v2b*r2
                    v2b = v2b_test(l,c,j,k)
                    r2 = r2_paeom_test(b,a,l)
                    r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k) + v2b*r2
                 end DO
                 
                 ! r3_paeom_test(a,b,c,j,k) = r3_paeom_test(a,b,c,j,k)
                 r3_paeom_test(b,c,a,j,k) = r3_paeom_test(a,b,c,j,k)
                 r3_paeom_test(c,a,b,j,k) = r3_paeom_test(a,b,c,j,k)
                 r3_paeom_test(b,a,c,j,k) = -1.d0 * r3_paeom_test(a,b,c,j,k)
                 r3_paeom_test(a,c,b,j,k) = -1.d0 * r3_paeom_test(a,b,c,j,k)
                 r3_paeom_test(c,b,a,j,k) = -1.d0 * r3_paeom_test(a,b,c,j,k)
                 r3_paeom_test(a,b,c,k,j) = -1.d0 * r3_paeom_test(a,b,c,j,k)
                 r3_paeom_test(b,c,a,k,j) = -1.d0 * r3_paeom_test(a,b,c,j,k)
                 r3_paeom_test(c,a,b,k,j) = -1.d0 * r3_paeom_test(a,b,c,j,k)
                 r3_paeom_test(b,a,c,k,j) = r3_paeom_test(a,b,c,j,k)
                 r3_paeom_test(a,c,b,k,j) = r3_paeom_test(a,b,c,j,k)
                 r3_paeom_test(c,b,a,k,j) = r3_paeom_test(a,b,c,j,k)

              end DO
           end DO
        end DO
     end DO
  end DO
  

  IF ( test == 4 ) then

     IF ( iam == 0 ) write(6,*) ' Testing PA-EOM right 3p2h...'
     DO ch3 = ch3_paeom_min, ch3_paeom_max
        DO aind = climits_paeom3(ch3,1), climits_paeom3(ch3,2)
           a         = clist_paeom3(ch3)%ival2(aind,1)
           ch1       = clist_paeom3(ch3)%ival2(aind,2)
           bra_min   = mapping_paeom3(ch3)%ival2(aind,1)
           bra_max   = mapping_paeom3(ch3)%ival2(aind,2)
           ch2       = ch2_paeom3(ch3)
           ket_confs = number_2b(1,ch2)
           DO bra = bra_min, bra_max
              b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
              c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
              ! IF ( a >= b ) cycle
              DO ket = 1, ket_confs
                 j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
                 k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
                 IF ( abs( r3_paeom(ch3)%val1(aind)%cval(bra,ket) - r3_paeom_test(a,b,c,j,k) ) > 1.e-5 ) then
                    WRITE(6,*) "r3_eqn ", a, b, c, j, k, r3_paeom(ch3)%val1(aind)%cval(bra,ket), &
                         r3_paeom_test(a,b,c,j,k)
                    stop
                 end IF
              end DO
           end DO
        end DO
     end DO
     
  end IF
  
end SUBROUTINE build_r3_paeom_eqn_test


! <ab|x|cd> <-- <ab|v|cd> + (1/2).<ab|t|kl>.<kl|v|cd>
SUBROUTINE build_hbar2b_I3pa_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch, bra_confs, bra,ket
  INTEGER :: a,b,c,d, k,l
  COMPLEX(dpc) :: v2, t2

  IF ( iam == 0 ) write(6,*) "Testing I3_2b..."
  
  ! <ab|x|cd>
  IF ( .not. ALLOCATED(hbar2b_I3_test) ) then
     ALLOCATE( hbar2b_I3_test(below_ef+1:tot_orbs, below_ef+1:tot_orbs, below_ef+1:tot_orbs, below_ef+1:tot_orbs) )
  end IF
  hbar2b_I3_test = 0.d0
  
  DO a = below_ef+1, tot_orbs-1
     DO b = a+1, tot_orbs ! a < b
        DO c = below_ef+1, tot_orbs-1
           DO d = c+1, tot_orbs ! c < d
              ! <ab|x|cd> <-- <ab|v|cd>
              hbar2b_I3_test(a,b,c,d) = hbar2b_I3_test(a,b,c,d) + v2b_test(a,b,c,d)
              ! <ab|x|cd> <-- +(1/2).<ab|t|kl>.<kl|v|cd>
              DO k = 1, below_ef-1
                 DO l = k+1, below_ef ! x 1/2
                    t2 = t2_ccm_test(a,b,k,l)
                    v2 = v2b_test(k,l,c,d)
                    hbar2b_I3_test(a,b,c,d) = hbar2b_I3_test(a,b,c,d) + t2*v2
                 end DO
              end DO
              hbar2b_I3_test(b,a,d,c) = hbar2b_I3_test(a,b,c,d)
              hbar2b_I3_test(b,a,c,d) = -1.d0 * hbar2b_I3_test(a,b,c,d)
              hbar2b_I3_test(a,b,d,c) = -1.d0 * hbar2b_I3_test(a,b,c,d)
           end DO
        end DO
     end DO
  end DO

  
  IF ( test == 4 ) then
     DO ch = 1, channels_2b%number_confs
        IF ( number_2b(3,ch) == 0 ) cycle
        IF ( r2_paeom_ind(2*ch) == 0 ) cycle
        bra_confs = number_2b(3,ch)
        DO bra = 1, bra_confs
           a   = lookup_2b_configs(3,ch)%ival2(1,bra)
           b   = lookup_2b_configs(3,ch)%ival2(2,bra)
           DO ket = 1, bra_confs
              c   = lookup_2b_configs(3,ch)%ival2(1,ket)
              d   = lookup_2b_configs(3,ch)%ival2(2,ket)
              IF ( abs( hbar2b_I3_test(a,b,c,d) - hbar2b_I3(ch)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "hbar2b I3: ", a, b, c, d, hbar2b_I3_test(a,b,c,d), hbar2b_I3(ch)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO
  end IF
  
end SUBROUTINE build_hbar2b_I3pa_test

! <ia|x|jb> <-- <ia|v|jb> - <ik|v|cb>.<ca|t|jk>
SUBROUTINE build_hbar2b_I5pa_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch2, ket_confs, bra,ket
  INTEGER :: i,a,j,b, k,c
  COMPLEX(dpc) :: v2, t2

  IF ( iam == 0 ) write(6,*) "Testing I5_2b..."
  
  ! <ia|x|jb>
  IF ( .not. ALLOCATED(hbar2b_I5_test) ) then
     ALLOCATE( hbar2b_I5_test(1:below_ef, below_ef+1:tot_orbs, 1:below_ef, below_ef+1:tot_orbs) )
  end IF
  hbar2b_I5_test = 0.d0
  
  DO i = 1, below_ef
     DO a = below_ef+1, tot_orbs
        DO j = 1, below_ef
           DO b = below_ef+1, tot_orbs
              ! <ia|x|jb> <-- <ia|v|jb>
              hbar2b_I5_test(i,a,j,b) = hbar2b_I5_test(i,a,j,b) + v2b_test(i,a,j,b)
              ! <ia|x|jb> <-- -<ik|v|cb>.<ca|t|jk>
              DO c = below_ef+1, tot_orbs
                 IF ( c == a .or. c == b ) cycle
                 DO k = 1, below_ef
                    IF ( k == i .or. k == j ) cycle
                    t2 = t2_ccm_test(c,a,j,k)
                    v2 = v2b_test(i,k,c,b)
                    hbar2b_I5_test(i,a,j,b) = hbar2b_I5_test(i,a,j,b) - t2*v2
                 end DO
              end DO
           end DO
        end DO
     end DO
  end DO

  
  IF ( test == 4 ) then
     DO ch2 = 1, channels_2bcross%number_confs
        IF ( number_2bcross(2,ch2) == 0 ) cycle
        IF ( r2_paeom_ind_cross(2*ch2) == 0 ) cycle
        ket_confs = number_2bcross(2,ch2)
        DO bra = 1, ket_confs
           i   = lookup_2bcross_configs(2,ch2)%ival2(1,bra)
           b   = lookup_2bcross_configs(2,ch2)%ival2(2,bra)
           DO ket = 1, ket_confs
              j   = lookup_2bcross_configs(2,ch2)%ival2(1,ket)
              a   = lookup_2bcross_configs(2,ch2)%ival2(2,ket)
              IF ( abs( hbar2b_I5_test(i,a,j,b) - hbar2b_I5_cross(ch2)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "hbar2b I5_cross: ", i, a, j, b, hbar2b_I5_test(i,a,j,b), hbar2b_I5_cross(ch2)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO
  end IF
  
end SUBROUTINE build_hbar2b_I5pa_test

! <ab|x|ic> <-- <ab|v|ic> + P(ab).<ad|t|ik>.<kb|v|dc> + (1/2).<ab|t|kl>.<kl|v|ic>
SUBROUTINE build_hbar2b_I6pa_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, bra,ket
  INTEGER :: a,b,i,c, k,l,d
  COMPLEX(dpc) :: v2, t2

  IF ( iam == 0 ) write(6,*) "Testing I6_2b..."
  
  ! <ab|x|ic>
  IF ( .not. ALLOCATED(hbar2b_I6_test) ) then
     ALLOCATE( hbar2b_I6_test(below_ef+1:tot_orbs, below_ef+1:tot_orbs, 1:below_ef, below_ef+1:tot_orbs) )
  end IF
  hbar2b_I6_test = 0.d0

  DO a = below_ef+1, tot_orbs-1
     DO b = a+1, tot_orbs ! a < b
        DO i = 1, below_ef
           DO c = below_ef+1, tot_orbs
              ! <ab|x|ic> <-- <ab|v|ic>
              hbar2b_I6_test(a,b,i,c) = hbar2b_I6_test(a,b,i,c) + v2b_test(a,b,i,c)
              ! <ab|x|ic> <-- +(1/2).<ab|t|kl>.<kl|v|ic>
              DO k = 1, below_ef-1
                 DO l = k+1, below_ef ! x 1/2
                    t2 = t2_ccm_test(a,b,k,l)
                    v2 = v2b_test(k,l,i,c)
                    hbar2b_I6_test(a,b,i,c) = hbar2b_I6_test(a,b,i,c) + t2*v2
                 end DO
              end DO
              ! <ab|x|ic> <-- +P(ab).<ad|t|ik>.<kb|v|dc>
              DO d = below_ef+1, tot_orbs
                 IF ( d == c ) cycle
                 DO k = 1, below_ef
                    IF ( k == i ) cycle
                    t2 = t2_ccm_test(a,d,i,k)
                    v2 = v2b_test(k,b,d,c)
                    hbar2b_I6_test(a,b,i,c) = hbar2b_I6_test(a,b,i,c) + t2*v2
                    t2 = t2_ccm_test(b,d,i,k)
                    v2 = v2b_test(k,a,d,c)
                    hbar2b_I6_test(a,b,i,c) = hbar2b_I6_test(a,b,i,c) - t2*v2
                 end DO
              end DO
              hbar2b_I6_test(b,a,i,c) = -1.d0 * hbar2b_I6_test(a,b,i,c)
           end DO
        end DO
     end DO
  end DO

  
  IF ( test == 4 ) then
     DO ch = 1, channels_2b%number_confs
        IF ( number_2b(3,ch) == 0 ) cycle
        IF ( r2_paeom_ind(2*ch) == 0 ) cycle
        bra_confs = number_2b(3,ch)
        ket_confs = number_2b(2,ch)
        IF ( ket_confs <= 0 ) cycle
        DO bra = 1, bra_confs
           a   = lookup_2b_configs(3,ch)%ival2(1,bra)
           b   = lookup_2b_configs(3,ch)%ival2(2,bra)
           DO ket = 1, ket_confs
              i   = lookup_2b_configs(2,ch)%ival2(1,ket)
              c   = lookup_2b_configs(2,ch)%ival2(2,ket)
              IF ( abs( hbar2b_I6_test(a,b,i,c) - hbar2b_I6(ch)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "hbar2b I6: ", a, b, i, c, hbar2b_I6_test(a,b,i,c), hbar2b_I6(ch)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO
  end IF
  
end SUBROUTINE build_hbar2b_I6pa_test

! <i|z3b|> <-- (1/2).<ik|v|cd>.<cd|r|k>
SUBROUTINE build_hbar3b_paeom_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: i, k,c,d
  COMPLEX(dpc) :: v2, r2

  IF ( iam == 0 ) write(6,*) "Testing I3_3b_paeom..."
  
  ! <i|z3b|>
  IF ( .not. ALLOCATED(hbar3b_paeom_I3_test) ) then
     ALLOCATE( hbar3b_paeom_I3_test(1:below_ef) )
  end IF
  hbar3b_paeom_I3_test = 0.d0

  DO i = 1, below_ef

     ! <i|z3b|> <-- (1/2).<ik|v|cd>.<cd|r|k>
     DO k = 1, below_ef
        IF ( k == i ) cycle
        DO c = below_ef+1, tot_orbs-1
           DO d = c+1, tot_orbs ! x 1/2
              r2 = r2_paeom_test(c,d,k)
              v2 = v2b_test(i,k,c,d)
              hbar3b_paeom_I3_test(i) = hbar3b_paeom_I3_test(i) + r2*v2
           end DO
        end DO
     end DO

  end DO
  
  
  IF ( test == 4 ) then
     DO i = 1, below_ef
        IF ( abs( hbar3b_paeom_I3_test(i) - hbar3b_paeom_I3(i) ) > 1.e-5 ) then
           WRITE(6,*) "hbar3b I3: ", i, hbar3b_paeom_I3_test(i), hbar3b_paeom_I3(i)
           stop
        end IF
     end DO
  end IF
  
end SUBROUTINE build_hbar3b_paeom_test



!
! PREOM
!
SUBROUTINE build_preom_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch,ch1,ch2,ch3, bra_confs,ket_confs
  INTEGER :: ket_min,ket_max, bra,ket
  INTEGER :: i,j,k,b,c, iind
  
  ! <|r|i>
  IF ( .not. ALLOCATED(r1_preom_test) ) then
     ALLOCATE( r1_preom_test(1:below_ef) )
  end IF
  r1_preom_test = 0.d0

  DO iind = 1, 2
     i    = r1_preom_ind(iind)
     r1_preom_test(i) = r1_preom(iind)
  end DO

  ! <b|r|ij>
  IF ( .not. ALLOCATED(r2_preom_test) ) then
     ALLOCATE( r2_preom_test(below_ef+1:tot_orbs,1:below_ef,1:below_ef) )
  end IF
  r2_preom_test = 0.d0

  DO ch = 1, channels_2b%number_confs
     IF ( number_2b(1,ch) == 0 ) cycle
     IF ( r2_preom_ind(2*ch) == 0 ) cycle
     ket_confs = number_2b(1,ch)
     DO ket = 1, ket_confs
        i   = lookup_2b_configs(1,ch)%ival2(1,ket)
        j   = lookup_2b_configs(1,ch)%ival2(2,ket)
        DO bra = 1, 2
           b   = r2_preom_ind(2*(ch-1)+bra)
           r2_preom_test(b,i,j) = r2_preom(ch)%cval(bra,ket)
           r2_preom_test(b,j,i) = -1.d0 * r2_preom(ch)%cval(bra,ket)
        end DO
     end DO
  end DO

  IF ( eom_approx > 0 ) then
     
     ! <bc|r|ijk>
     IF ( .not. ALLOCATED( r3_preom_test ) ) ALLOCATE( r3_preom_test(below_ef+1:tot_orbs, below_ef+1:tot_orbs, &
          1:below_ef, 1:below_ef, 1:below_ef) )
     r3_preom_test = 0.d0
     
     DO ch3 = ch3_preom_min, ch3_preom_max
        DO iind = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
           i         = klist_preom3(ch3)%ival2(iind,1)
           ch2       = klist_preom3(ch3)%ival2(iind,2)
           ket_min   = mapping_preom3(ch3)%ival2(iind,1)
           ket_max   = mapping_preom3(ch3)%ival2(iind,2)
           ch1       = ch1_preom3(ch3)
           bra_confs = number_2b(3,ch1)
           DO bra = 1, bra_confs
              b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
              c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
              DO ket = ket_min, ket_max
                 j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
                 k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
                 r3_preom_test(b,c,i,j,k) = r3_preom(ch3)%val1(iind)%cval(bra,ket)
                 r3_preom_test(c,b,i,j,k) = -1.d0 * r3_preom(ch3)%val1(iind)%cval(bra,ket)
                 r3_preom_test(b,c,i,k,j) = -1.d0 * r3_preom(ch3)%val1(iind)%cval(bra,ket)
                 r3_preom_test(c,b,i,k,j) = r3_preom(ch3)%val1(iind)%cval(bra,ket)
              end DO
           end DO
        end DO
     end DO
     
  end IF
  
end SUBROUTINE build_preom_test

! <|r|i> <-- -<k|x|i>.<|r|k> - (1/2).<kl|v|id>.<d|r|kl> + (1/4).<kl|v|cd>.<cd|r|ikl>
SUBROUTINE build_r1_preom_eqn_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel
  
  IMPLICIT NONE
  INTEGER :: i, iind, c,d,k,l
  COMPLEX(dpc) :: v1b, v2b, r1, r2, r3
  
  ! <|r|i>
  IF ( .not. ALLOCATED( r1_preom_eqn_test ) ) ALLOCATE( r1_preom_eqn_test(1:below_ef) )
  r1_preom_eqn_test = 0.d0

  DO i = 1, below_ef
     ! <|r|i> <-- - <k|x|i>.<|r|k>
     DO k = 1, below_ef
        v1b = hbar1b_I3_test(k,i)
        r1 = r1_preom_test(k)
        r1_preom_eqn_test(i) = r1_preom_eqn_test(i) - v1b*r1
     end DO

     ! <|r|i> <-- - (1/2).<kl|v|id>.<d|r|kl>
     DO d = below_ef+1, tot_orbs
        DO k = 1, below_ef-1
           DO l = k+1, below_ef ! x 1/2
              v2b = v2b_test(k,l,i,d)
              r2 = r2_preom_test(d,k,l)
              r1_preom_eqn_test(i) = r1_preom_eqn_test(i) - v2b*r2
           end DO
        end DO
     end DO
     
     ! <|r|i> <-- + (1/4).<kl|v|cd>.<cd|r|ikl>
     IF ( eom_approx > 0 ) then
        DO c = below_ef+1, tot_orbs-1
           DO d = c+1, tot_orbs ! x 1/2
              DO k = 1, below_ef-1
                 DO l = k+1, below_ef ! x 1/2
                    v2b = v2b_test(k,l,c,d)
                    r3 = r3_preom_test(c,d,i,k,l)
                    r1_preom_eqn_test(i) = r1_preom_eqn_test(i) + v2b*r3
                 end DO
              end DO
           end DO
        end DO
     end IF
  end DO
  
  
  IF ( test == 4 ) then
     
     IF ( iam == 0 ) write(6,*) ' Testing PR-EOM right 1h...'
     DO iind = 1, 2
        i    = r1_preom_ind(iind)
        r1_preom_test(i) = r1_preom(iind)
        IF ( abs( r1_preom_eqn(iind) - r1_preom_eqn_test(i) ) > 1.e-5 ) then
           WRITE(6,*) "r1_eqn ", i, r1_preom_eqn(iind), r1_preom_eqn_test(i)
           stop
        end IF
     end DO
     
  end IF
  
end SUBROUTINE build_r1_preom_eqn_test

! <b|r|ij> <-- - <kb|x|ij>.<|r|k> - P(ij).<k|x|i>.<b|r|kj> + <b|x|c>.<c|r|ij> + (1/2).<kl|x|ij>.<b|r|kl> - P(ij).<kb|x|ic>.<c|r|kj> + <|z3b|c>.<cb|t|ij>
!              - (1/2).P(ij).<kl|v|jc>.<bc|r|ikl> - (1/2).<kb|v|cd>.<cd|r|ijk>
SUBROUTINE build_r2_preom_eqn_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch,ket_confs,bra,ket
  INTEGER :: i,j,b, c,d,k,l
  COMPLEX(dpc) :: v1b, v2b, r1, r2, r3, t2
  
  ! <b|r|ij>
  IF ( .not. ALLOCATED( r2_preom_eqn_test ) ) ALLOCATE( r2_preom_eqn_test(below_ef+1:tot_orbs,1:below_ef,1:below_ef) )
  r2_preom_eqn_test = 0.d0

  DO b = below_ef+1, tot_orbs
     DO i = 1, below_ef-1
        DO j = i+1, below_ef ! i <-> j

           ! <b|r|ij> <-- - <kb|x|ij>.<|r|k>
           DO k = 1, below_ef
              v2b = hbar2b_I7_test(k,b,i,j)
              r1 = r1_preom_test(k)
              r2_preom_eqn_test(b,i,j) = r2_preom_eqn_test(b,i,j) - v2b*r1
           end DO

           ! <b|r|ij> <-- - P(ij).<k|x|i>.<b|r|kj>
           DO k = 1, below_ef
              v1b = hbar1b_I3_test(k,i)
              r2 = r2_preom_test(b,k,j)
              r2_preom_eqn_test(b,i,j) = r2_preom_eqn_test(b,i,j) - v1b*r2
              v1b = hbar1b_I3_test(k,j)
              r2 = r2_preom_test(b,k,i)
              r2_preom_eqn_test(b,i,j) = r2_preom_eqn_test(b,i,j) + v1b*r2
           end DO

           ! <b|r|ij> <-- + <b|x|c>.<c|r|ij>
           DO c = below_ef+1, tot_orbs
              v1b = hbar1b_I2_test(b,c)
              r2 = r2_preom_test(c,i,j)
              r2_preom_eqn_test(b,i,j) = r2_preom_eqn_test(b,i,j) + v1b*r2
           end DO

           ! <b|r|ij> <-- + (1/2).<kl|x|ij>.<b|r|kl>
           DO k = 1, below_ef-1
              DO l = k+1, below_ef ! x 1/2
                 v2b = hbar2b_I4_test(k,l,i,j)
                 r2 = r2_preom_test(b,k,l)
                 r2_preom_eqn_test(b,i,j) = r2_preom_eqn_test(b,i,j) + v2b*r2
              end DO
           end DO

           ! <b|r|ij> <-- - P(ij).<kb|x|ic>.<c|r|kj>
           DO c = below_ef+1, tot_orbs
              DO k = 1, below_ef
                 v2b = hbar2b_I5_test(k,b,i,c)
                 r2 = r2_preom_test(c,k,j)
                 r2_preom_eqn_test(b,i,j) = r2_preom_eqn_test(b,i,j) - v2b*r2
                 v2b = hbar2b_I5_test(k,b,j,c)
                 r2 = r2_preom_test(c,k,i)
                 r2_preom_eqn_test(b,i,j) = r2_preom_eqn_test(b,i,j) + v2b*r2
              end DO
           end DO

           ! <b|r|ij> <-- + <|z3b|c>.<cb|t|ij>
           DO c = below_ef+1, tot_orbs
              v1b = hbar3b_preom_I2_test(c)
              t2 = t2_ccm_test(c,b,i,j)
              r2_preom_eqn_test(b,i,j) = r2_preom_eqn_test(b,i,j) + v1b*t2
           end DO
                      
           IF ( eom_approx > 0 ) then
              ! <b|r|ij> <-- - (1/2).P(ij).<kl|v|jc>.<bc|r|ikl>
              DO c = below_ef+1, tot_orbs
                 DO k = 1, below_ef-1
                    DO l = k+1, below_ef ! x 1/2
                       v2b = v2b_test(k,l,j,c)
                       r3 = r3_preom_test(b,c,i,k,l)
                       r2_preom_eqn_test(b,i,j) = r2_preom_eqn_test(b,i,j) - v2b*r3
                       v2b = v2b_test(k,l,i,c)
                       r3 = r3_preom_test(b,c,j,k,l)
                       r2_preom_eqn_test(b,i,j) = r2_preom_eqn_test(b,i,j) + v2b*r3
                    end DO
                 end DO
              end DO
              
              ! <b|r|ij> <-- - (1/2).<kb|v|cd>.<cd|r|ijk>
              DO c = below_ef+1, tot_orbs-1
                 DO d = c+1, tot_orbs ! x 1/2
                    DO k = 1, below_ef
                       v2b = v2b_test(k,b,c,d)
                       r3 = r3_preom_test(c,d,i,j,k)
                       r2_preom_eqn_test(b,i,j) = r2_preom_eqn_test(b,i,j) - v2b*r3
                    end DO
                 end DO
              end DO
           end IF
           
           r2_preom_eqn_test(b,j,i) = -1.d0 * r2_preom_eqn_test(b,i,j)
           
        end DO
     end DO
  end DO
  

  IF ( test == 4 ) then

     IF ( iam == 0 ) write(6,*) ' Testing PR-EOM right 1p2h...'
     DO ch = 1, channels_2b%number_confs
        IF ( number_2b(1,ch) == 0 ) cycle
        IF ( r2_preom_ind(2*ch) == 0 ) cycle
        ket_confs = number_2b(1,ch)
        DO ket = 1, ket_confs
           i   = lookup_2b_configs(1,ch)%ival2(1,ket)
           j   = lookup_2b_configs(1,ch)%ival2(2,ket)
           DO bra = 1, 2
              b   = r2_preom_ind(2*(ch-1)+bra)
              IF ( abs( r2_preom_eqn(ch)%cval(bra,ket) - r2_preom_eqn_test(b,i,j) ) > 1.e-5 ) then
                 WRITE(6,*) "r2_eqn ", b, i, j, r2_preom_eqn(ch)%cval(bra,ket), r2_preom_eqn_test(b,i,j)
                 stop
              end IF
           end DO
        end DO
     end DO
     
  end IF

end SUBROUTINE build_r2_preom_eqn_test

! <bc|r|ijk> <-- - P(i/jk).<l|x|i>.<bc|r|ljk> + P(bc).<c|x|d>.<bd|r|ijk> - P(i/jk|bc).<lc|v|jk>.<b|r|il> - P(i/jk).<bc|v|id>.<d|r|jk>
SUBROUTINE build_r3_preom_eqn_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch3,ch1,ch2, bra_confs
  INTEGER :: ket_min,ket_max, bra,ket, iind
  INTEGER :: Nx2,Ny2,Nz2,Tz2
  INTEGER :: i,j,k,b,c, l,d
  COMPLEX(dpc) :: v2b, r2, r3
  
  ! <bc|r|ijk>
  DO i = 1, below_ef-2
     DO j = i+1, below_ef-1 ! i < j
        DO k = j+1, below_ef ! j < k
           Nx2 = all_orbit%nx(i) + all_orbit%nx(j) + all_orbit%nx(k)
           Ny2 = all_orbit%ny(i) + all_orbit%ny(j) + all_orbit%ny(k)
           Nz2 = all_orbit%nz(i) + all_orbit%nz(j) + all_orbit%nz(k)
           Tz2 = all_orbit%tz(i) + all_orbit%tz(j) + all_orbit%tz(k)
           DO b = below_ef+1, tot_orbs-1
              DO c = b+1, tot_orbs ! b < c
                 IF ( Nx2 - all_orbit%nx(b) - all_orbit%nx(c) /= nx_preom ) cycle
                 IF ( Ny2 - all_orbit%ny(b) - all_orbit%ny(c) /= ny_preom ) cycle
                 IF ( Nz2 - all_orbit%nz(b) - all_orbit%nz(c) /= nz_preom ) cycle
                 IF ( Tz2 - all_orbit%tz(b) - all_orbit%tz(c) /= tz_preom ) cycle

                 ! <bc|r|ijk> <-- - P(i/jk).<l|x|i>.<bc|r|ljk> + P(bc).<c|x|d>.<bd|r|ijk>
                 r3 = r3_preom_test(b,c,i,j,k)
                 r3_preom_test(b,c,i,j,k) = 0.d0
                 r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) + hbar1b_I2_test(b,b) * r3
                 r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) + hbar1b_I2_test(c,c) * r3
                 r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) - hbar1b_I3_test(i,i) * r3
                 r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) - hbar1b_I3_test(j,j) * r3
                 r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) - hbar1b_I3_test(k,k) * r3

                 ! <bc|r|ijk> <-- - P(i/jk|bc).<lc|v|jk>.<b|r|il>
                 DO l = 1, below_ef
                    v2b = v2b_test(l,c,j,k)
                    r2 = r2_preom_test(b,i,l)
                    r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) - v2b*r2
                    v2b = v2b_test(l,c,i,k)
                    r2 = r2_preom_test(b,j,l)
                    r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) + v2b*r2
                    v2b = v2b_test(l,c,j,i)
                    r2 = r2_preom_test(b,k,l)
                    r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) + v2b*r2
                    v2b = v2b_test(l,b,j,k)
                    r2 = r2_preom_test(c,i,l)
                    r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) + v2b*r2
                    v2b = v2b_test(l,b,i,k)
                    r2 = r2_preom_test(c,j,l)
                    r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) - v2b*r2
                    v2b = v2b_test(l,b,j,i)
                    r2 = r2_preom_test(c,k,l)
                    r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) - v2b*r2
                 end DO

                 ! <bc|r|ijk> <-- - P(i/jk).<bc|v|id>.<d|r|jk>
                 DO d = below_ef+1, tot_orbs
                    v2b = v2b_test(b,c,i,d)
                    r2 = r2_preom_test(d,j,k)
                    r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) - v2b*r2
                    v2b = v2b_test(b,c,j,d)
                    r2 = r2_preom_test(d,i,k)
                    r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) + v2b*r2
                    v2b = v2b_test(b,c,k,d)
                    r2 = r2_preom_test(d,j,i)
                    r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k) + v2b*r2
                 end DO
                 
                 ! r3_preom_test(b,c,i,j,k) = r3_preom_test(b,c,i,j,k)
                 r3_preom_test(b,c,j,k,i) = r3_preom_test(b,c,i,j,k)
                 r3_preom_test(b,c,k,i,j) = r3_preom_test(b,c,i,j,k)
                 r3_preom_test(b,c,j,i,k) = -1.d0 * r3_preom_test(b,c,i,j,k)
                 r3_preom_test(b,c,i,k,j) = -1.d0 * r3_preom_test(b,c,i,j,k)
                 r3_preom_test(b,c,k,j,i) = -1.d0 * r3_preom_test(b,c,i,j,k)
                 r3_preom_test(c,b,i,j,k) = -1.d0 * r3_preom_test(b,c,i,j,k)
                 r3_preom_test(c,b,j,k,i) = -1.d0 * r3_preom_test(b,c,i,j,k)
                 r3_preom_test(c,b,k,i,j) = -1.d0 * r3_preom_test(b,c,i,j,k)
                 r3_preom_test(c,b,j,i,k) = r3_preom_test(b,c,i,j,k)
                 r3_preom_test(c,b,i,k,j) = r3_preom_test(b,c,i,j,k)
                 r3_preom_test(c,b,k,j,i) = r3_preom_test(b,c,i,j,k)

              end DO
           end DO
        end DO
     end DO
  end DO
  

  IF ( test == 4 ) then

     IF ( iam == 0 ) write(6,*) ' Testing PR-EOM right 2p3h...'
     DO ch3 = ch3_preom_min, ch3_preom_max
        DO iind = klimits_preom3(ch3,1), klimits_preom3(ch3,2)
           i         = klist_preom3(ch3)%ival2(iind,1)
           ch2       = klist_preom3(ch3)%ival2(iind,2)
           ket_min   = mapping_preom3(ch3)%ival2(iind,1)
           ket_max   = mapping_preom3(ch3)%ival2(iind,2)
           ch1       = ch1_preom3(ch3)
           bra_confs = number_2b(3,ch1)
           DO bra = 1, bra_confs
              b   = lookup_2b_configs(3,ch1)%ival2(1,bra)
              c   = lookup_2b_configs(3,ch1)%ival2(2,bra)
              DO ket = ket_min, ket_max
                 j   = lookup_2b_configs(1,ch2)%ival2(1,ket)
                 k   = lookup_2b_configs(1,ch2)%ival2(2,ket)
                 IF ( i >= j ) cycle
                 IF ( abs( r3_preom(ch3)%val1(iind)%cval(bra,ket) - r3_preom_test(b,c,i,j,k) ) > 1.e-5 ) then
                    WRITE(6,*) "r3_eqn ", b, c, i, j, k, r3_preom(ch3)%val1(iind)%cval(bra,ket), r3_preom_test(b,c,i,j,k)
                    stop
                 end IF
              end DO
           end DO
        end DO
     end DO
     
  end IF
  
end SUBROUTINE build_r3_preom_eqn_test


! <ij|x|kl> <-- <ij|v|kl> + (1/2).<ij|v|cd>.<cd|t|kl>
SUBROUTINE build_hbar2b_I4pr_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch, bra_confs, bra,ket
  INTEGER :: i,j,k,l, c,d
  COMPLEX(dpc) :: v2, t2

  IF ( iam == 0 ) write(6,*) "Testing I4_2b..."
  
  ! <ij|x|kl>
  IF ( .not. ALLOCATED(hbar2b_I4_test) ) then
     ALLOCATE( hbar2b_I4_test(1:below_ef, 1:below_ef, 1:below_ef, 1:below_ef) )
  end IF
  hbar2b_I4_test = 0.d0

  DO i = 1, below_ef-1
     DO j = i+1, below_ef ! i < j
        DO k = 1, below_ef-1
           DO l = k+1, below_ef ! k < l
              ! <ij|x|kl> <-- <ij|v|kl>
              hbar2b_I4_test(i,j,k,l) = hbar2b_I4_test(i,j,k,l) + v2b_test(i,j,k,l)
              ! <ij|x|kl> <-- +(1/2).<ij|v|cd>.<cd|t|kl>
              DO c = below_ef+1, tot_orbs-1
                 DO d = c+1, tot_orbs ! x 1/2
                    t2 = t2_ccm_test(c,d,k,l)
                    v2 = v2b_test(i,j,c,d)
                    hbar2b_I4_test(i,j,k,l) = hbar2b_I4_test(i,j,k,l) + t2*v2
                 end DO
              end DO
              hbar2b_I4_test(j,i,l,k) = hbar2b_I4_test(i,j,k,l)
              hbar2b_I4_test(j,i,k,l) = -1.d0 * hbar2b_I4_test(i,j,k,l)
              hbar2b_I4_test(i,j,l,k) = -1.d0 * hbar2b_I4_test(i,j,k,l)
           end DO
        end DO
     end DO
  end DO  

  
  IF ( test == 4 ) then
     DO ch = 1, channels_2b%number_confs
        bra_confs = number_2b(1,ch)
        IF ( bra_confs <= 0 ) cycle
        DO bra = 1, bra_confs
           i   = lookup_2b_configs(1,ch)%ival2(1,bra)
           j   = lookup_2b_configs(1,ch)%ival2(2,bra)
           DO ket = 1, bra_confs
              k   = lookup_2b_configs(1,ch)%ival2(1,ket)
              l   = lookup_2b_configs(1,ch)%ival2(2,ket)
              IF ( abs( hbar2b_I4_test(i,j,k,l) - hbar2b_I4(ch)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "hbar2b I4: ", i, j, k, l, hbar2b_I4_test(i,j,k,l), hbar2b_I4(ch)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO
  end IF
  
end SUBROUTINE build_hbar2b_I4pr_test

! <ia|x|jb> <-- <ia|v|jb> - <ik|v|cb>.<ca|t|jk>
SUBROUTINE build_hbar2b_I5pr_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch2, ket_confs, bra,ket
  INTEGER :: i,a,j,b, k,c
  COMPLEX(dpc) :: v2, t2

  IF ( iam == 0 ) write(6,*) "Testing I5_2b..."
  
  ! <ia|x|jb>
  IF ( .not. ALLOCATED(hbar2b_I5_test) ) then
     ALLOCATE( hbar2b_I5_test(1:below_ef, below_ef+1:tot_orbs, 1:below_ef, below_ef+1:tot_orbs) )
  end IF
  hbar2b_I5_test = 0.d0
  
  DO i = 1, below_ef
     DO a = below_ef+1, tot_orbs
        DO j = 1, below_ef
           DO b = below_ef+1, tot_orbs
              ! <ia|x|jb> <-- <ia|v|jb>
              hbar2b_I5_test(i,a,j,b) = hbar2b_I5_test(i,a,j,b) + v2b_test(i,a,j,b)
              ! <ia|x|jb> <-- -<ik|v|cb>.<ca|t|jk>
              DO c = below_ef+1, tot_orbs
                 IF ( c == a .or. c == b ) cycle
                 DO k = 1, below_ef
                    IF ( k == i .or. k == j ) cycle
                    t2 = t2_ccm_test(c,a,j,k)
                    v2 = v2b_test(i,k,c,b)
                    hbar2b_I5_test(i,a,j,b) = hbar2b_I5_test(i,a,j,b) - t2*v2
                 end DO
              end DO
           end DO
        end DO
     end DO
  end DO

  
  IF ( test == 4 ) then
     DO ch2 = 1, channels_2bcross%number_confs
        IF ( number_2bcross(2,ch2) == 0 ) cycle
        IF ( r2_preom_ind_cross(2*ch2) == 0 ) cycle
        ket_confs = number_2bcross(2,ch2)
        DO bra = 1, ket_confs
           i   = lookup_2bcross_configs(2,ch2)%ival2(1,bra)
           b   = lookup_2bcross_configs(2,ch2)%ival2(2,bra)
           DO ket = 1, ket_confs
              j   = lookup_2bcross_configs(2,ch2)%ival2(1,ket)
              a   = lookup_2bcross_configs(2,ch2)%ival2(2,ket)
              IF ( abs( hbar2b_I5_test(i,a,j,b) - hbar2b_I5_cross(ch2)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "hbar2b I5_cross: ", i, a, j, b, hbar2b_I5_test(i,a,j,b), hbar2b_I5_cross(ch2)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO
  end IF
  
end SUBROUTINE build_hbar2b_I5pr_test

! <ia|x|jk> <-- <ia|v|jk> + P(jk).<ca|t|lk>.<il|v|jc> + (1/2).<ia|v|cd>.<cd|t|jk>
SUBROUTINE build_hbar2b_I7pr_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: ch, bra_confs,ket_confs, bra,ket
  INTEGER :: i,a,j,k, c,d,l
  COMPLEX(dpc) :: v2, t2

  IF ( iam == 0 ) write(6,*) "Testing I7_2b..."
  
  ! <ia|x|jk>
  IF ( .not. ALLOCATED(hbar2b_I7_test) ) then
     ALLOCATE( hbar2b_I7_test(1:below_ef, below_ef+1:tot_orbs, 1:below_ef, 1:below_ef) )
  end IF
  hbar2b_I7_test = 0.d0

  DO i = 1, below_ef
     DO a = below_ef+1, tot_orbs
        DO j = 1, below_ef-1
           DO k = j+1, below_ef ! j <-> k
              ! <ia|x|jk> <-- <ia|v|jk>
              hbar2b_I7_test(i,a,j,k) = hbar2b_I7_test(i,a,j,k) + v2b_test(i,a,j,k)
              ! <ia|x|jk> <-- + (1/2).<ia|v|cd>.<cd|t|jk>
              DO c = below_ef+1, tot_orbs-1
                 DO d = c+1, tot_orbs ! x 1/2
                    t2 = t2_ccm_test(c,d,j,k)
                    v2 = v2b_test(i,a,c,d)
                    hbar2b_I7_test(i,a,j,k) = hbar2b_I7_test(i,a,j,k) + t2*v2
                 end DO
              end DO
              ! <ia|x|jk> <-- + P(jk).<ca|t|lk>.<il|v|jc>
              DO l = 1, below_ef
                 IF ( l == i ) cycle
                 DO c = below_ef+1, tot_orbs
                    IF ( c == a ) cycle
                    t2 = t2_ccm_test(c,a,l,k)
                    v2 = v2b_test(i,l,j,c)
                    hbar2b_I7_test(i,a,j,k) = hbar2b_I7_test(i,a,j,k) + t2*v2
                    t2 = t2_ccm_test(c,a,l,j)
                    v2 = v2b_test(i,l,k,c)
                    hbar2b_I7_test(i,a,j,k) = hbar2b_I7_test(i,a,j,k) - t2*v2
                 end DO
              end DO
              hbar2b_I7_test(i,a,k,j) = -1.d0 * hbar2b_I7_test(i,a,j,k)
           end DO
        end DO
     end DO
  end DO

  
  IF ( test == 4 ) then
     DO ch = 1, channels_2b%number_confs
        IF ( number_2b(1,ch) == 0 ) cycle
        IF ( r2_preom_ind(2*ch) == 0 ) cycle
        bra_confs = number_2b(2,ch)
        ket_confs = number_2b(1,ch)
        IF ( bra_confs <= 0 ) cycle
        DO bra = 1, bra_confs
           i   = lookup_2b_configs(2,ch)%ival2(1,bra)
           a   = lookup_2b_configs(2,ch)%ival2(2,bra)
           DO ket = 1, ket_confs
              j   = lookup_2b_configs(1,ch)%ival2(1,ket)
              k   = lookup_2b_configs(1,ch)%ival2(2,ket)
              IF ( abs( hbar2b_I7_test(i,a,j,k) - hbar2b_I7(ch)%cval(bra,ket) ) > 1.e-5 ) then
                 WRITE(6,*) "hbar2b I7: ", i, a, j, k, hbar2b_I7_test(i,a,j,k), hbar2b_I7(ch)%cval(bra,ket)
                 stop
              end IF
           end DO
        end DO
     end DO
  end IF
  
end SUBROUTINE build_hbar2b_I7pr_test

! <|z3b|a> <-- -(1/2).<kl|v|ad>.<d|r|kl>
SUBROUTINE build_hbar3b_preom_test
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  INTEGER :: a, k,l,d
  COMPLEX(dpc) :: v2, r2

  IF ( iam == 0 ) write(6,*) "Testing I2_3b_preom..."
  
  ! <|z3b|a>
  IF ( .not. ALLOCATED(hbar3b_preom_I2_test) ) then
     ALLOCATE( hbar3b_preom_I2_test(below_ef+1:tot_orbs) )
  end IF
  hbar3b_preom_I2_test = 0.d0

  DO a = below_ef+1, tot_orbs

     ! <|z3b|a> <-- -(1/2).<kl|v|ad>.<d|r|kl>
     DO d = below_ef+1, tot_orbs
        IF ( d == a ) cycle
        DO k = 1, below_ef-1
           DO l = k+1, below_ef ! x 1/2
              r2 = r2_preom_test(d,k,l)
              v2 = v2b_test(k,l,a,d)
              hbar3b_preom_I2_test(a) = hbar3b_preom_I2_test(a) - r2*v2
           end DO
        end DO
     end DO

  end DO
  
  
  IF ( test == 4 ) then
     DO a = below_ef+1, tot_orbs
        IF ( abs( hbar3b_preom_I2_test(a) - hbar3b_preom_I2(a) ) > 1.e-5 ) then
           WRITE(6,*) "hbar3b I2: ", a, hbar3b_preom_I2_test(a), hbar3b_preom_I2(a)
           stop
        end IF
     end DO
  end IF
  
end SUBROUTINE build_hbar3b_preom_test


! E = (1/4).<ij|v|ab><ab|v|ij>/e_abij
SUBROUTINE mbpt2_test(e2)
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(IN) :: e2
  INTEGER :: a,b,i,j
  COMPLEX(dpc) :: etest, denom
  
  IF ( iam == 0 ) write(6,*) "Testing MBPT2..."

  etest = 0.d0
  DO a = below_ef+1, tot_orbs-1
     DO b = a+1, tot_orbs ! a < b
        DO i = 1, below_ef-1
           DO j = i+1, below_ef ! i < j
              denom = fock_test(i,i) + fock_test(j,j) - fock_test(a,a) - fock_test(b,b)
              etest = etest + v2b_test(i,j,a,b)*v2b_test(a,b,i,j)/denom
           end DO
        end DO
     end DO
  end DO

  IF ( abs( etest - e2 ) > 1.e-5 ) then
     WRITE(6,*) "MBPT2: ", etest, e2
     stop
  end IF
  
end SUBROUTINE mbpt2_test

! E = (1/8).<kl|v|ab>.<ij|v|kl>.<ab|v|ij>/(e_abij.e_abkl)
SUBROUTINE mbpt3_hh_test(e3)
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(IN) :: e3
  INTEGER :: a,b,i,j, k,l
  COMPLEX(dpc) :: etest, v2, denom1, denom2
  
  IF ( iam == 0 ) write(6,*) "Testing MBPT3_hh..."

  etest = 0.d0
  DO i = 1, below_ef-1
     DO j = i+1, below_ef ! i < j
        DO k = 1, below_ef-1
           DO l = k+1, below_ef ! k < l
              v2 = v2b_test(i,j,k,l)
              DO a = below_ef+1, tot_orbs-1
                 DO b = a+1, tot_orbs ! a < b
                    denom1 = fock_test(i,i) + fock_test(j,j) - fock_test(a,a) - fock_test(b,b)
                    denom2 = fock_test(k,k) + fock_test(l,l) - fock_test(a,a) - fock_test(b,b)
                    etest = etest + v2b_test(k,l,a,b)*v2b_test(a,b,i,j)*v2/(denom1*denom2)
                 end DO
              end DO
           end DO
        end DO
     end DO
  end DO
  
  IF ( abs( etest - e3 ) > 1.e-5 ) then
     WRITE(6,*) "MBPT3_hh: ", etest, e3
     stop
  end IF
  
end SUBROUTINE mbpt3_hh_test

! E = (1/8).<ij|v|ab>.<ab|v|cd>.<cd|v|ij>/(e_abij.e_cdij)
SUBROUTINE mbpt3_pp_test(e3)
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(IN) :: e3
  INTEGER :: a,b,i,j, c,d
  COMPLEX(dpc) :: etest, v2, denom1, denom2
  
  IF ( iam == 0 ) write(6,*) "Testing MBPT3_pp..."

  etest = 0.d0
  DO a = below_ef+1, tot_orbs-1
     DO b = a+1, tot_orbs ! a < b
        DO c = below_ef+1, tot_orbs-1
           DO d = c+1, tot_orbs ! c < d
              v2 = v2b_test(a,b,c,d)
              DO i = 1, below_ef-1
                 DO j = i+1, below_ef ! i < j
                    denom1 = fock_test(i,i) + fock_test(j,j) - fock_test(a,a) - fock_test(b,b)
                    denom2 = fock_test(i,i) + fock_test(j,j) - fock_test(c,c) - fock_test(d,d)
                    etest = etest + v2b_test(i,j,a,b)*v2b_test(c,d,i,j)*v2/(denom1*denom2)
                 end DO
              end DO
           end DO
        end DO
     end DO
  end DO
  
  IF ( abs( etest - e3 ) > 1.e-5 ) then
     WRITE(6,*) "MBPT3_pp: ", etest, e3
     stop
  end IF
  
end SUBROUTINE mbpt3_pp_test

! E = -<kj|v|ac>.<ic|v|kb>.<ab|v|ij>/(e_ackj.e_abij)
SUBROUTINE mbpt3_ph_test(e3)
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(IN) :: e3
  INTEGER :: a,b,i,j, k,c
  COMPLEX(dpc) :: etest, v2, denom1, denom2
  
  IF ( iam == 0 ) write(6,*) "Testing MBPT3_ph..."

  etest = 0.d0
  DO i = 1, below_ef
     DO c = below_ef+1, tot_orbs
        DO k = 1, below_ef
           DO b = below_ef+1, tot_orbs
              v2 = v2b_test(i,c,k,b)
              DO j = 1, below_ef
                 DO a = below_ef+1, tot_orbs
                    denom1 = fock_test(i,i) + fock_test(j,j) - fock_test(a,a) - fock_test(b,b)
                    denom2 = fock_test(k,k) + fock_test(j,j) - fock_test(a,a) - fock_test(c,c)
                    etest = etest - v2b_test(k,j,a,c)*v2b_test(a,b,i,j)*v2/(denom1*denom2)
                 end DO
              end DO
           end DO
        end DO
     end DO
  end DO
  
  IF ( abs( etest - e3 ) > 1.e-5 ) then
     WRITE(6,*) "MBPT3_ph: ", etest, e3
     stop
  end IF
  
end SUBROUTINE mbpt3_ph_test

! E = +(1/2).<ab|v|ij>.<ij|v|ac>.<c|v|b>/(e_abij.e_acij)
!     -(1/2).<ab|v|ij>.<ik|v|ab>.<j|v|k>/(e_abij.e_abik)
SUBROUTINE mbpt3_1b_test(e3)
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE parallel

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(IN) :: e3
  INTEGER :: a,b,i,j, k,c
  COMPLEX(dpc) :: etest, v2, denom1, denom2
  
  IF ( iam == 0 ) write(6,*) "Testing MBPT3_1b..."
  
  etest = 0.d0
  DO a = below_ef+1, tot_orbs-1
     DO b = a+1, tot_orbs ! a < b
        DO i = 1, below_ef-1
           DO j = i+1, below_ef ! i < j
              denom1 = fock_test(i,i) + fock_test(j,j) - fock_test(a,a) - fock_test(b,b)
              v2 = v2b_test(a,b,i,j)
              DO c = below_ef+1, tot_orbs
                 IF ( c == b ) cycle
                 denom2 = fock_test(i,i) + fock_test(j,j) - fock_test(a,a) - fock_test(c,c)
                 etest = etest + v2*v2b_test(i,j,a,c)*fock_test(c,b)/(denom1*denom2)
              end DO
              DO k = 1, below_ef
                 IF ( k == j ) cycle
                 denom2 = fock_test(i,i) + fock_test(k,k) - fock_test(a,a) - fock_test(b,b)
                 etest = etest - v2*v2b_test(i,k,a,b)*fock_test(j,k)/(denom1*denom2)
              end DO
           end DO
        end DO
     end DO
  end DO
  
  IF ( abs( etest - e3 ) > 1.e-5 ) then
     WRITE(6,*) "MBPT3_1b: ", etest, e3
     stop
  end IF
  
end SUBROUTINE mbpt3_1b_test

! E = (1/36).<ijk|w|abc>.<abc|w|ijk>/e_abcijk
SUBROUTINE mbpt2_3b_test(e2)
  USE single_particle_orbits
  USE constants
  USE operator_storage
  USE configurations
  USE chiral_potentials
  USE parallel

  IMPLICIT NONE
  COMPLEX(dpc), INTENT(IN) :: e2
  INTEGER :: a,b,c,i,j,k
  COMPLEX(dpc) :: v3b, etest, denom
  
  IF ( iam == 0 ) write(6,*) "Testing MBPT2_3b..."

  etest = 0.d0
  DO a = below_ef+1, tot_orbs-2
     DO b = a+1, tot_orbs-1 ! a < b
        DO c = b+1, tot_orbs ! b < c
           DO i = 1, below_ef-2
              DO j = i+1, below_ef-1 ! i < j
                 DO k = j+1, below_ef ! j < k
                    denom = fock_test(i,i) + fock_test(j,j) + fock_test(k,k) &
                         - fock_test(a,a) - fock_test(b,b) - fock_test(c,c)
                    v3b = v3int(a,b,c,i,j,k)
                    etest = etest + conjg(v3b)*v3b/denom
                 end DO
              end DO
           end DO
        end DO
     end DO
  end DO

  IF ( abs( etest - e2 ) > 1.e-5 ) then
     WRITE(6,*) "MBPT2_3b: ", etest, e2
     stop
  end IF
  
end SUBROUTINE mbpt2_3b_test
