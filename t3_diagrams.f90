
MODULE t3_diagrams
  
CONTAINS

  !
  ! T3 diagram 1: <abc|t|ijk> <-- -<cd|t|ij>.<ab|v|kd>
  ! Loops over particle index d, contracts T2 with V2(pphp).
  !
  ! cval is passed explicitly (not accessed via t3_ccm module variable)
  ! so that each OpenMP thread gets its own array descriptor on the
  ! stack. bra_lo preserves the original lower bound of the pointer-
  ! remapped array, since assumed-shape arguments rebase to 1.
  !
  SUBROUTINE t3_diag1(ch3, ch1,ch2, cind1,kind1, c1,k, bra,bra1,bra_min, phase1, cval)
    USE single_particle_orbits
    USE configurations
    USE constants
    USE operator_storage
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ch3, ch1,ch2, cind1,kind1, c1,k, bra,bra1,bra_min, phase1
    COMPLEX(dpc), INTENT(INOUT) :: cval(bra_min:,:)
    INTEGER :: bra0,ket0, id,d, ket,ket1,ket_confs, phase
    COMPLEX(dpc) :: factor
    
    ket_confs = number_2b_t3(ch3)%ival2(2,ch2)  
    IF ( ket_confs <= 0 ) return
    
    DO id = 1, t3_hp_inv(ch1, k)%n
       d = t3_hp_inv(ch1, k)%idx(id)
       IF ( d == c1 ) cycle
       IF ( ch2 /= pp_channel_2b%ival2(min(d,c1), max(d,c1)) ) cycle
       ket0 = hp_config_2b%ival2(k,d)
       bra0 = pp_config_2b%ival2(d,c1)
       factor = phase1 * sign(1, c1 - d) * v2b_pphp(ch1)%cval(bra1,ket0)
       DO ket = 1, ket_confs
          ket1 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
          cval(bra,ket) = cval(bra,ket) + factor * t2_ccm(ch2)%cval(bra0,ket1)
       end DO
    end DO
    
  end SUBROUTINE t3_diag1
  
  !
  ! T3 diagram 2: <abc|t|ijk> <-- +<ab|t|lk>.<lc|v|ij>
  ! Loops over hole index l, contracts T2 with V2(hphh).
  !
  SUBROUTINE t3_diag2(ch3, ch1,ch2, cind1,kind1, c1,k, bra,bra1,bra_min, phase1, cval)
    USE single_particle_orbits
    USE configurations
    USE constants
    USE operator_storage
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ch3, ch1,ch2, cind1,kind1, c1,k, bra,bra1,bra_min, phase1
    COMPLEX(dpc), INTENT(INOUT) :: cval(bra_min:,:)
    INTEGER :: bra0,ket0, il,l, ket,ket1,ket_confs
    COMPLEX(dpc) :: factor
    
    ket_confs = number_2b_t3(ch3)%ival2(2,ch2)
    IF ( ket_confs <= 0 ) return

    DO il = 1, t3_hh_inv(ch1, k)%n
       l = t3_hh_inv(ch1, k)%idx(il)
       IF ( l == k ) cycle
       IF ( ch2 /= hp_channel_2b%ival2(l,c1) ) cycle
       bra0 = hp_config_2b%ival2(l,c1)
       ket0 = hh_config_2b%ival2(l,k)
       factor = phase1 * sign(1, k - l) * t2_ccm(ch1)%cval(bra1,ket0)
       DO ket = 1, ket_confs
          ket1 = hh_config_t3(ch3)%ival1(ch2)%ival1(ket)
          cval(bra,ket) = cval(bra,ket) + factor * v2b_hphh(ch2)%cval(bra0,ket1)
       end DO
    end DO
    
  end SUBROUTINE t3_diag2  
  
end MODULE t3_diagrams
