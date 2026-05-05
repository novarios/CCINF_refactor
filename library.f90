
SUBROUTINE sort_eigs(a,n)
  USE kind_params
  IMPLICIT NONE  
  INTEGER :: i, j, n
  COMPLEX(dpc), DIMENSION(n), INTENT(INOUT) :: a
  COMPLEX(dpc) :: temp1
  DO i = 1, n
     DO j = 1, n
        IF ( real( a(i) ) < real( a(j) ) ) THEN
           temp1 = a(i)
           a(i) = a(j) 
           a(j) = temp1
        END IF
     END DO
  END DO
END SUBROUTINE sort_eigs

! Sort routine
SUBROUTINE sort_channel( c, b, n )
  USE kind_params
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n
  INTEGER, INTENT(inout) :: c(n), b(n)
  INTEGER :: i,j, tmp1, tmp2
  DO i = 1, n
     DO j = i+1, n
        IF ( c(j) < c(i) ) THEN
           tmp1 = c(j)
           c(j) = c(i)
           c(i) = tmp1
       
           tmp2 = b(j) 
           b(j) = b(i)
           b(i) = tmp2
        END IF
     END DO
  END DO
END SUBROUTINE sort_channel

! Swaps values of 2 integers
SUBROUTINE swap_ab(a,b)
  USE kind_params
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: a, b
  INTEGER :: c
  c = a; a = b; b = c
END SUBROUTINE SWAP_AB

! Function to check triangular relations
LOGICAL FUNCTION triag(i,j,k)
  USE kind_params
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i, j, k
  triag = ((i-j-k)*(i-ABS(j-k)) > 0)
END FUNCTION triag

! Function to calculate norm of g-mat
REAL(dp) FUNCTION dij(ja,jb)
  USE kind_params
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ja, jb
  IF(ja == jb ) THEN
     dij=sqrt(2.d0)
  ELSE
     dij=1.0
  ENDIF
END FUNCTION dij

! Function to calculate phase factors (-1)**(n)
INTEGER FUNCTION iph(n)
  USE kind_params
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  iph=(-1)**n
END FUNCTION iph

FUNCTION delta(a,b)
  USE kind_params
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: a,b
  REAL(dp) :: delta
  delta=0.0
  if(a == b)delta=1.0
end FUNCTION delta

SUBROUTINE zero(n,a)
  USE kind_params
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: i
  REAL(dp), INTENT(OUT) :: a(n)
  do i=1,n
    a(i)=0.0
  end do
end subroutine zero


SUBROUTINE lapack_diag_real(h, vec, eig_r, n )
  USE kind_params
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n
  REAL(dp), DIMENSION(n,n), INTENT(in) :: h
  REAL(dp), DIMENSION(n,n), INTENT(out) :: vec
  REAL(dp), DIMENSION(n), INTENT(out) :: eig_r
  REAL(dp), DIMENSION(n) :: eig_i
  REAL(dp), DIMENSION(n,n) ::  vl
  REAL(dp), DIMENSION(20000) :: work
  DOUBLE PRECISION, DIMENSION(2*n) :: rwork
  INTEGER :: i, lda, ldvl, ldvr, info, lwork
  CHARACTER*1 :: jobvl, jobvr
  REAL(dp) :: norm

  jobvl = 'N' ;  jobvr = 'V';  lda = n
  ldvl = 1;  ldvr = n;  lwork = 20000
  eig_r = 0.; eig_i = 0.; vec = 0.
  CALL dgeev( jobvl, jobvr, n, h, lda, eig_r, eig_i, &
       vl, ldvl, vec, ldvr, work, lwork, rwork, info )
  
  ! normalization
  do i = 1, n
     norm = sum( vec(:,i)*vec(:,i) )
     vec(:,i) = vec(:,i)/sqrt(norm)
  end do
  
  call hf_sort_real(eig_r, vec, n)

end SUBROUTINE lapack_diag_real

!
! eigenvalue sort
! sort vector a(1) < a(2) < ... < a(n)
SUBROUTINE hf_sort_real(a,b, n)
  USE kind_params
  IMPLICIT NONE  
  INTEGER :: i, j, n
  REAL(dp) :: temp1
  REAL(dp), DIMENSION(n), INTENT(INOUT) :: a
  REAL(dp), DIMENSION(n,n), INTENT(INOUT) :: b
  REAL(dp), DIMENSION(n) :: temp3
  
  DO i = 1, n
     DO j = 1, n
        IF ( a(i) < a(j) ) THEN           
           temp1 = a(i)
           a(i) = a(j) 
           a(j) = temp1
           
           temp3(:) = b(:,i)
           b(:,i) = b(:,j) 
           b(:,j) = temp3(:)
        END IF
     END DO
  END DO
  
END SUBROUTINE hf_sort_real


SUBROUTINE lapack_diag_cmplx(h, cvec, ceig, n )
  USE kind_params
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n
  COMPLEX(dpc), DIMENSION(n,n), INTENT(in) :: h
  COMPLEX(dpc), DIMENSION(n,n), INTENT(out) :: cvec
  COMPLEX(dpc), DIMENSION(n), INTENT(out) :: ceig
  COMPLEX(dpc), DIMENSION(n,n) ::  vl
  COMPLEX(dpc), DIMENSION(20000) :: work
  DOUBLE PRECISION, DIMENSION(2*n) :: rwork
  INTEGER :: i, lda, ldvl, ldvr, info, lwork
  CHARACTER*1 :: jobvl, jobvr
  COMPLEX(dpc) :: norm

  jobvl = 'N' ;  jobvr = 'V';  lda = n
  ldvl = 1;  ldvr = n;  lwork = 20000
  ceig = 0.; cvec = 0.
  CALL zgeev( jobvl, jobvr, n, h, lda, ceig, &
       vl, ldvl, cvec, ldvr, work, lwork, rwork, info )
  
  ! normalization
  do i = 1, n
     norm = sum( cvec(:,i)*cvec(:,i) )
     cvec(:,i) = cvec(:,i)/sqrt(norm)
  end do
  
  call hf_sort_cmplx(ceig, cvec, n)

end SUBROUTINE lapack_diag_cmplx

!
! eigenvalue sort
! sort vector re[a(1)] < re[a(2)] < ... < re[a(n)]
SUBROUTINE hf_sort_cmplx(a,b, n)
  USE kind_params
  IMPLICIT NONE  
  INTEGER :: i, j, n
  COMPLEX(dpc) :: temp1
  COMPLEX(dpc), DIMENSION(n), INTENT(INOUT) :: a
  COMPLEX(dpc), DIMENSION(n,n), INTENT(INOUT) :: b
  COMPLEX(dpc), DIMENSION(n) :: temp3
  
  DO i = 1, n
     DO j = 1, n
        IF ( real( a(i) ) < real( a(j) ) ) THEN
           temp1 = a(i)
           a(i) = a(j) 
           a(j) = temp1
           
           temp3(:) = b(:,i)
           b(:,i) = b(:,j) 
           b(:,j) = temp3(:)
        END IF
     END DO
  END DO
  
END SUBROUTINE hf_sort_cmplx
