MODULE hdf5_wrapper
  USE hdf5
  USE h5lt
  USE h5l
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: read_double_matrix
  PUBLIC :: write_double_matrix
  PUBLIC :: read_integer_vector
  PUBLIC :: write_integer_vector
  PUBLIC :: read_double_vector
  PUBLIC :: write_double_vector
  ! PUBLIC :: write_cc
  ! PUBLIC :: read_cc

CONTAINS
  
  FUNCTION write_double_matrix(hid, matrix, dset_name) RESULT(error)
    INTEGER(HID_T) :: hid, did
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:,:) :: matrix
    CHARACTER(LEN=*) :: dset_name
    INTEGER(HSIZE_T) :: dims(2)
    INTEGER :: error, rank
    
    rank = 2
    dims = SHAPE(matrix)
    IF (h5ltfind_dataset_f(hid, dset_name) /= 1) then
       CALL h5ltmake_dataset_double_f(hid, dset_name, rank, dims, matrix, error)
    ELSE
       call h5dopen_f(hid, dset_name, did, error)
       call h5dwrite_f(did, H5T_NATIVE_DOUBLE, matrix, dims, error)
       call h5dclose_f(did, error)
    END IF
  END FUNCTION write_double_matrix

  
  FUNCTION read_double_matrix(hid, matrix, dset_name) RESULT(error)
    INTEGER(HID_T) :: hid
    DOUBLE PRECISION, POINTER, INTENT(INOUT), DIMENSION(:,:) :: matrix
    CHARACTER(LEN=*) :: dset_name
    INTEGER(HSIZE_T) :: dims(2)
    INTEGER(SIZE_T) :: type_size
    INTEGER :: type_class, error
    
    IF (ASSOCIATED(matrix)) DEALLOCATE(matrix)
    
    IF (h5ltfind_dataset_f(hid, dset_name) /= 1) RETURN
    CALL h5ltget_dataset_info_f(hid, dset_name, dims, type_class, type_size, error)
    
    ALLOCATE(matrix(dims(1), dims(2)))
    CALL h5ltread_dataset_double_f(hid, dset_name, matrix, dims, error)
  END FUNCTION read_double_matrix


  FUNCTION write_integer_vector(hid, matrix, dset_name) RESULT(error)
    INTEGER(HID_T) :: hid, did
    INTEGER, INTENT(IN), DIMENSION(:) :: matrix
    CHARACTER(LEN=*) :: dset_name
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER :: error, rank
    
    rank = 1
    dims = SHAPE(matrix)
    IF (h5ltfind_dataset_f(hid, dset_name) /= 1) then
       CALL h5ltmake_dataset_int_f(hid, dset_name, rank, dims, matrix, error)
    ELSE
       call h5dopen_f(hid, dset_name, did, error)
       call h5dwrite_f(did, H5T_NATIVE_INTEGER, matrix, dims, error)
       call h5dclose_f(did, error)
    END IF
  END FUNCTION write_integer_vector
  
  
  FUNCTION read_integer_vector(hid, matrix, dset_name) RESULT(error)
    INTEGER(HID_T) :: hid
    INTEGER, POINTER, INTENT(INOUT), DIMENSION(:) :: matrix
    CHARACTER(LEN=*) :: dset_name
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER(SIZE_T) :: type_size
    INTEGER :: type_class, error
    
    IF (ASSOCIATED(matrix)) DEALLOCATE(matrix)
    
    IF (h5ltfind_dataset_f(hid, dset_name) /= 1) RETURN
    CALL h5ltget_dataset_info_f(hid, dset_name, dims, type_class, type_size, error)
    
    ALLOCATE(matrix(dims(1)))
    CALL h5ltread_dataset_int_f(hid, dset_name, matrix, dims, error)
  END FUNCTION read_integer_vector


  FUNCTION write_double_vector(hid, matrix, dset_name) RESULT(error)
    INTEGER(HID_T) :: hid, did
    DOUBLE PRECISION, INTENT(IN), DIMENSION(:) :: matrix
    CHARACTER(LEN=*) :: dset_name
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER :: error, rank
    
    rank = 1
    dims = SHAPE(matrix)
    IF (h5ltfind_dataset_f(hid, dset_name) /= 1) then
       CALL h5ltmake_dataset_double_f(hid, dset_name, rank, dims, matrix, error)
    ELSE
       call h5dopen_f(hid, dset_name, did, error)
       call h5dwrite_f(did, H5T_NATIVE_DOUBLE, matrix, dims, error)
       call h5dclose_f(did, error)
    END IF
  END FUNCTION write_double_vector

  
  FUNCTION read_double_vector(hid, matrix, dset_name) RESULT(error)
    INTEGER(HID_T) :: hid
    DOUBLE PRECISION, POINTER, INTENT(INOUT), DIMENSION(:) :: matrix
    CHARACTER(LEN=*) :: dset_name
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER(SIZE_T) :: type_size
    INTEGER :: type_class, error
    
    IF (ASSOCIATED(matrix)) DEALLOCATE(matrix)
    
    IF (h5ltfind_dataset_f(hid, dset_name) /= 1) RETURN
    CALL h5ltget_dataset_info_f(hid, dset_name, dims, type_class, type_size, error)
    
    ALLOCATE(matrix(dims(1)))
    CALL h5ltread_dataset_double_f(hid, dset_name, matrix, dims, error)
  END FUNCTION read_double_vector
  
  
  ! FUNCTION write_cc(filename) RESULT(error)
  !   use constants
  !   use operator_storage
  !   use configurations
  !   use parallel
    
  !   implicit none
  !   character(len=*), intent(in) :: filename
  !   character(len=13), parameter :: dsetname0 = "CC_parameters"
  !   character(len=16), parameter :: dsetname1 = "CC_energies_real"
  !   character(len=16), parameter :: dsetname2 = "CC_energies_imag"
  !   character(len=19), parameter :: dsetname30 = "T2_Amplitudes_real_"
  !   character(len=19), parameter :: dsetname40 = "T2_Amplitudes_imag_"
  !   character(len=100) :: dsetname3, dsetname4
    
  !   integer(hid_t) :: file_id
  !   integer(hid_t) :: dset_id0, dset_id1, dset_id2, dset_id3, dset_id4
  !   integer(hid_t) :: filespace0, filespace1, filespace2
  !   integer(hid_t) :: memspace0, memspace1, memspace2, memspace3, memspace4
  !   integer(hid_t) :: plist_id
  !   integer(hsize_t), dimension(1) :: p_dims
  !   integer(hsize_t), dimension(1) :: e_dims
  !   integer(hsize_t), dimension(2) :: vec_dims
  !   integer :: comm, info, error
  !   integer :: rank1 = 1
  !   integer :: rank2 = 2
  !   integer :: params(6)
  !   real(dp) :: energies(4)
  !   integer(hsize_t) :: pdim = 6 ! Nmax, #N, #P, density*1000, CC_approx, tnf_approx
  !   integer(hsize_t) :: edim = 4 ! E_CC, E0, Ecorr2, Ecorr3
    
  !   integer :: ch, bra_confs,ket_confs
    
  !   comm = MPI_COMM_WORLD
  !   info = MPI_INFO_NULL

  !   p_dims = (/ pdim /)
  !   e_dims = (/ edim /)

  !   call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  !   call h5pset_fapl_mpio_f(plist_id, comm, info, error)
  !   call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
  !   call h5pclose_f(plist_id, error)

  !   ! Write CC Parameters
  !   params = 0
  !   params(1) = Nmax
  !   params(2) = Nocc
  !   params(3) = Pocc
  !   params(4) = int(1000.d0*rho)
  !   params(5) = cc_approx
  !   params(6) = tnf_approx
  !   call h5screate_simple_f(rank1, p_dims, filespace0, error)
  !   call h5dcreate_f(file_id, trim(dsetname0), H5T_NATIVE_INTEGER, filespace0, dset_id0, error)
  !   call h5sclose_f(filespace0, error)
  !   if ( iam == 0 ) then
  !      call h5screate_simple_f(rank1, p_dims, memspace0, error)
  !      call h5dget_space_f(dset_id0, filespace0, error)
  !      call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
  !      call h5dwrite_f(dset_id0, H5T_NATIVE_INTEGER, params(1:pdim), p_dims, error, &
  !           file_space_id = filespace0, mem_space_id = memspace0, xfer_prp = plist_id)
  !      call h5pclose_f(plist_id, error)
  !      call h5sclose_f(filespace0, error)
  !      call h5sclose_f(memspace0, error)
  !   end if
  !   call h5dclose_f(dset_id0, error)
    
  !   ! Write CC Energies
  !   call h5screate_simple_f(rank1, e_dims, filespace1, error)
  !   call h5dcreate_f(file_id, trim(dsetname1), H5T_NATIVE_DOUBLE, filespace1, dset_id1, error)
  !   call h5dcreate_f(file_id, trim(dsetname2), H5T_NATIVE_DOUBLE, filespace1, dset_id2, error)
  !   call h5sclose_f(filespace1, error)
  !   energies = 0.d0
  !   energies(1) = real(eccdt)
  !   energies(2) = real(e0)
  !   energies(3) = real(ecorr2)
  !   if ( cc_approx > 1 .and. tnf_approx > 1 ) then
  !      energies(4) = real(ecorr3)
  !   end if
  !   if ( iam == 0 ) then
  !      call h5screate_simple_f(rank1, e_dims, memspace1, error)
  !      call h5dget_space_f(dset_id1, filespace1, error)
  !      call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
  !      call h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, energies(1:edim), e_dims, error, &
  !           file_space_id = filespace1, mem_space_id = memspace1, xfer_prp = plist_id)
  !      call h5pclose_f(plist_id, error)
  !      call h5sclose_f(filespace1, error)
  !      call h5sclose_f(memspace1, error)
  !   end if
  !   energies = 0.d0
  !   energies(1) = aimag(eccdt)
  !   energies(2) = aimag(e0)
  !   energies(3) = aimag(ecorr2)
  !   if ( cc_approx > 1 .and. tnf_approx > 1 ) then
  !      energies(4) = aimag(ecorr3)
  !   end if
  !   if ( iam == 0 ) then
  !      call h5screate_simple_f(rank1, e_dims, memspace2, error)
  !      call h5dget_space_f(dset_id2, filespace1, error)
  !      call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
  !      call h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, energies(1:edim), e_dims, error, &
  !           file_space_id = filespace1, mem_space_id = memspace2, xfer_prp = plist_id)
  !      call h5pclose_f(plist_id, error)
  !      call h5sclose_f(filespace1, error)
  !      call h5sclose_f(memspace2, error)
  !   end if
  !   call h5dclose_f(dset_id1, error)
  !   call h5dclose_f(dset_id2, error)

  !   ! Write T2 Amplitudes
  !   do ch = 1, channels_2b%number_confs
  !      if ( number_2b(3,ch) == 0 ) cycle
  !      if ( number_2b(1,ch) == 0 ) cycle
  !      bra_confs = number_2b(3,ch)
  !      ket_confs = number_2b(1,ch)
  !      vec_dims = (/ int(bra_confs,8), int(ket_confs,8) /)
  !      write(dsetname3,'(a,i4.4)') dsetname30, ch
  !      write(dsetname4,'(a,i4.4)') dsetname40, ch
  !      call h5screate_simple_f(rank2, vec_dims, filespace2, error)    
  !      call h5dcreate_f(file_id, trim(dsetname3), H5T_NATIVE_DOUBLE, filespace2, dset_id3, error)
  !      call h5dcreate_f(file_id, trim(dsetname4), H5T_NATIVE_DOUBLE, filespace2, dset_id4, error)
  !      call h5sclose_f(filespace2, error)
  !      if ( iam == 0 ) then
  !         call h5screate_simple_f(rank2, vec_dims, memspace3, error)
  !         call h5dget_space_f(dset_id3, filespace2, error)
  !         call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
  !         call h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, real(t2_ccm(ch)%cval(1:bra_confs,1:ket_confs)), vec_dims, error, &
  !              file_space_id = filespace2, mem_space_id = memspace3, xfer_prp = plist_id)
  !         call h5pclose_f(plist_id, error)
  !         call h5sclose_f(filespace2, error)
  !         call h5sclose_f(memspace3, error)
          
  !         call h5screate_simple_f(rank2, vec_dims, memspace4, error)
  !         call h5dget_space_f(dset_id4, filespace2, error)
  !         call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
  !         call h5dwrite_f(dset_id4, H5T_NATIVE_DOUBLE, aimag(t2_ccm(ch)%cval(1:bra_confs,1:ket_confs)), vec_dims, error, &
  !              file_space_id = filespace2, mem_space_id = memspace4, xfer_prp = plist_id)
  !         call h5pclose_f(plist_id, error)
  !         call h5sclose_f(filespace2, error)
  !         call h5sclose_f(memspace4, error)
  !      end if
  !      call h5dclose_f(dset_id3, error)
  !      call h5dclose_f(dset_id4, error)
  !   end do
    
  !   call h5fclose_f(file_id, error)    
  !   if ( iam == 0 ) write(6,*) '...Done writing CC amplitudes!'
    
  ! END FUNCTION write_cc


  ! FUNCTION read_cc(filename) RESULT(error)
  !   use constants
  !   use operator_storage
  !   use configurations
  !   use parallel
    
  !   implicit none
  !   character(len=*), intent(in) :: filename
  !   character(len=13), parameter :: dsetname0 = "CC_parameters"
  !   character(len=16), parameter :: dsetname1 = "CC_energies_real"
  !   character(len=16), parameter :: dsetname2 = "CC_energies_imag"
  !   character(len=19), parameter :: dsetname30 = "T2_Amplitudes_real_"
  !   character(len=19), parameter :: dsetname40 = "T2_Amplitudes_imag_"
  !   character(len=100) :: dsetname3, dsetname4
    
  !   integer(hid_t) :: file_id
  !   integer(hid_t) :: dset_id0, dset_id1, dset_id2, dset_id3, dset_id4
  !   integer(hid_t) :: filespace0, filespace1, filespace2
  !   integer(hid_t) :: memspace0, memspace1, memspace2, memspace3, memspace4
  !   integer(hid_t) :: plist_id
  !   integer(hsize_t), dimension(1) :: p_dims
  !   integer(hsize_t), dimension(1) :: e_dims
  !   integer(hsize_t), dimension(2) :: vec_dims
  !   integer :: comm, info, error
  !   integer :: rank1 = 1
  !   integer :: rank2 = 2
  !   integer :: params(6)
  !   real(dp) :: energies_r(4), energies_i(4)
  !   integer(hsize_t) :: pdim = 6 ! #Nocc, #Pocc, density*1000, Nmax, CC_approx, tnf_approx
  !   integer(hsize_t) :: edim = 4 ! E_CC, E0, Ecorr2, Ecorr3

  !   integer :: ch, bra,ket,bra_confs,ket_confs
  !   real(dp), allocatable :: mat_r(:,:), mat_i(:,:)
    
  !   comm = MPI_COMM_WORLD
  !   info = MPI_INFO_NULL

  !   p_dims = (/ pdim /)
  !   e_dims = (/ edim /)
    
  !   call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  !   call h5pset_fapl_mpio_f(plist_id, comm, info, error)
  !   call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
  !   IF (error .ne. 0) then
  !      IF ( iam == 0 ) write(*,*)'ERROR: HDF5 File Not Found!'
  !      stop
  !   end IF
  !   call h5pclose_f(plist_id, error)
    
  !   ! Read CC Parameters
  !   call h5screate_simple_f(rank1, p_dims, memspace0, error)
  !   call h5dopen_f(file_id, trim(dsetname0), dset_id0, error)
  !   call h5dget_space_f(dset_id0, filespace0, error)
  !   call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
  !   call h5dread_f(dset_id0, H5T_NATIVE_INTEGER, params(1:pdim), p_dims, error, &
  !        file_space_id = filespace0, mem_space_id = memspace0, xfer_prp = plist_id)
  !   call h5pclose_f(plist_id, error)
  !   call h5sclose_f(filespace0, error)
  !   call h5sclose_f(memspace0, error)
  !   call h5dclose_f(dset_id0, error)
  !   if ( iam == 0 ) write(6,*) params(1), params(2), params(3), params(4), params(5), params(6)
  !   if ( params(1) /= Nmax .or. params(2) /= Nocc .or. params(3) /= Pocc .or. params(4) /= int(1000.d0*rho) &
  !        .or. params(5) /= cc_approx .or. params(6) /= tnf_approx ) then
  !      if ( iam == 0 ) write(6,*)
  !      if ( iam == 0 ) write(6,*) ' Mismatch between CC parameters from file and input'
  !      if ( iam == 0 ) write(6,*)
  !      if ( iam == 0 ) write(6,*) ' CC parameters from file: '
  !      if ( iam == 0 ) write(6,'(A6,I4,3x,A6,I4,3x,A6,I4,3x,A5,f5.2,3x,A11,I4,3x,A12,I4)') &
  !           'Nmax =',params(1),'Nocc =',params(2),'Pocc =',params(3),'rho =',real(params(4))/1000.d0,&
  !           'cc_approx =',params(5),'tnf_approx =',params(6)
  !      if ( iam == 0 ) write(6,*) ' CC parameters from input: '
  !      if ( iam == 0 ) write(6,'(A6,I4,3x,A6,I4,3x,A6,I4,3x,A5,f5.2,3x,A11,I4,3x,A12,I4)') &
  !           'Nmax =',Nmax,'Nocc =',Nocc,'Pocc =',Pocc,'rho =',rho,&
  !           'cc_approx =',cc_approx,'tnf_approx =',tnf_approx
  !      if ( iam == 0 ) write(6,*)
  !      call h5fclose_f(file_id, error)
  !      stop
  !   end if
        
  !   ! Read CC Energies
  !   call h5screate_simple_f(rank1, e_dims, memspace1, error)
  !   call h5dopen_f(file_id, trim(dsetname1), dset_id1, error)
  !   call h5dget_space_f(dset_id1, filespace1, error)
  !   call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
  !   call h5dread_f(dset_id1, H5T_NATIVE_DOUBLE, energies_r(1:edim), e_dims, error, &
  !        file_space_id = filespace1, mem_space_id = memspace1, xfer_prp = plist_id)
  !   call h5pclose_f(plist_id, error)
  !   call h5sclose_f(filespace1, error)
  !   call h5dclose_f(dset_id1, error)
  !   call h5sclose_f(memspace1, error)

  !   call h5screate_simple_f(rank1, e_dims, memspace2, error)
  !   call h5dopen_f(file_id, trim(dsetname2), dset_id2, error)
  !   call h5dget_space_f(dset_id2, filespace1, error)
  !   call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
  !   call h5dread_f(dset_id2, H5T_NATIVE_DOUBLE, energies_i(1:edim), e_dims, error, &
  !        file_space_id = filespace1, mem_space_id = memspace2, xfer_prp = plist_id)
  !   call h5pclose_f(plist_id, error)
  !   call h5sclose_f(filespace1, error)
  !   call h5dclose_f(dset_id2, error)
  !   call h5sclose_f(memspace2, error)
  !   if ( iam == 0 ) then
  !      write(6,*)
  !      write(6,*) ' ECC = ', energies_r(1), energies_i(1)
  !      write(6,*) ' E0  = ', energies_r(2), energies_i(2)
  !      write(6,*) ' E2  = ', energies_r(3), energies_i(3)
  !      write(6,*) ' E3  = ', energies_r(4), energies_i(4)
  !      write(6,*)
  !      e0 = dcmplx(energies_r(2), energies_i(2))
  !      ecorr2 = dcmplx(energies_r(3), energies_i(3))
  !      ecorr3 = dcmplx(energies_r(4), energies_i(4))
  !      ecorr = ecorr2 + ecorr3
  !      eccdt = e0 + ecorr
  !   end if    

  !   ! Read T2 Amplitudes
  !   do ch = 1, channels_2b%number_confs
  !      if ( number_2b(3,ch) == 0 ) cycle
  !      if ( number_2b(1,ch) == 0 ) cycle
  !      bra_confs = number_2b(3,ch)
  !      ket_confs = number_2b(1,ch)
  !      allocate( mat_r(bra_confs,ket_confs) )
  !      allocate( mat_i(bra_confs,ket_confs) )
  !      mat_r = 0.d0
  !      mat_i = 0.d0
  !      vec_dims = (/ int(bra_confs,8), int(ket_confs,8) /)
  !      write(dsetname3,'(a,i4.4)') dsetname30, ch
  !      write(dsetname4,'(a,i4.4)') dsetname40, ch

  !      call h5screate_simple_f(rank2, vec_dims, memspace3, error)
  !      call h5dopen_f(file_id, trim(dsetname3), dset_id3, error)
  !      call h5dget_space_f(dset_id3, filespace2, error)
  !      call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
  !      call h5dread_f(dset_id3, H5T_NATIVE_DOUBLE, mat_r(1:bra_confs,1:ket_confs), vec_dims, error, &
  !           file_space_id = filespace2, mem_space_id = memspace3, xfer_prp = plist_id)
  !      call h5pclose_f(plist_id, error)
  !      call h5sclose_f(filespace2, error)
  !      call h5dclose_f(dset_id3, error)
  !      call h5sclose_f(memspace3, error)
       
  !      call h5screate_simple_f(rank2, vec_dims, memspace4, error)
  !      call h5dopen_f(file_id, trim(dsetname4), dset_id4, error)
  !      call h5dget_space_f(dset_id4, filespace2, error)
  !      call h5pcreate_f(h5p_dataset_xfer_f, plist_id, error)
  !      call h5dread_f(dset_id4, H5T_NATIVE_DOUBLE, mat_i(1:bra_confs,1:ket_confs), vec_dims, error, &
  !           file_space_id = filespace2, mem_space_id = memspace4, xfer_prp = plist_id)
  !      call h5pclose_f(plist_id, error)
  !      call h5sclose_f(filespace2, error)
  !      call h5dclose_f(dset_id4, error)
  !      call h5sclose_f(memspace4, error)

  !      do bra = 1, bra_confs
  !         do ket = 1, ket_confs
  !            t2_ccm(ch)%cval(bra,ket) = dcmplx(mat_r(bra,ket), mat_i(bra,ket))
  !         end do
  !      end do
  !      deallocate( mat_r )
  !      deallocate( mat_i )
  !   end do
    
  !   call h5fclose_f(file_id, error)    
  !   if ( iam == 0 ) write(6,*) '...Done reading CC amplitudes!'
    
  ! END FUNCTION read_cc    
  
END MODULE hdf5_wrapper


