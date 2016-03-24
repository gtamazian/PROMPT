#include "fintrf.h"

      ! Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      implicit none

      mwPointer     :: plhs(*), prhs(*)
      integer       :: nlhs, nrhs

      ! Function declarations
      mwPointer     :: mxCreateDoubleMatrix, mxCreateDoubleScalar
      mwPointer     :: mxGetPr, mxGetField, mxGetCell
      mwPointer     :: mxGetM, mxGetN
      integer       :: mxIsStruct, mxIsNumeric

      ! Pointers to input data
      mwPointer      :: model_ptr, index_ptr, angle_ptr
      mwPointer      :: masses_ptr, r_ptr, alpha_ptr, psi_ptr, 
     +  str_crd_ptr, end_crd_ptr, rot_mat_ptr
      
      ! Pointers to output data
      mwPointer      :: fval_ptr, grad_ptr

      ! Arguments for computational routines
      mwPointer                                   :: atom_num, conf_num
      mwPointer                                   :: angle_num
      real(kind=8), allocatable, dimension(:)     :: masses
      real(kind=8), allocatable, dimension(:,:)   :: r, alpha, psi
      real(kind=8), allocatable, dimension(:,:)   :: str_crd, end_crd
      real(kind=8), allocatable, dimension(:,:,:) :: rot_mat
      real(kind=8), allocatable, dimension(:)     :: angle_values
      integer(kind=4), allocatable, dimension(:)  :: angle_indices
      
      ! Check for the proper number of function arguments
      if (nrhs .ne. 3) then
        call mexErrMsgIdAndTxt('MATLAB:trmobjfunc:nInput',
     +                         'Three input arguments required.')
      endif

      ! Check input types
      if (mxIsStruct(prhs(1)) .eq. 0) then
        call mexErrMsgIdAndTxt('MATLAB:trmobjfunc:NonStructure',
     +                         'The first argument must be a structure
     +                          representing a transformation model.')
      endif

      if (mxIsNumeric(prhs(2)) .eq. 0) then
        call mexErrMsgIdAndTxt('MATLAB:trmobjfunc:NonNumeric',
     +                         'The second argument must be a numeric 
     +                          array.')
      endif

      if (mxIsNumeric(prhs(3)) .eq. 0) then
        call mexErrMsgIdAndTxt('MATLAB:trmobjfunc:NonNumeric',
     +                         'The third argument must be a numeric 
     +                          array.')
      endif

      ! Get pointers for input data
      model_ptr   = prhs(1)
      index_ptr   = prhs(2)
      angle_ptr   = prhs(3)

      masses_ptr  = mxGetPr(mxGetField(model_ptr, 1, 'm'))
      r_ptr       = mxGetField(model_ptr, 1, 'r')
      alpha_ptr   = mxGetField(model_ptr, 1, 'alpha')
      psi_ptr     = mxGetField(model_ptr, 1, 'psi')
      str_crd_ptr = mxGetField(model_ptr, 1, 'StartCoords')
      end_crd_ptr = mxGetField(model_ptr, 1, 'FinishCoords')
      rot_mat_ptr = mxGetField(model_ptr, 1, 'U')

      ! Determine the model parameters: the number of atoms and the
      ! number of configurations
      atom_num  = mxGetM(r_ptr) + 1
      conf_num  = mxGetN(r_ptr)
      angle_num = mxGetM(angle_ptr) * mxGetN(angle_ptr)

      ! Create result arrays
      plhs(1) = mxCreateDoubleMatrix(1, 1, 0.0)
      plhs(2) = mxCreateDoubleMatrix(angle_num, 1, 0.0)

      ! Get pointers to output data
      fval_ptr = mxGetPr(plhs(1))
      grad_ptr = mxGetPr(plhs(2))

      ! Allocate Fortran arrays
      allocate(masses(atom_num))
      allocate(r(atom_num - 1, conf_num))
      allocate(alpha(atom_num - 2, conf_num))
      allocate(psi(atom_num - 3, conf_num))
      allocate(str_crd(atom_num, 3))
      allocate(end_crd(atom_num, 3))
      allocate(rot_mat(conf_num, 3, 3))
      allocate(angle_values(angle_num))
      allocate(angle_indices(angle_num))

      ! Copy the input data to Fortran arrays
      call mxCopyPtrToReal8(masses_ptr, masses, atom_num)
      call mxCopyPtrToReal8(r_ptr, r, (atom_num - 1) * conf_num)
      call mxCopyPtrToReal8(alpha_ptr, alpha, (atom_num - 2) * 
     +  conf_num)
      call mxCopyPtrToReal8(psi_ptr, psi, (atom_num - 3) * conf_num)
      call mxCopyPtrToReal8(str_crd_ptr, str_crd, atom_num * 3)
      call mxCopyPtrToReal8(end_crd_ptr, end_crd, atom_num * 3)
      call mxCopyPtrToReal8(rot_mat_ptr, rot_mat, conf_num * 3 * 3)
      call mxCopyPtrToReal8(angle_ptr, angle_values, angle_num)
      call mxCopyPtrToInteger4(index_ptr, angle_indices, angle_num)

      ! Launch the computational routines
      call mexPrint("TODO")

      ! Deallocate Fortran arrays
      deallocate(masses, r, alpha, psi, str_crd, end_crd, rot_mat,
     +  angle_values, angle_indices)

      return
      end subroutine mexFunction

      ! Procedure implementing the cross product of a pair of vectors
      subroutine cross3d(a, b, c)
        real(kind=8), intent(in),  dimension(3) :: a, b
        real(kind=8), intent(out), dimension(3) :: c

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)

      return
      end subroutine cross3d

      ! Procedure to restore Cartesian coordinates from internal ones
      ! for a single configuration
      subroutine restoreCoords(r, alpha, psi, atom_num, coords)
        real(kind=8), intent(in),  allocatable :: r(:)
        real(kind=8), intent(in),  allocatable :: alpha(:)
        real(kind=8), intent(in),  allocatable :: psi(:)
        integer,      intent(in)               :: atom_num
        real(kind=8), intent(out), allocatable :: coords(:,:)

        real(kind=8), dimension(3)             :: bc, n
        real(kind=8), dimension(3, 3)          :: M

        allocate(coords(atom_num, 3))

        coords(1, :) = [0.0d0, 0.0d0, 0.0d0]
        coords(2, :) = [r(1), 0.0d0, 0.0d0]
        coords(3, :) = coords(2, :) + [r(2) * dcos(alpha(1)), 
     +    r(2) * dsin(alpha(1)), 0.0d0]

        coords(4:, 1) = r(3:) * dcos(alpha(2:))
        coords(4:, 2) = r(3:) * dsin(alpha(2:)) * dcos(psi)
        coords(4:, 3) = r(3:) * dsin(alpha(2:)) * dsin(psi)

        do i = 4, atom_num
          ! calculate the bc vector
          bc = coords(i-1, :) - coords(i-2, :)
          bc = bc / norm2(bc)
          ! calculate the n vector
          call cross3d(coords(i-2, :) - coords(i-3, :), bc, n)
          n = n / norm2(n)
          ! calculate the M matrix
          M(:, 1) = bc
          call cross3d(n, bc, M(:, 2))
          M(:, 3) = n
          ! get the point coordinates
          coords(i, :) = transpose(matmul(M, transpose(coords(i, :))) +
     +      coords(i-1, :)
        end do
      
      return
      end subroutine restoreCoords

      ! Procedure to restore Cartesian coordinates for all
      ! configurations of a transformation.
      subroutine trRestoreCoords(r, alpha, psi, atom_num, conf_num,
     +  rot_mat, start_coords, coords)
        real(kind=8), intent(in),  allocatable :: r(:,:)
        real(kind=8), intent(in),  allocatable :: alpha(:,:)
        real(kind=8), intent(in),  allocatable :: psi(:,:)
        integer,      intent(in)               :: atom_num
        integer,      intent(in)               :: conf_num
        real(kind=8), intent(in),  allocatable :: rot_mat(:,:,:)
        real(kind=8), intent(in),  allocatable :: start_coords(:,:) 
        real(kind=8), intent(out), allocatable, target :: coords(:,:,:)

        integer                    :: i, j
        real(kind=8), dimension(3) :: curr_trans, first_trans
        real(kind=8), pointer      :: curr_conf(:,:)

        allocate(coords(conf_num, atom_num, 3))

        coords(1, :, :) = start_coords
        first_trans = sum(start_coords, 1) / atom_num

        do i = 2, conf_num
          curr_conf => coords(i, :, :)
          call restoreCoords(r(:, i), alpha(:, i), psi(:, i), conf_num,
     +      curr_conf)
          ! apply the rotation and the transformation
          curr_conf = matmul(curr_conf, rot_mat(i, :, :))
          curr_trans = sum(curr_conf, 1) / atom_num
          do j = 1, atom_num
            curr_conf(j, :) = curr_conf(j, :) - curr_trans + 
     +        first_trans
          end do
        end do

      return
      end subroutine trRestoreCoords


