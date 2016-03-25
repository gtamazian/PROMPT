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


