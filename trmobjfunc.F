#include "fintrf.h"

      ! Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      use prompt, only: objFunc, TrModel
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
      type(TrModel)   :: m
      mwPointer       :: atom_num, conf_num, angle_num
      real(kind=8), allocatable, dimension(:) :: angle_values
      real(kind=8), allocatable, dimension(:) :: angle_indices
      real(kind=8)                            :: func_val
      real(kind=8), allocatable, dimension(:) :: grad_vec
      
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
      index_ptr   = mxGetPr(prhs(2))
      angle_ptr   = mxGetPr(prhs(3))

      masses_ptr  = mxGetPr(mxGetField(model_ptr, 1, 'm'))
      r_ptr       = mxGetPr(mxGetField(model_ptr, 1, 'r'))
      alpha_ptr   = mxGetPr(mxGetField(model_ptr, 1, 'alpha'))
      psi_ptr     = mxGetPr(mxGetField(model_ptr, 1, 'psi'))
      str_crd_ptr = mxGetPr(mxGetField(model_ptr, 1, 'StartCoords'))
      rot_mat_ptr = mxGetPr(mxGetField(model_ptr, 1, 'U'))

      ! Determine the model parameters: the number of atoms and the
      ! number of configurations
      atom_num  = mxGetM(mxGetField(model_ptr, 1, 'r')) + 1
      conf_num  = mxGetN(mxGetField(model_ptr, 1, 'r'))
      angle_num = mxGetM(prhs(2)) * mxGetN(prhs(2)) 

      ! Create result arrays
      plhs(1) = mxCreateDoubleMatrix(1, 1, 0.0)
      plhs(2) = mxCreateDoubleMatrix(angle_num, 1, 0.0)

      ! Get pointers to output data
      fval_ptr = mxGetPr(plhs(1))
      grad_ptr = mxGetPr(plhs(2))

      ! Allocate Fortran arrays
      allocate(m%atom_masses(atom_num))
      allocate(m%r(atom_num - 1, conf_num))
      allocate(m%alpha(atom_num - 2, conf_num))
      allocate(m%psi(atom_num - 3, conf_num))
      allocate(m%start_coords(atom_num, 3))
      allocate(m%rot_mat(3, 3, conf_num))
      allocate(angle_values(angle_num * (conf_num - 2)))
      allocate(angle_indices(angle_num))
      allocate(grad_vec(angle_num))

      m%conf_num = int(conf_num, 4)
      m%atom_num = int(atom_num, 4)

      ! Copy the input data to Fortran arrays
      call mxCopyPtrToReal8(masses_ptr, m%atom_masses, atom_num)
      call mxCopyPtrToReal8(r_ptr, m%r, (atom_num - 1) * conf_num)
      call mxCopyPtrToReal8(alpha_ptr, m%alpha, (atom_num - 2) * 
     +  conf_num)
      call mxCopyPtrToReal8(psi_ptr, m%psi, (atom_num - 3) * 
     +  conf_num)
      call mxCopyPtrToReal8(str_crd_ptr, m%start_coords, 
     +  atom_num * 3)
      call mxCopyPtrToReal8(rot_mat_ptr, m%rot_mat, conf_num * 3 * 3)
      call mxCopyPtrToReal8(index_ptr, angle_indices, angle_num)
      call mxCopyPtrToReal8(angle_ptr, angle_values, 
     +  angle_num * (conf_num - 2))

      ! Launch the computational routines
      call objFunc(m, int(angle_num, 4), int(angle_indices, 4), 
     +  angle_values, 2, func_val, grad_vec)

      ! Copy the output data back to MATLAB data structures
      call mxCopyReal8ToPtr(func_val, fval_ptr, 1)
      call mxCopyReal8ToPtr(grad_vec, grad_ptr, angle_num)

      ! Deallocate Fortran arrays
      deallocate(m%atom_masses)
      deallocate(m%r, m%alpha, m%psi)
      deallocate(m%start_coords, m%rot_mat)
      deallocate(angle_values, angle_indices)

      return
      end subroutine mexFunction


