! PROMPT Toolbox for MATLAB
!
! By Gaik Tamazian, 2016.
! gaik (dot) tamazian (at) gmail (dot) com

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
      mwPointer      :: model_ptr, p_index_ptr, t_index_ptr, angle_ptr
      mwPointer      :: masses_ptr, r_ptr, alpha_ptr, psi_ptr, 
     +  str_crd_ptr, end_crd_ptr, rot_mat_ptr
      
      ! Pointers to output data
      mwPointer      :: fval_ptr, grad_ptr

      ! Arguments for computational routines
      type(TrModel)   :: m
      mwPointer       :: atom_num, conf_num, p_num, t_num
      real(kind=8), allocatable, dimension(:) :: p_indices
      real(kind=8), allocatable, dimension(:) :: t_indices
      real(kind=8), allocatable, dimension(:) :: angle_values
      real(kind=8)                            :: func_val
      real(kind=8), allocatable, dimension(:) :: grad_vec
      
      ! Check for the proper number of function arguments
      if (nrhs .ne. 4) then
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

      if (mxIsNumeric(prhs(4)) .eq. 0) then
        call mexErrMsgIdAndTxt('MATLAB:trmobjfunc:NonNumeric',
     +                         'The fourth argument must be a numeric 
     +                          array.')
      endif

      ! Get pointers for input data
      model_ptr     = prhs(1)
      p_index_ptr   = mxGetPr(prhs(2))
      t_index_ptr   = mxGetPr(prhs(3))
      angle_ptr     = mxGetPr(prhs(4))

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
      p_num = mxGetM(prhs(2)) * mxGetN(prhs(2)) 
      t_num = mxGetM(prhs(3)) * mxGetN(prhs(3)) 

      ! Create result arrays
      plhs(1) = mxCreateDoubleMatrix(1, 1, 0.0)
      plhs(2) = mxCreateDoubleMatrix((p_num + t_num) * 
     +  (conf_num - 2), 1, 0.0)

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
      allocate(p_indices(p_num))
      allocate(angle_values((p_num + t_num) * (conf_num - 2)))
      allocate(t_indices(t_num))
      allocate(grad_vec((p_num + t_num) * (conf_num - 2)))

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
      call mxCopyPtrToReal8(p_index_ptr, p_indices, p_num)
      call mxCopyPtrToReal8(t_index_ptr, t_indices, t_num)
      call mxCopyPtrToReal8(angle_ptr, angle_values, 
     +  (p_num + t_num) * (conf_num - 2))

      ! Launch the computational routines
      call objFunc(m, int(p_num, 4), int(p_indices, 4),
     +  int(t_num, 4), int(t_indices, 4),
     +  angle_values, 2, nlhs > 1, func_val, grad_vec)

      ! Copy the output data back to MATLAB data structures
      call mxCopyReal8ToPtr(func_val, fval_ptr, 1)
      call mxCopyReal8ToPtr(grad_vec, grad_ptr, (p_num + t_num) *
     +  (conf_num - 2))

      ! Deallocate Fortran arrays
      deallocate(m%atom_masses)
      deallocate(m%r, m%alpha, m%psi)
      deallocate(m%start_coords, m%rot_mat)
      deallocate(angle_values, p_indices, t_indices)

      return
      end subroutine mexFunction


