! PROMPT Toolbox for MATLAB
!
! By Gaik Tamazian, 2016.
! mail (at) gtamazian (dot) com

#include "fintrf.h"

      ! Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      use prompt, only: gradDesc, TrModel
      implicit none

      mwPointer   :: plhs(*), prhs(*)
      integer     :: nlhs, nrhs

      ! Function declarations
      mwPointer   :: mxCreateDoubleMatrix, mxCreateDoubleScalar,
     +  mxCreateNumericArray
      integer*4   :: mxClassIDFromClassName
      mwPointer   :: mxGetPr, mxGetField
      mwPointer   :: mxGetM, mxGetN
      integer     :: mxIsStruct, mxIsNumeric

      ! Pointers to input data
      mwPointer   :: model_ptr, p_index_ptr, t_index_ptr, gradtol_ptr,
     +  maxsteps_ptr
      mwPointer   :: masses_ptr, r_ptr, alpha_ptr, psi_ptr, 
     +  str_crd_ptr, rot_mat_ptr
      
      ! Pointers to output data
      mwPointer   :: x_ptr, fval_ptr, gnorm_ptr, exitflag_ptr

      ! Arguments for compurational routines
      type(TrModel)   :: m
      mwPointer       :: atom_num, conf_num, p_num, t_num
      real(kind=8), allocatable, dimension(:) :: p_indices
      real(kind=8), allocatable, dimension(:) :: t_indices
      real(kind=8)                            :: gradtol
      real(kind=8)                            :: maxsteps
      real(kind=8), allocatable, dimension(:) :: x
      real(kind=8)                            :: fval, gnorm
      integer                                 :: exitflag
      real(kind=8), allocatable, dimension(:,:) :: temp_coords

      ! Check for the proper number of function arguments
      if (nrhs .ne. 5) then
        call mexErrMsgIdAndTxt('MATLAB:trmoptim:nInput',
     +                         'Three input arguments required.')
      endif

      ! Check input types
      if (mxIsStruct(prhs(1)) .eq. 0) then
        call mexErrMsgIdAndTxt('MATLAB:trmoptim:NonStructure',
     +                         'The first argument must be a structure
     +                          representing a transformation model.')
      endif

      if (mxIsNumeric(prhs(2)) .eq. 0) then
        call mexErrMsgIdAndTxt('MATLAB:trmoptim:NonNumeric',
     +                         'The second argument must be a numeric 
     +                          array.')
      endif

      if (mxIsNumeric(prhs(3)) .eq. 0) then
        call mexErrMsgIdAndTxt('MATLAB:trmoptim:NonNumeric',
     +                         'The third argument must be a numeric 
     +                          array.')
      endif

      if (mxIsNumeric(prhs(4)) .eq. 0) then
        call mexErrMsgIdAndTxt('MATLAB:trmoptim:NonNumeric',
     +                         'The fourth argument must be a real 
     +                          value.')
      endif

      if (mxIsNumeric(prhs(5)) .eq. 0) then
        call mexErrMsgIdAndTxt('MATLAB:trmoptim:NonNumeric',
     +                         'The fifth argument must be an integer 
     +                          value.')
      endif

      ! Get pointers for input data
      model_ptr     = prhs(1)
      p_index_ptr   = mxGetPr(prhs(2))
      t_index_ptr   = mxGetPr(prhs(3))
      gradtol_ptr   = mxGetPr(prhs(4))
      maxsteps_ptr  = mxGetPr(prhs(5))

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
      plhs(1) = mxCreateDoubleMatrix((p_num + t_num) *
     +  (conf_num - 2), 1, 0.0)
      plhs(2) = mxCreateDoubleMatrix(1, 1, 0.0)
      plhs(3) = mxCreateDoubleMatrix(1, 1, 0.0)
      plhs(4) = mxCreateNumericArray(1, 1, 
     +  mxClassIDFromClassName('int32'), 0);

      ! Get pointers to output data
      x_ptr        = mxGetPr(plhs(1))
      fval_ptr     = mxGetPr(plhs(2))
      gnorm_ptr    = mxGetPr(plhs(3))
      exitflag_ptr = mxGetPr(plhs(4))

      ! Allocate Fortran arrays
      allocate(m%atom_masses(atom_num))
      allocate(m%r(atom_num - 1, conf_num))
      allocate(m%alpha(atom_num - 2, conf_num))
      allocate(m%psi(atom_num - 3, conf_num))
      allocate(m%start_coords(3, atom_num))
      allocate(temp_coords(atom_num, 3))
      allocate(m%rot_mat(3, 3, conf_num))
      allocate(p_indices(p_num))
      allocate(t_indices(t_num))
      allocate(x((p_num + t_num) * (conf_num - 2)))

      m%conf_num = int(conf_num, 4)
      m%atom_num = int(atom_num, 4)

      ! Copy the input data to Fortran arrays
      call mxCopyPtrToReal8(masses_ptr, m%atom_masses, atom_num)
      call mxCopyPtrToReal8(r_ptr, m%r, (atom_num - 1) * conf_num)
      call mxCopyPtrToReal8(alpha_ptr, m%alpha, (atom_num - 2) * 
     +  conf_num)
      call mxCopyPtrToReal8(psi_ptr, m%psi, (atom_num - 3) * 
     +  conf_num)
      call mxCopyPtrToReal8(str_crd_ptr, temp_coords,
     +  atom_num * 3)
      m%start_coords = transpose(temp_coords)
      call mxCopyPtrToReal8(rot_mat_ptr, m%rot_mat, conf_num * 3 * 3)
      call mxCopyPtrToReal8(p_index_ptr, p_indices, p_num)
      call mxCopyPtrToReal8(t_index_ptr, t_indices, t_num)
      call mxCopyPtrToReal8(gradtol_ptr, gradtol, 1)
      call mxCopyPtrToReal8(maxsteps_ptr, maxsteps, 1)

      ! Launch the cimputational routine
      call gradDesc(m, int(p_num, 4), int(p_indices, 4),
     +  int(t_num, 4), int(t_indices, 4),
     +  gradtol, int(maxsteps, 4), x, fval, gnorm, exitflag)

      ! Copy the output data back to MATLAB data structures
      call mxCopyReal8ToPtr(x, x_ptr, (p_num + t_num) *
     +  (conf_num - 2))
      call mxCopyReal8ToPtr(fval, fval_ptr, 1)
      call mxCopyReal8ToPtr(gnorm, gnorm_ptr, 1)
      call mxCopyInteger4ToPtr(exitflag, exitflag_ptr, 1);

      ! Deallocate Fortran arrays
      deallocate(m%atom_masses)
      deallocate(m%r, m%alpha, m%psi)
      deallocate(m%start_coords, m%rot_mat)
      deallocate(temp_coords)
      deallocate(x, p_indices, t_indices)

      return
      end subroutine mexFunction


