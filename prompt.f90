! PROMPT Toolbox for MATLAB
!
! By Gaik Tamazian, 2016.
! gaik (dot) tamazian (at) gmail (dot) com

module prompt

  implicit none
  
  type :: TrModel
    real(kind=8), pointer :: r(:,:), alpha(:,:), psi(:,:)
    real(kind=8), pointer :: start_coords(:,:)
    real(kind=8), pointer :: atom_masses(:)
    real(kind=8), pointer :: rot_mat(:,:,:)
    integer(kind=4)       :: atom_num, conf_num
  end type TrModel

contains
  
  ! Procedure implementing the cross product of a pair of vectors
  pure function cross3d(a, b) result(c)
    real(kind=8), intent(in) :: a(3), b(3)
    real(kind=8) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end function cross3d

  pure function cross(a, b, n) result(c)
    integer(kind=4), intent(in)               :: n
    real(kind=8), dimension(n, 3), intent(in) :: a, b
    real(kind=8), dimension(n, 3) :: c

    integer(kind=4) :: i

    do i = 1, n
      c(i, 1) = a(i, 2) * b(i, 3) - a(i, 3) * b(i, 2)
      c(i, 2) = a(i, 3) * b(i, 1) - a(i, 1) * b(i, 3)
      c(i, 3) = a(i, 1) * b(i, 2) - a(i, 2) * b(i, 1)
    end do

  end function cross

  pure function dot(a, b, n) result(c)
    integer(kind=4), intent(in) :: n
    real(kind=8), dimension(n, 3), intent(in) :: a, b
    real(kind=8), dimension(n) :: c

    integer(kind=4) :: i

    do i = 1, n
      c(i) = dot_product(a(i, :), b(i, :))
    end do

  end function dot
 
  ! Procedure to restore Cartesian coordinates from internal ones
  ! for a single configuration
  pure function restoreCoords(r, alpha, psi, atom_num) result(coords)
    integer(kind=4), intent(in)  :: atom_num
    real(kind=8),    intent(in)  :: r(atom_num - 1)
    real(kind=8),    intent(in)  :: alpha(atom_num -2)
    real(kind=8),    intent(in)  :: psi(atom_num - 3)
    real(kind=8)                 :: coords(atom_num, 3)

    integer(kind=4)               :: i
    real(kind=8), dimension(3)    :: bc, n
    real(kind=8), dimension(3, 3) :: M

    coords(1, :) = [0.0d0, 0.0d0, 0.0d0]
    coords(2, :) = [r(1), 0.0d0, 0.0d0]
    coords(3, :) = coords(2, :) + [r(2) * cos(alpha(1)), r(2) * &
      sin(alpha(1)), 0.0d0]

    coords(4:, 1) = r(3:) * cos(alpha(2:))
    coords(4:, 2) = r(3:) * sin(alpha(2:)) * cos(psi)
    coords(4:, 3) = r(3:) * sin(alpha(2:)) * sin(psi)

    do i = 4, atom_num
      ! calculate the bc vector
      bc = coords(i-1, :) - coords(i-2, :)
      bc = bc / sqrt(sum(bc**2))
      ! calculate the n vector
      n = cross3d(coords(i-2, :) - coords(i-3, :), bc)
      n = n / sqrt(sum(n**2))
      ! calculate the M matrix
      M(:, 1) = bc
      M(:, 2) =  cross3d(n, bc)
      M(:, 3) = n
      ! get the point coordinates
      coords(i, :) = matmul(M, coords(i, :)) + coords(i-1, :)
    end do
  
  return
  end function restoreCoords

  ! Procedure to restore Cartesian coordinates for all
  ! configurations of a transformation.
  subroutine trRestoreCoords(m, coords, conf_coords)
    type(TrModel), intent(in)          :: m
    real(kind=8),  intent(out), target :: coords(m%atom_num, 3, &
      m%conf_num)
    real(kind=8),  intent(out)         :: conf_coords(m%atom_num, 3, &
      m%conf_num)

    integer(kind=4)            :: i, j
    real(kind=8), dimension(3) :: curr_trans, first_trans
    real(kind=8), pointer      :: curr_conf(:,:)

    coords(:, :, 1) = m%start_coords
    first_trans = sum(m%start_coords, 1) / m%atom_num

    do i = 2, m%conf_num
      curr_conf => coords(:, :, i)
      curr_conf =  restoreCoords(m%r(:, i), m%alpha(:, i), &
        m%psi(:, i), m%atom_num)
      conf_coords(:, :, i) = coords(:, :, i)
      ! apply the rotation and the transformation
      curr_conf = matmul(curr_conf, m%rot_mat(:, :, i))
      curr_trans = sum(curr_conf, 1) / m%atom_num
      do j = 1, m%atom_num
        curr_conf(j, :) = curr_conf(j, :) - curr_trans + &
          first_trans
      end do
    end do

  return
  end subroutine trRestoreCoords

  ! Procedure to calculate transtormation cost; note that it also
  ! returns Cartesian coordinates of configuration atoms restored from
  ! their internal coordinates
  subroutine trCost(m, p, cost_val, tr_coords, conf_coords)
    type(TrModel),   intent(in)  :: m
    integer(kind=4), intent(in)  :: p
    real(kind=8),    intent(out) :: cost_val
    real(kind=8),    intent(out) :: tr_coords(m%atom_num, 3, m%conf_num)
    real(kind=8),    intent(out) :: conf_coords(m%atom_num, 3, &
      m%conf_num)

    integer(kind=4)                     :: i
    real(kind=8), dimension(m%atom_num) :: temp_dist
    
    call trRestoreCoords(m, tr_coords, conf_coords)
      
    cost_val = 0
    do i = 1, m%conf_num - 1
      temp_dist = sqrt(sum((tr_coords(:, :, i + 1) - &
        tr_coords(:, :, i)) ** 2, 2))
      cost_val = cost_val + sum(m%atom_masses * (temp_dist ** p))
    end do

  return
  end subroutine trCost

  ! Procedure implementing the objective function
  subroutine objFunc(m, p_num, p_indices, t_num, t_indices, &
    angle_values, p, calc_grad, func_val, grad_vec)
    type(TrModel),   intent(in)                     :: m
    integer(kind=4), intent(in)                     :: p_num
    integer(kind=4), intent(in), dimension(p_num)   :: p_indices
    integer(kind=4), intent(in)                     :: t_num
    integer(kind=4), intent(in), dimension(t_num)   :: t_indices
    real(kind=8),    intent(in),  &
      dimension((p_num + t_num) * (m%conf_num - 2)) :: angle_values
    integer(kind=4), intent(in)                     :: p
    logical,         intent(in)                     :: calc_grad
    real(kind=8),    intent(out)                    :: func_val
    real(kind=8),    intent(out), &
      dimension((p_num + t_num) * (m%conf_num - 2)) :: grad_vec

    type(TrModel) :: temp_m
    integer       :: i
    real(kind=8), dimension(m%atom_num, 3, m%conf_num) :: coords
    real(kind=8), dimension(m%atom_num, 3, m%conf_num) :: conf_coords
    real(kind=8), dimension(p_num, m%conf_num - 2)     :: temp_p_angles
    real(kind=8), dimension(t_num, m%conf_num - 2)     :: temp_t_angles

    temp_m = m

    if (p_num > 0) then
      temp_p_angles = reshape( &
        angle_values(:p_num * (m%conf_num - 2)), &
        [p_num, m%conf_num - 2])
      do i = 1, p_num
        temp_m%alpha(p_indices(i), 2:(m%conf_num - 1)) = &
          temp_p_angles(i, :)
      end do
    end if

    if (t_num > 0) then
      temp_t_angles = reshape( &
        angle_values(p_num * (m%conf_num - 2) + 1:), &
        [t_num, m%conf_num - 2])
      do i = 1, t_num
        temp_m%psi(t_indices(i), 2:(m%conf_num - 1)) = &
          temp_t_angles(i, :)
      end do
    end if

    call trCost(temp_m, p, func_val, coords, conf_coords)

    if (calc_grad) then
      grad_vec = g(temp_m, p_num, p_indices, t_num, t_indices, coords, &
        conf_coords) 
    else
      grad_vec = 0
    end if
  
  return
  end subroutine objFunc

  ! The gradient of the transformation cost objective function
  pure function g(m, p_num, p_indices, t_num, t_indices, coords, &
    conf_coords)
    type(TrModel), intent(in)   :: m
    integer(kind=4), intent(in) :: p_num, t_num
    integer(kind=4), dimension(p_num), intent(in) :: p_indices
    integer(kind=4), dimension(t_num), intent(in) :: t_indices
    real(kind=8), dimension(m%atom_num, 3, m%conf_num), &
      intent(in) :: coords
    real(kind=8), dimension(m%atom_num, 3, m%conf_num), &
      intent(in) :: conf_coords
    real(kind=8), dimension((p_num + t_num) * (m%conf_num - 2)) :: g

    integer(kind=4) :: i, j, angle_num
    real(kind=8), dimension(m%atom_num, 3, m%conf_num) :: vS
    real(kind=8), dimension(m%atom_num - 1, 3, m%conf_num) :: vR
    real(kind=8), dimension(m%atom_num - 2, 3, m%conf_num) :: vN, vP
    real(kind=8), dimension(m%atom_num, 3, m%atom_num - 1) :: vQ
    real(kind=8), dimension(m%atom_num, 3) :: mean_vQ, v
    real(kind=8), dimension(p_num, m%conf_num - 2) :: grad_p
    real(kind=8), dimension(t_num, m%conf_num - 2) :: grad_t

    vS = s(coords, m%atom_num, m%conf_num)
    vR = r(conf_coords, m%atom_num, m%conf_num)
    vN = n(vR, m%atom_num, m%conf_num)
    vP = p(vR, m%atom_num, m%conf_num)

    do j = 2, m%conf_num - 1
      vQ = q(conf_coords(:, :, j), m%atom_num, 3)

      if (p_num > 0) then
        do i = 1, m%atom_num - 1
          vQ(i + 2, :, i) = conf_coords(i + 2, :, j) - &
            conf_coords(i + 1, :, j)
        end do
      end if

      do angle_num = 1, p_num
        i = p_indices(angle_num)
        mean_vQ = spread(sum(vQ(:, :, i) , 1) / m%atom_num, &
          1, m%atom_num)
        v = cross(spread(vP(i, :, j), 1, m%atom_num), vQ(:, :, i) - &
          mean_vQ, m%atom_num)
        grad_p(angle_num, j - 1) = 2 * sum(m%atom_masses * &
          dot(vS(:, :, j), matmul(v, m%rot_mat(:, :, j)), m%atom_num))
      end do

      if (p_num > 0) then
        do i = 1, m%atom_num - 1
          vQ(i + 2, :, i) = 0
        end do
      end if

      do angle_num = 1, t_num
        i = t_indices(angle_num)
        mean_vQ = spread(sum(vQ(:, :, i) , 1) / m%atom_num, &
          1, m%atom_num)
        v = cross(spread(vN(i, :, j), 1, m%atom_num), vQ(:, :, i) - &
          mean_vQ, m%atom_num)
        grad_t(angle_num, j - 1) = 2 * sum(m%atom_masses * &
          dot(vS(:, :, j), matmul(v, m%rot_mat(:, :, j)), m%atom_num))
      end do
    end do

    if (p_num > 0) then
      g(:p_num * (m%conf_num - 2)) = &
        reshape(grad_p, [p_num * (m%conf_num - 2)])
    end if
    g(p_num * (m%conf_num - 2) + 1:) = &
      reshape(grad_t, [t_num * (m%conf_num - 2)])

  end function g 

  pure function r(x, atom_num, conf_num)
    integer(kind=4), intent(in) :: atom_num, conf_num
    real(kind=8), dimension(atom_num, 3, conf_num), intent(in) :: x
    real(kind=8), dimension(atom_num - 1, 3, conf_num) :: r

    r = x(2:atom_num, :, :) - x(1:(atom_num - 1), :, :)

  end function r

  pure function n(r, atom_num, conf_num)
    integer(kind=4), intent(in) :: atom_num, conf_num
    real(kind=8), dimension(atom_num - 1, 3, conf_num), &
      intent(in) :: r
    real(kind=8), dimension(atom_num - 2, 3, conf_num) :: n

    real(kind=8), dimension(atom_num - 1, 3, conf_num) :: temp_n

    temp_n = r / spread(sqrt(sum(r**2, 2)), 2, 3)
    n = temp_n(2:, :, :) 

  end function n

  pure function p(r, atom_num, conf_num)
    integer(kind=4), intent(in) :: atom_num, conf_num
    real(kind=8), dimension(atom_num - 1, 3, conf_num), &
      intent(in) :: r
    real(kind=8), dimension(atom_num - 2, 3, conf_num) :: p

    integer(kind=4) :: i, j

    do i = 1, atom_num - 2
      do j = 1, conf_num
        p(i, :, j) = cross3d(r(i, :, j), r(i + 1, :, j)) 
        p(i, :, j) = p(i, :, j) / sqrt(sum(p(i, :, j) ** 2))
      end do
    end do

  end function p

  pure function s(x, atom_num, conf_num)
    integer(kind=4), intent(in) :: atom_num, conf_num
    real(kind=8), dimension(atom_num, 3, conf_num), intent(in) :: x
    real(kind=8), dimension(atom_num, 3, conf_num) :: s

    s(:, :, 2:(conf_num - 1)) = 2*x(:, :, 2:(conf_num - 1)) - &
      x(:, :, :(conf_num - 2)) - x(:, :, 3:)

  end function s

  pure function q(conf_x, atom_num, shift)
    integer(kind=4), intent(in) :: atom_num, shift
    real(kind=8), dimension(atom_num, 3), intent(in)   :: conf_x
    real(kind=8), dimension(atom_num, 3, atom_num - 1) :: q
    
    integer(kind=4) :: l, i
    
    q = 0
    do i = 1, atom_num - 1
      do l = i + shift,  atom_num 
        q(l, :, i) = conf_x(l, :) - conf_x(i + 1, :)
      end do
    end do

  end function q

end module prompt

