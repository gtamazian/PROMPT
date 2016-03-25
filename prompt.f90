module prompt

  implicit none
  
  type :: TrModel
    real(kind=8), allocatable :: r(:,:), alpha(:,:), psi(:,:)
    real(kind=8), allocatable :: start_coords(:,:)
    real(kind=8), allocatable :: atom_masses(:)
    real(kind=8), allocatable :: rot_mat(:,:,:)
    integer                   :: atom_num, conf_num
  end type TrModel

contains
  
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
    real(kind=8), intent(in)  :: r(atom_num - 1)
    real(kind=8), intent(in)  :: alpha(atom_num -2)
    real(kind=8), intent(in)  :: psi(atom_num - 3)
    integer,      intent(in)  :: atom_num
    real(kind=8), intent(out) :: coords(atom_num, 3)

    integer                       :: i
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
      bc = bc / norm2(bc)
      ! calculate the n vector
      call cross3d(coords(i-2, :) - coords(i-3, :), bc, n)
      n = n / norm2(n)
      ! calculate the M matrix
      M(:, 1) = bc
      call cross3d(n, bc, M(:, 2))
      M(:, 3) = n
      ! get the point coordinates
      coords(i, :) = matmul(M, coords(i, :)) + coords(i-1, :)
    end do
  
  return
  end subroutine restoreCoords

  ! Procedure to restore Cartesian coordinates for all
  ! configurations of a transformation.
  subroutine trRestoreCoords(m, coords)
    type(TrModel), intent(in)          :: m
    real(kind=8),  intent(out), target :: &
        coords(m%conf_num, m%atom_num, 3)

    integer                            :: i, j
    real(kind=8), dimension(3)         :: curr_trans, first_trans
    real(kind=8), pointer :: curr_conf(:,:)

    coords(1, :, :) = m%start_coords
    first_trans = sum(m%start_coords, 1) / m%atom_num

    do i = 2, m%conf_num
      curr_conf => coords(i, :, :)
      call restoreCoords(m%r(:, i), m%alpha(:, i), &
        m%psi(:, i), m%conf_num, curr_conf)
      ! apply the rotation and the transformation
      curr_conf = matmul(curr_conf, m%rot_mat(i, :, :))
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
  subroutine trCost(m, p, cost_val, tr_coords)
    type(TrModel), intent(in) :: m
    integer, intent(in)       :: p
    real(kind=8), intent(out) :: cost_val
    real(kind=8), intent(out) :: tr_coords(m%conf_num, m%atom_num, 3)

    integer                                         :: i
    real(kind=8), dimension(m%atom_num, m%conf_num) :: distances
    
    call trRestoreCoords(m, tr_coords)
      
    cost_val = 0
    do i = 2, m%conf_num
      distances(:, i) = sqrt(sum((tr_coords(i, :, :) - &
        tr_coords(i - 1, :, :)) ** 2, 2))
    end do

    cost_val = sum(m%atom_masses * sum(distances ** p, 2))
  
  return
  end subroutine trCost

end module prompt

