module prompt

  implicit none
  
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
    real(kind=8), intent(in)  :: r(:)
    real(kind=8), intent(in)  :: alpha(:)
    real(kind=8), intent(in)  :: psi(:)
    integer,      intent(in)  :: atom_num
    real(kind=8), intent(out) :: coords(:,:)

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
  subroutine trRestoreCoords(r, alpha, psi, atom_num, conf_num, &
    rot_mat, start_coords, coords)
    real(kind=8), intent(in),  allocatable :: r(:,:)
    real(kind=8), intent(in),  allocatable :: alpha(:,:)
    real(kind=8), intent(in),  allocatable :: psi(:,:)
    integer,      intent(in)               :: atom_num
    integer,      intent(in)               :: conf_num
    real(kind=8), intent(in),  allocatable :: rot_mat(:,:,:)
    real(kind=8), intent(in),  allocatable :: start_coords(:,:) 
    real(kind=8), intent(out), allocatable, target :: coords(:,:,:)

    integer                            :: i, j
    real(kind=8), dimension(3)         :: curr_trans, first_trans
    real(kind=8), pointer :: curr_conf(:,:)

    allocate(coords(conf_num, atom_num, 3))

    coords(1, :, :) = start_coords
    first_trans = sum(start_coords, 1) / atom_num

    do i = 2, conf_num
      curr_conf => coords(i, :, :)
      call restoreCoords(r(:, i), alpha(:, i), psi(:, i), conf_num, &
        curr_conf)
      ! apply the rotation and the transformation
      curr_conf = matmul(curr_conf, rot_mat(i, :, :))
      curr_trans = sum(curr_conf, 1) / atom_num
      do j = 1, atom_num
        curr_conf(j, :) = curr_conf(j, :) - curr_trans + &
          first_trans
      end do
    end do

  return
  end subroutine trRestoreCoords

  ! Procedure to calculate transtormation cost; note that it also
  ! returns Cartesian coordinates of configuration atoms restored from
  ! their internal coordinates
  subroutine trCost(r, alpha, psi, atom_masses, atom_num, conf_num, &
    rot_mat, start_coords, p, cost_val, tr_coords)
    real(kind=8), intent(in), allocatable  :: r(:,:)
    real(kind=8), intent(in), allocatable  :: alpha(:,:)
    real(kind=8), intent(in), allocatable  :: psi(:,:)
    real(kind=8), intent(in), allocatable  :: atom_masses(:)
    integer, intent(in)                    :: atom_num, conf_num, p
    real(kind=8), intent(in),  allocatable :: rot_mat(:,:,:)
    real(kind=8), intent(in),  allocatable :: start_coords(:,:) 
    real(kind=8), intent(out)              :: cost_val
    real(kind=8), intent(out), allocatable :: tr_coords(:,:,:)

    integer                                     :: i
    real(kind=8), dimension(atom_num, conf_num) :: distances
    
    call trRestoreCoords(r, alpha, psi, atom_num, conf_num, rot_mat, &
      start_coords, tr_coords)
      
    cost_val = 0
    do i = 2, conf_num
      distances(:, i) = sqrt(sum((tr_coords(i, :, :) - &
        tr_coords(i - 1, :, :)) ** 2, 2))
    end do

    cost_val = sum(atom_masses * sum(distances ** p, 2))
  
  return
  end subroutine trCost

end module prompt

