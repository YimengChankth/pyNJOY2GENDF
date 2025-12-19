module csv_nd_with_index
    use, intrinsic :: iso_fortran_env
    implicit none
    contains

    subroutine write_csv_nd_with_index(unit, A, shp)
        integer, intent(in) :: unit
        real(real64), intent(in) :: A(*)       ! 1D linearized array
        integer, intent(in) :: shp(:)  ! shape of original N-D array
        integer :: n, r, i, k
        integer, allocatable :: idx(:)

        r = size(shp)
        n = product(shp)  ! total number of elements

        allocate(idx(r))

        do i = 1, n
            call ind2sub(shp, i, idx)
            ! write(*,*)idx
            do k = 1, size(idx)
                write(unit, '(i0, ",")', advance="no") idx(k)
            end do
            write(unit, '(1pe14.7)') A(i)
        end do

        deallocate(idx)
    end subroutine write_csv_nd_with_index

    !------------------------------------------------------------
    subroutine ind2sub(shp, lin, idx)
        integer, intent(in) :: shp(:)
        integer, intent(in) :: lin
        integer, intent(out) :: idx(:)
        integer :: k, r, m, t

        r = size(shp)
        t = lin - 1
        do k = 1, r
            m = 1
            if (k < r) m = product(shp(k+1:r))
            idx(k) = t / m + 1
            t = mod(t, m)
        end do
    end subroutine ind2sub


    ! ym note: seems rather silly to do  this
    subroutine write_csv_3D(unit, A)
        integer, intent(in) :: unit
        real(real64), intent(in) :: A(:,:,:)       
        call write_csv_nd_with_index(unit, reshape(A, [size(A)]), shape(A))
    end subroutine write_csv_3D


    subroutine write_csv_2D(unit, A)
        integer, intent(in) :: unit
        real(real64), intent(in) :: A(:,:)       
        call write_csv_nd_with_index(unit, reshape(A, [size(A)]), shape(A))
    end subroutine write_csv_2D

    subroutine write_csv_1D(unit, A)
        integer, intent(in) :: unit
        real(real64), intent(in) :: A(:)       
        call write_csv_nd_with_index(unit, reshape(A, [size(A)]), shape(A))
    end subroutine write_csv_1D

end module csv_nd_with_index
