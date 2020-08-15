module cg

    real(8)::TOLERANCE=1.0d-6

contains

function solve_cg(A, b) result(x)
    implicit none
    ! Dummy Arguments
    real(8),dimension(:,:),intent(in)::A
    real(8),dimension(:),intent(in)::b
    real(8),dimension(size(b))::x
    ! Local variables
    integer(4)::k, maxit, t
    real(8),dimension(:),allocatable::r,z,p,Ap,d
    real(8)::alpha, beta, r_x_z, r_x_z_next
    real(8)::ti, tf
    t = size(b)
    maxit = int(1.1*t)
    allocate(r(t),d(t))
    ! First steps
    d = diag(A)
    x = b / d
    r = b - matmul(A,x)
    if (maxval(dabs(r)).lt.TOLERANCE) return
    allocate(z(t),p(t),Ap(t))
    k = 0
    z = r / d
    p = z
    call cpu_time(ti)
    r_x_z = dot_product(r,z)
    do while (k < maxit)
        Ap = matmul(A,p)
        alpha = r_x_z / dot_product(p,Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        if (maxval(dabs(r)).lt.TOLERANCE) exit
        z = r / d
        r_x_z_next = dot_product(r,z)
        beta = r_x_z_next / r_x_z
        p = z + beta * p
        k = k + 1
        r_x_z = r_x_z_next
    enddo
    call cpu_time(tf)
    ! print *, "residual vector -> ",r
    deallocate(r,p,z,d)
    if (k < maxit) then
        print 10, k,(tf-ti)
    else
        print 20, k
    endif
    10 format("Conjugate gradient solver converged after ",i5," iterations in ",f9.3," seconds!")
    20 format("Conjugate gradient solver did not converge in ",i5," iterations!")
end function solve_cg

function diag(A) result(d)
    implicit none
    real(8),dimension(:,:),intent(in)::A
    real(8),dimension(size(A,dim=1))::d
    integer(4)::i
    do i = 1,size(A,dim=1)
        d(i) = A(i,i)
    enddo
end function diag

end module cg
