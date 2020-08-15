module truss
    use cg
    implicit none

    integer(4)::nnodes, nelem
    ! Nodes' data
    real(8),dimension(:),allocatable::x,y,fx,fy
    integer(4),dimension(:),allocatable::rx,ry
    ! Elements' data
    real(8),dimension(:),allocatable::area,young
    integer(4),dimension(:,:),allocatable::connect
    ! Global system of equation
    real(8),dimension(:,:),allocatable::Kg
    real(8),dimension(:),allocatable::fg,disp
    ! Results
    real(8),dimension(:),allocatable::normal_force

contains

    subroutine read_input(filename)
        implicit none 
        character(*),intent(in)::filename
        character(20)::key
        integer(4)::i,index
        print "(2a)","Input file: ",filename
        open(unit=20,file=filename)
        read(20,*)key,nnodes
        read(20,*)key,nelem
        ! Allocating arrays
        allocate(x(nnodes),y(nnodes),rx(nnodes),ry(nnodes))
        allocate(fx(nnodes),fy(nnodes))
        allocate(connect(nelem,2),area(nelem),young(nelem))
        ! Reading nodes
        read(20,*)key
        if (key.eq."nodes") then
            write(*,"(a)",advance="no")"Reading nodes... "
            do i = 1,nnodes
                read(20,*)index,x(index),y(index),rx(index),ry(index), &
                fx(index),fy(index)
            enddo
            write(*,"(a)")"DONE!"
        else
            print *, "Error reading nodes data!"
        endif
        read(20,*)key
        if (key.eq."elements") then
            write(*,"(a)",advance="no")"Reading elements... "
            do i = 1,nelem
                read(20,*)index,connect(index,1),connect(index,2), &
                area(index),young(index)
            enddo
            write(*,"(a)")"DONE!"
        else
            print *, "Error reading elements data!"
        endif
    end subroutine read_input

    subroutine print_data()
        implicit none 
        integer(4)::i
        print *, "Truss data:"
        print *, "Nodes"
        do i=1,nnodes
            print 10, i, x(i), y(i), rx(i), ry(i), fx(i), fy(i)
        enddo
        print *, "Elements: "
        do i=1,nelem
            print 20, i, connect(i,1),connect(i,2),area(i),young(i)
        enddo

        10 format("Node # ",i2,", x=",f9.3,", y=",f9.3, &
        ", rx=",i2,", ry=",i2,", fx=",f9.3,", fy=",f9.3)
        20 format("Element # ",i2," from node # ",i2," to node # ",i2, &
        ", area = ",f9.3, ", Young modulus = ",es9.3)
    end subroutine print_data

    function length(element_number) result(l)
        implicit none
        integer(4),intent(in):: element_number 
        real(8)::l
        real(8)::xi,yi,xf,yf
        integer(4)::ni,nf
        ni = connect(element_number,1)
        nf = connect(element_number,2)
        xi = x(ni); yi = y(ni)
        xf = x(nf); yf = y(nf)        
        l = dsqrt((xf - xi) ** 2 + (yf - yi) ** 2)
    end function length

    subroutine sine_cosine(element_number, sine, cosine)
        implicit none
        integer(4), intent(in):: element_number
        real(8), intent(out):: sine, cosine
        real(8)::xi,yi,xf,yf,l
        integer(4)::ni,nf
        ni = connect(element_number,1)
        nf = connect(element_number,2)
        xi = x(ni); yi = y(ni)
        xf = x(nf); yf = y(nf)        
        l = length(element_number)
        sine = (yf - yi) / l
        cosine = (xf - xi) / l
    end subroutine sine_cosine

    function local_stiffness_matrix(element_number) result(k)
        implicit none
        integer(4), intent(in)::element_number
        real(8), dimension(4,4)::k
        real(8)::E,A,L,c,s,k_
        E = young(element_number)
        A = area(element_number)
        L = length(element_number)
        k_ = (E * A)/L
        call sine_cosine(element_number,s,c)
        k(1,1) = c ** 2
        k(1,2) = c * s; k(2,1) = k(1,2)
        k(1,3) = -(c ** 2); k(3,1) = k(1,3)
        k(1,4) = -c * s; k(4,1) = k(1,4)
        k(2,2) = s ** 2
        k(2,3) = -c * s; k(3,2) = k(2,3)
        k(2,4) = -(s ** 2); k(4,2) = k(2,4)
        k(3,3) = c ** 2
        k(3,4) = c * s; k(4,3) = k(3,4)
        k(4,4) = s ** 2
        k = k * k_
    end function local_stiffness_matrix

    subroutine assemble_system_of_equation()
        implicit none
        integer(4)::i,j,k
        integer(4),dimension(4)::indexes
        integer(4)::ni,nf,num
        real(8),dimension(4,4)::kl
        ! Allocating global stiffness matrix
        allocate(Kg(2*nnodes,2*nnodes),fg(2*nnodes),disp(2*nnodes))
        Kg = 0.0d0; fg = 0.0d0; disp = 0.0d0
        ! Assembling the global stiffness matrix
        do k = 1,nelem
            ni = connect(k, 1); nf = connect(k, 2)
            indexes(1) = 2 * ni - 1; indexes(2) = 2 * ni
            indexes(3) = 2 * nf - 1; indexes(4) = 2 * nf
            kl = local_stiffness_matrix(k)
            do i = 1,4
                do j = 1,4
                    Kg(indexes(i),indexes(j)) = Kg(indexes(i),indexes(j)) + kl(i,j) 
                enddo
            enddo
        enddo
        ! Assembling the global force vector
        do k = 1,nnodes
            fg(2*k-1) = fg(2*k-1) + fx(k)
            fg(2*k) = fg(2*k) + fy(k)
        enddo
        ! Applying Dirichlet Boundary Conditions
        do k = 1,nnodes
            if (rx(k).eq.1) then
                num = 2 * k - 1
                do i=1,2*nnodes
                    Kg(i,num) = 0.0d0
                    kg(num,i) = 0.0d0
                enddo
                Kg(num,num) = 1.0d0
                fg(num) = 0.0d0
            endif
            if (ry(k).eq.1) then
                num = 2 * k
                do i=1,2*nnodes
                    Kg(i,num) = 0.0d0
                    kg(num,i) = 0.0d0
                enddo
                Kg(num,num) = 1.0d0
                fg(num) = 0.0d0
            endif
         enddo
    end subroutine assemble_system_of_equation

    subroutine solve_system_of_equation()
        implicit none
        disp = solve_cg(Kg,fg)
    end subroutine solve_system_of_equation

    subroutine compute_normal_forces()
        implicit none
        integer(4)::k,ni,nf
        real(8)::E,A,L,c,s,ui,vi,uf,vf,di,df
        if (allocated(normal_force)) deallocate(normal_force)
        allocate(normal_force(nelem))
        do k = 1, nelem
            E = young(k)
            A = area(k)
            L = length(k)
            call sine_cosine(k,s,c)
            ni = connect(k,1)
            nf = connect(k,2)
            ui = disp(2*ni-1);vi = disp(2*ni)
            uf = disp(2*nf-1);vf = disp(2*nf)
            di = c * ui + s * vi
            df = c * uf + s * vf
            normal_force(k) =  E * A * (df - di) / L
        enddo

    end subroutine compute_normal_forces

    subroutine print_results()
        implicit none
        integer(4)::i
        print *," Results:"
        do i=1,nnodes
            print 10, i, disp(2*i-1), disp(2*i)
        enddo
        print *,""
        do i=1,nelem
            print 20, i, normal_force(i)
        enddo
        10 format(" Node #",i3,", u =",f10.6,", v = ",f10.6)
        20 format(" Element #",i3,", N = ",f9.3)

    end subroutine print_results

end module truss
