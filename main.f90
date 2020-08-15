program main
    use truss
    implicit none
    
    call read_input("example01.txt")
    !call print_data()
    call assemble_system_of_equation()
    call solve_system_of_equation()
    call compute_normal_forces()
    call print_results()

end program main
