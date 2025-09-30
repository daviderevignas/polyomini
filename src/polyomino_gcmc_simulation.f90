subroutine polyomino_gcmc_simulation(n_mc_steps, n_savings, &
    & n_omino_types, ominos_sizes, max_omino_size, omino_shapes, &
    & initial_n_ominos_per_shape, initial_n_ominos, initial_positions,initial_orientations,&
    & max_n_ominos, chemical_potentials, &
    & move_probabilities, box_size, prng_seed, &
    & system_snapshots, history_positions, history_orientations, history_n_ominos,&
    & history_vec_omino_types)
    implicit none

    ! Inputs
    integer, intent(in) :: n_mc_steps, n_savings
    integer, intent(in) :: prng_seed
    integer, intent(in) :: n_omino_types
    integer, intent(in) :: ominos_sizes(n_omino_types)
    integer, intent(in) :: max_omino_size
    integer, intent(in) :: omino_shapes(max_omino_size,2,n_omino_types)
    integer, intent(in) :: initial_n_ominos_per_shape(n_omino_types)
    integer, intent(in) :: initial_n_ominos
    integer, intent(in) :: max_n_ominos
    double precision, intent(in) :: chemical_potentials(n_omino_types)
    double precision, intent(in) :: move_probabilities(4)
    integer, intent(in) :: box_size
    integer, intent(in) :: initial_positions(initial_n_ominos,2)
    integer, intent(in) :: initial_orientations(initial_n_ominos)

    ! Outputs
    integer, intent(out) :: system_snapshots(box_size,box_size,n_savings)
    integer, intent(out) :: history_positions(max_n_ominos,2,n_savings)
    integer, intent(out) :: history_orientations(max_n_ominos,n_savings)
    integer, intent(out) :: history_n_ominos(n_savings)
    integer, intent(out) :: history_vec_omino_types(max_n_ominos,n_savings)
    
    ! Local variables
    integer :: i_step, i_omino_type, i_omino, i_omino_of_type, i_cell
    integer :: ix_cell, iy_cell
    integer :: rand_omino, rand_x, rand_y, rand_rotation, rand_type
    integer :: n_seed
    integer :: n_lattice_sites 
    integer :: current_system_state(box_size,box_size)

    integer, allocatable :: seed(:)

    integer, allocatable :: current_positions(:,:)
    integer, allocatable :: current_orientations(:)
    integer, allocatable :: vec_omino_types(:)

    integer :: current_n_ominos
    integer :: current_n_ominos_per_type(n_omino_types)

    integer :: save_every
    logical :: has_overlap

    integer :: rot_omino_shapes(max_omino_size,2,n_omino_types,4)
    integer :: orig_x, orig_y

    double precision :: rand_real, cum_prob
    double precision :: acc_prob, rand_accept

    ! Initialize variables
    save_every = n_mc_steps/n_savings

    ! PRNG initialization
    call random_seed(size=n_seed)
    allocate(seed(n_seed))
    seed = prng_seed                  
    call random_seed(put=seed) 
    deallocate(seed)

    n_lattice_sites = box_size**2

    ! Allocate dynamic arrays
    allocate(current_positions(max_n_ominos,2))
    allocate(current_orientations(max_n_ominos))
    allocate(vec_omino_types(max_n_ominos))

    ! Initialize based on input
    current_n_ominos = initial_n_ominos
    current_n_ominos_per_type = initial_n_ominos_per_shape
    history_vec_omino_types = -1

    ! Copy initial positions and orientations, set up vec_omino_types
    i_omino = 0
    do i_omino_type = 1, n_omino_types
        do i_omino_of_type = 1, initial_n_ominos_per_shape(i_omino_type)
            i_omino = i_omino + 1
            if (i_omino > initial_n_ominos) stop 'Error: initial_n_ominos mismatch'
            vec_omino_types(i_omino) = i_omino_type
            current_positions(i_omino,:) = initial_positions(i_omino,:)
            current_orientations(i_omino) = initial_orientations(i_omino)
        end do
    end do
    if (i_omino /= initial_n_ominos) stop 'Error: initial_n_ominos does not match sum of initial_n_ominos_per_shape'

    ! Precompute rotated shapes
    rot_omino_shapes = 0
    do i_omino_type = 1, n_omino_types
        do i_cell = 1, ominos_sizes(i_omino_type)
            orig_x = omino_shapes(i_cell, 1, i_omino_type)
            orig_y = omino_shapes(i_cell, 2, i_omino_type)
            
            ! Rotation 1: 0 degrees clockwise (x,y) -> (x,y)
            rot_omino_shapes(i_cell, 1, i_omino_type, 1) = orig_x
            rot_omino_shapes(i_cell, 2, i_omino_type, 1) = orig_y

            ! Rotation 2: 90 degrees clockwise (x,y) -> (-y,x)
            rot_omino_shapes(i_cell, 1, i_omino_type, 2) = -orig_y
            rot_omino_shapes(i_cell, 2, i_omino_type, 2) = orig_x
            
            ! Rotation 3: 180 degrees (x,y) -> (-x,-y)
            rot_omino_shapes(i_cell, 1, i_omino_type, 3) = -orig_x
            rot_omino_shapes(i_cell, 2, i_omino_type, 3) = -orig_y
            
            ! Rotation 4: 270 degrees clockwise (x,y) -> (y,-x)
            rot_omino_shapes(i_cell, 1, i_omino_type, 4) = orig_y
            rot_omino_shapes(i_cell, 2, i_omino_type, 4) = -orig_x
        end do
    end do

    ! Initialize system state
    current_system_state = 0
    do i_step = 0, n_mc_steps-1

        ! Save snapshot
        if (mod(i_step, save_every) == 0) then
            current_system_state = 0
            
            do i_omino = 1, current_n_ominos
                i_omino_type = vec_omino_types(i_omino)
                do i_cell = 1, ominos_sizes(i_omino_type)
                    ix_cell = modulo(current_positions(i_omino,1) + & 
                        & rot_omino_shapes(i_cell,1,i_omino_type,current_orientations(i_omino)) - 1, box_size) + 1
                    iy_cell = modulo(current_positions(i_omino,2) + &
                        & rot_omino_shapes(i_cell,2,i_omino_type,current_orientations(i_omino)) - 1, box_size) + 1
                    current_system_state(ix_cell, iy_cell) = 1
                end do
            end do

            system_snapshots(:,:,i_step/save_every+1) = current_system_state
            history_positions(:,:,i_step/save_every+1) = -1
            history_positions(1:current_n_ominos,:,i_step/save_every+1) = current_positions(1:current_n_ominos,:)
            history_orientations(:,i_step/save_every+1) = -1
            history_orientations(1:current_n_ominos,i_step/save_every+1) = current_orientations(1:current_n_ominos)
            history_n_ominos(i_step/save_every+1) = current_n_ominos
            history_vec_omino_types(1:current_n_ominos,i_step/save_every+1)=vec_omino_types(1:current_n_ominos)
        end if

        ! Select move type
        call random_number(rand_real)
        cum_prob = 0.0d0
        if (rand_real < move_probabilities(1)) then
            ! Translation
            if (current_n_ominos == 0) cycle
            call generate_randi(1, current_n_ominos, rand_omino)
            call generate_randi(1, box_size, rand_x)
            call generate_randi(1, box_size, rand_y)

            call check_overlap_after_translation(rand_omino, rand_x, rand_y, has_overlap)

            if (.not. has_overlap) then
                current_positions(rand_omino,1) = rand_x
                current_positions(rand_omino,2) = rand_y
            end if

        else
            cum_prob = move_probabilities(1)
            if (rand_real < cum_prob + move_probabilities(2)) then
                ! Rotation
                if (current_n_ominos == 0) cycle
                call generate_randi(1, current_n_ominos, rand_omino)
                call generate_randi(1, 4, rand_rotation)

                call check_overlap_after_rotation(rand_omino, rand_rotation, has_overlap)

                if (.not. has_overlap) then
                    current_orientations(rand_omino) = rand_rotation
                end if

            else
                cum_prob = cum_prob + move_probabilities(2)
                if (rand_real < cum_prob + move_probabilities(3)) then
                    ! Insertion
                    call generate_randi(1, n_omino_types, rand_type)
                    if (current_n_ominos >= max_n_ominos) cycle
                    call generate_randi(1, box_size, rand_x)
                    call generate_randi(1, box_size, rand_y)
                    call generate_randi(1, 4, rand_rotation)

                    call check_overlap_for_insertion(rand_type, rand_x, rand_y, rand_rotation, has_overlap)

                    if (.not. has_overlap) then
                        
                        acc_prob = exp(chemical_potentials(rand_type)) / dble(current_n_ominos_per_type(rand_type) + 1)
                        call random_number(rand_accept)
                        if (rand_accept < min(1.0d0, acc_prob)) then
                            current_n_ominos = current_n_ominos + 1
                            current_positions(current_n_ominos,1) = rand_x
                            current_positions(current_n_ominos,2) = rand_y
                            current_orientations(current_n_ominos) = rand_rotation
                            vec_omino_types(current_n_ominos) = rand_type
                            current_n_ominos_per_type(rand_type) = current_n_ominos_per_type(rand_type) + 1
                        end if
                    end if

                else
                    ! Deletion
                    if (current_n_ominos == 0) cycle
                    call generate_randi(1, current_n_ominos, rand_omino)
                    rand_type = vec_omino_types(rand_omino)
                    acc_prob = dble(current_n_ominos_per_type(rand_type)) / exp(chemical_potentials(rand_type))
                    call random_number(rand_accept)
                    if (rand_accept < min(1.0d0, acc_prob)) then
                        ! Remove by swapping with last and decrementing
                        current_positions(rand_omino,:) = current_positions(current_n_ominos,:)
                        current_orientations(rand_omino) = current_orientations(current_n_ominos)
                        vec_omino_types(rand_omino) = vec_omino_types(current_n_ominos)
                        current_n_ominos = current_n_ominos - 1
                        current_n_ominos_per_type(rand_type) = current_n_ominos_per_type(rand_type) - 1
                    end if
                end if
            end if
        end if
    end do

    ! Deallocate dynamic arrays
    deallocate(current_positions)
    deallocate(current_orientations)
    deallocate(vec_omino_types)

    contains

    subroutine generate_randi(min_val, max_val, randi)
        implicit none
        integer, intent(in) :: min_val, max_val
        integer, intent(out) :: randi
        double precision :: rand_real

        call random_number(rand_real)
        randi = min_val + int(rand_real * (max_val - min_val + 1)) ! [[min, max]] inclusive
    end subroutine generate_randi

    subroutine check_overlap_after_translation(omino_id, new_x, new_y, has_overlap)
        implicit none
        integer, intent(in) :: omino_id, new_x, new_y
        logical, intent(out) :: has_overlap
        
        integer :: i_omino, i_omino_type, i_cell, j_cell
        integer :: ix_cell, iy_cell, jx_cell, jy_cell
        integer :: other_omino_type
        integer :: delta_x, delta_y

        has_overlap = .false.
        
        i_omino_type = vec_omino_types(omino_id)
        
        do i_omino = 1, current_n_ominos
            if (i_omino == omino_id) cycle

            delta_x = current_positions(i_omino,1) - new_x
            if (delta_x > box_size/2) then
                delta_x = delta_x - box_size
            else if (delta_x < -box_size/2) then
                delta_x = delta_x + box_size
            end if

            if (abs(delta_x) < 2*max_omino_size) then 
                delta_y = current_positions(i_omino,2) - new_y
                if (delta_y > box_size/2) then
                    delta_y = delta_y - box_size
                else if (delta_y < -box_size/2) then
                    delta_y = delta_y + box_size
                end if

                if (abs(delta_y) < 2*max_omino_size) then 
                    other_omino_type = vec_omino_types(i_omino)
                    
                    do i_cell = 1, ominos_sizes(i_omino_type)
                        ix_cell = modulo(new_x + & 
                            & rot_omino_shapes(i_cell,1,i_omino_type,current_orientations(omino_id)) - 1, box_size) + 1
                        iy_cell = modulo(new_y + &
                            & rot_omino_shapes(i_cell,2,i_omino_type,current_orientations(omino_id)) - 1, box_size) + 1
                        
                        do j_cell = 1, ominos_sizes(other_omino_type)
                            jx_cell = modulo(current_positions(i_omino,1) + &
                                & rot_omino_shapes(j_cell,1,other_omino_type,current_orientations(i_omino)) - 1, box_size) + 1
                            jy_cell = modulo(current_positions(i_omino,2) + & 
                                & rot_omino_shapes(j_cell,2,other_omino_type,current_orientations(i_omino)) - 1, box_size) + 1
                            
                            if (ix_cell == jx_cell .and. iy_cell == jy_cell) then
                                has_overlap = .true.
                                return
                            end if
                        end do
                    end do
                end if
            end if
        end do
    end subroutine check_overlap_after_translation

    subroutine check_overlap_after_rotation(omino_id, new_rotation, has_overlap)
        implicit none
        integer, intent(in) :: omino_id, new_rotation
        logical, intent(out) :: has_overlap
        
        integer :: i_omino, i_omino_type, i_cell, j_cell
        integer :: ix_cell, iy_cell, jx_cell, jy_cell
        integer :: other_omino_type
        integer :: delta_x, delta_y
        
        has_overlap = .false.
        
        i_omino_type = vec_omino_types(omino_id)
        
        do i_omino = 1, current_n_ominos
            if (i_omino == omino_id) cycle
            
            other_omino_type = vec_omino_types(i_omino)

            delta_x = current_positions(i_omino,1) - current_positions(omino_id,1)
            if (delta_x > box_size/2) then
                delta_x = delta_x - box_size
            else if (delta_x < -box_size/2) then
                delta_x = delta_x + box_size
            end if

            if (abs(delta_x) < 2*max_omino_size) then 
                delta_y = current_positions(i_omino,2) - current_positions(omino_id,2)
                if (delta_y > box_size/2) then
                    delta_y = delta_y - box_size
                else if (delta_y < -box_size/2) then
                    delta_y = delta_y + box_size
                end if

                if (abs(delta_y) < 2*max_omino_size) then 
                    do i_cell = 1, ominos_sizes(i_omino_type)
                        ix_cell = modulo(current_positions(omino_id,1) + rot_omino_shapes(i_cell,1,i_omino_type,new_rotation) - 1, box_size) + 1
                        iy_cell = modulo(current_positions(omino_id,2) + rot_omino_shapes(i_cell,2,i_omino_type,new_rotation) - 1, box_size) + 1
                        
                        do j_cell = 1, ominos_sizes(other_omino_type)
                            jx_cell = modulo(current_positions(i_omino,1) + rot_omino_shapes(j_cell,1,other_omino_type,current_orientations(i_omino)) - 1, box_size) + 1
                            jy_cell = modulo(current_positions(i_omino,2) + rot_omino_shapes(j_cell,2,other_omino_type,current_orientations(i_omino)) - 1, box_size) + 1
                            
                            if (ix_cell == jx_cell .and. iy_cell == jy_cell) then
                                has_overlap = .true.
                                return
                            end if
                        end do
                    end do
                end if
            end if
        end do
    end subroutine check_overlap_after_rotation

    subroutine check_overlap_for_insertion(omino_type, new_x, new_y, new_rotation, has_overlap)
        implicit none
        integer, intent(in) :: omino_type, new_x, new_y, new_rotation
        logical, intent(out) :: has_overlap
        
        integer :: i_omino, i_cell, j_cell
        integer :: ix_cell, iy_cell, jx_cell, jy_cell
        integer :: other_omino_type
        integer :: delta_x, delta_y
        
        has_overlap = .false.
        
        do i_omino = 1, current_n_ominos
            other_omino_type = vec_omino_types(i_omino)

            delta_x = current_positions(i_omino,1) - new_x
            if (delta_x > box_size/2) then
                delta_x = delta_x - box_size
            else if (delta_x < -box_size/2) then
                delta_x = delta_x + box_size
            end if

            if (abs(delta_x) < 2*max_omino_size) then 
                delta_y = current_positions(i_omino,2) - new_y
                if (delta_y > box_size/2) then
                    delta_y = delta_y - box_size
                else if (delta_y < -box_size/2) then
                    delta_y = delta_y + box_size
                end if

                if (abs(delta_y) < 2*max_omino_size) then 
                    do i_cell = 1, ominos_sizes(omino_type)
                        ix_cell = modulo(new_x + rot_omino_shapes(i_cell,1,omino_type,new_rotation) - 1, box_size) + 1
                        iy_cell = modulo(new_y + rot_omino_shapes(i_cell,2,omino_type,new_rotation) - 1, box_size) + 1
                        
                        do j_cell = 1, ominos_sizes(other_omino_type)
                            jx_cell = modulo(current_positions(i_omino,1) + rot_omino_shapes(j_cell,1,other_omino_type,current_orientations(i_omino)) - 1, box_size) + 1
                            jy_cell = modulo(current_positions(i_omino,2) + rot_omino_shapes(j_cell,2,other_omino_type,current_orientations(i_omino)) - 1, box_size) + 1
                            
                            if (ix_cell == jx_cell .and. iy_cell == jy_cell) then
                                has_overlap = .true.
                                return
                            end if
                        end do
                    end do
                end if
            end if
        end do
    end subroutine check_overlap_for_insertion

end subroutine polyomino_gcmc_simulation