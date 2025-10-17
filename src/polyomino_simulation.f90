subroutine polyomino_simulation(n_mc_steps, n_savings, &
    & n_omino_types, ominos_sizes, max_omino_size, omino_shapes, &
    & n_ominos_per_shape, n_ominos, initial_positions, initial_orientations, &
    & box_size, &
    & prng_seed, &
    & system_snapshots, history_positions, history_orientations)
    implicit none
    
    ! Inputs
    integer, intent(in) :: n_mc_steps, n_savings
    integer, intent(in) :: prng_seed
    integer, intent(in) :: n_omino_types
    integer, intent(in) :: ominos_sizes(n_omino_types)
    integer, intent(in) :: max_omino_size
    integer, intent(in) :: omino_shapes(max_omino_size,2,n_omino_types)
    integer, intent(in) :: n_ominos_per_shape(n_omino_types)

    integer, intent(in) :: n_ominos
    integer, intent(in) :: initial_positions(n_ominos,2)
    integer, intent(in) :: initial_orientations(n_ominos)

    integer, intent(in) :: box_size


    integer, intent(out) :: system_snapshots(box_size,box_size,n_savings)
    integer, intent(out) :: history_positions(n_ominos,2,n_savings)
    integer, intent(out) :: history_orientations(n_ominos,n_savings)
    
    integer :: i_step, i_omino_type, i_omino, i_omino_of_type, i_cell
    integer :: ix_cell, iy_cell
    integer :: rand_omino, rand_x, rand_y, rand_rotation
    integer :: n_seed
    integer :: n_lattice_sites 
    integer :: current_system_state(box_size,box_size)

    integer, allocatable :: seed(:)

    integer :: current_positions(n_ominos,2)
    integer :: current_orientations(n_ominos)

    integer :: vec_omino_types(n_ominos)

    integer :: save_every
    logical :: has_overlap

    integer :: rot_omino_shapes(max_omino_size,2,n_omino_types,4)
    integer :: orig_x, orig_y




    save_every = n_mc_steps/n_savings

    !PRNG initialization
    call random_seed(size=n_seed)
    allocate(seed(n_seed))
    seed = prng_seed                  
    call random_seed(put=seed) 
    deallocate(seed)

    n_lattice_sites = box_size**2
    current_positions=initial_positions
    current_orientations=initial_orientations

    i_omino = 0
    do i_omino_type = 1, n_omino_types
        do i_omino_of_type = 1,  n_ominos_per_shape(i_omino_type)
            i_omino = i_omino + 1
            vec_omino_types(i_omino)=i_omino_type
        end do
    end do

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

    current_system_state = 0
    do i_step = 0, n_mc_steps-1

        !---------------- SAVE SNAPSHOT ----------------
        if (mod(i_step, save_every) == 0) then
            current_system_state = 0
            
            do i_omino = 1, n_ominos
                i_omino_type=vec_omino_types(i_omino)
                do i_cell=1,ominos_sizes(i_omino_type)

                    ! ix_cell = current_positions(i_omino,1)+omino_shapes(i_cell,1,i_omino_type)
                    ! iy_cell = current_positions(i_omino,2)+omino_shapes(i_cell,2,i_omino_type)

                    ix_cell = modulo(current_positions(i_omino,1) + & 
                        & rot_omino_shapes(i_cell,1,i_omino_type,current_orientations(i_omino)) - 1, box_size) + 1
                    iy_cell = modulo(current_positions(i_omino,2) + &
                        & rot_omino_shapes(i_cell,2,i_omino_type,current_orientations(i_omino)) - 1, box_size) + 1

                    current_system_state(&
                        & ix_cell, &
                        & iy_cell &
                        & )=i_omino_type
                    
                end do
            end do

            system_snapshots(:,:,i_step/save_every+1) = current_system_state
            history_positions(:,:,i_step/save_every+1)=current_positions
            history_orientations(:,i_step/save_every+1)=current_orientations

        end if



        !---------------- TRANSLATION ----------------
        call generate_randi(1,n_ominos,rand_omino)
        call generate_randi(1,box_size,rand_x)
        call generate_randi(1,box_size,rand_y)

        call check_overlap_after_translation(rand_omino, rand_x, rand_y, has_overlap)

        ! Only accept the move if there's no overlap
        if (.not. has_overlap) then
            current_positions(rand_omino,1) = rand_x
            current_positions(rand_omino,2) = rand_y
        end if


        !---------------- ROTATION ----------------
        call generate_randi(1,n_ominos,rand_omino)
        call generate_randi(1,4,rand_rotation)

        call check_overlap_after_rotation(rand_omino, rand_rotation, has_overlap)

        if (.not. has_overlap) then
            current_orientations(rand_omino) = rand_rotation
        end if

    end do

    

    contains


    subroutine generate_randi(min,max,randi)
        implicit none
        integer, intent(in) :: min, max
        integer, intent(out) :: randi
        double precision :: rand_real

        call random_number(rand_real)
        randi = min + int(rand_real * (max - min + 1)) ! [[min, max]] inclusive
        
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
        
        do i_omino = 1, n_ominos
            if (i_omino == omino_id) cycle

            delta_x = current_positions(i_omino,1) - new_x
            if (delta_x > box_size/2) then
                delta_x = delta_x - box_size
            else if (delta_x < -box_size/2) then
                delta_x = delta_x + box_size
            end if

            if (delta_x < 2*max_omino_size) then 
                delta_y = current_positions(i_omino,2) - new_y
                if (delta_y > box_size/2) then
                    delta_y = delta_y - box_size
                else if (delta_y < -box_size/2) then
                    delta_y = delta_y + box_size
                end if

                if (delta_y < 2*max_omino_size) then 
                    
                    other_omino_type = vec_omino_types(i_omino)
                    
                    do i_cell = 1, ominos_sizes(i_omino_type)
                        ! Use rotated shape for the moved omino
                        ix_cell = modulo(new_x + & 
                            & rot_omino_shapes(i_cell,1,i_omino_type,current_orientations(omino_id)) - 1, box_size) + 1
                        iy_cell = modulo(new_y + &
                            & rot_omino_shapes(i_cell,2,i_omino_type,current_orientations(omino_id)) - 1, box_size) + 1
                        
                        do j_cell = 1, ominos_sizes(other_omino_type)
                            ! Use rotated shape for the other omino
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
        
        do i_omino = 1, n_ominos
            if (i_omino == omino_id) cycle
            
            other_omino_type = vec_omino_types(i_omino)

            delta_x = current_positions(i_omino,1) - current_positions(omino_id,1)
            if (delta_x > box_size/2) then
                delta_x = delta_x - box_size
            else if (delta_x < -box_size/2) then
                delta_x = delta_x + box_size
            end if

            if (delta_x < 2*max_omino_size) then 
                delta_y = current_positions(i_omino,2) - current_positions(omino_id,2)
                if (delta_y > box_size/2) then
                    delta_y = delta_y - box_size
                else if (delta_y < -box_size/2) then
                    delta_y = delta_y + box_size
                end if

                if (delta_y < 2*max_omino_size) then 
            
                    do i_cell = 1, ominos_sizes(i_omino_type)
                        ! Use the NEW rotation for the omino being rotated
                        ix_cell = modulo(current_positions(omino_id,1) + rot_omino_shapes(i_cell,1,i_omino_type,new_rotation) - 1, box_size) + 1
                        iy_cell = modulo(current_positions(omino_id,2) + rot_omino_shapes(i_cell,2,i_omino_type,new_rotation) - 1, box_size) + 1
                        
                        do j_cell = 1, ominos_sizes(other_omino_type)
                            ! Use current rotation for other ominos
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


end subroutine polyomino_simulation



