
program test

    use libflit
    
    real, allocatable, dimension(:, :) :: w
	
	call mpistart
	
	print *, num2str(rankid), ' of ', num2str(nrank)
	
	call mpibarrier
	
	ngroup = 5
	call mpistart_group
	call mpibarrier_group
	
	print *, num2str(rankid_group), ' of ', num2str(nrank_group), ' @ group = ', num2str(groupid), ', global id = ', num2str(rankid), ' of ', num2str(nrank)
	
	if (rankid_group == 0) then
		w = ones(4, 4)
	else
		w = zeros(4, 4)
	end if
	
	print *, '@group id = ', rankid_group, ' with max = ', maxval(w)
	
	call allreduce_array_group(w)
	call mpibarrier_group
	
	print *, '@group id = ', rankid_group, ' after reduction with max = ', maxval(w)
	
	call mpibarrier_group
	call mpiend_group

	call mpiend

end program test
