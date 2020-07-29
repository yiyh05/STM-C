	subroutine compute_dzmm(zmm, lx, dzmm)
	implicit none
	integer i, lx
	double precision zmm(-5:23), dzmm(23)

	dzmm(1) = zmm(2)/2.0d0
	dzmm(2) = (zmm(3) + zmm(2))/2.0d0 - zmm(2)/2.0d0
	do i = 3, lx
         if (i /= lx) then
            dzmm(i) = (zmm(i+1) - zmm(i-1))/2.0d0 ! thickness of first layer
         else if (i == lx) then  
            dzmm(i) = (zmm(i)-zmm(i-1))/2.0d0
         endif
	enddo
	
	end subroutine compute_dzmm
