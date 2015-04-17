function wavg_3vars_rank3_realv_realw(w1, v1, w2, v2, w3, v3) result(wavg)
    !-------------------------------------------------------------------
    ! Weighted average of 3 variables each a real of rank 3 with real 
    ! weights.
    !-------------------------------------------------------------------
    implicit none
    real, dimension(:, :, :), intent(in) :: v1, v2, v3
    real, intent(in) :: w1, w2, w3
    real, dimension(shape(v1)) :: wavg
    wavg = w1 * v1 + w2 * v2 + w3 * v3
end function wavg_3vars_rank3_realv_realw

