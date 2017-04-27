module math_utils

    implicit none
    private
    
    public :: wavg
    interface wavg
        !---------------------------------------------------------------
        ! wavg (Weighted average) is a set of overloaded functions that 
        ! can be used to find the weighted average of the passed 
        ! arguments. The general syntax of invoking this function is:
        !   wavg(w1, val1, w2, val2, w3, val3, ...)
        !---------------------------------------------------------------
        module procedure wavg_3vars_rank3_realv_realw
    end interface wavg

    contains
        
        include "wavg_implementation.inc"

end module math_utils
