  !< Change the full string to particular case: upper or lower.
module str_case
  !< Change the full string to particular case: upper or lower.
  use global, only: STRING_BUFFER_LENGTH
  implicit none
  private
  
  character(len=STRING_BUFFER_LENGTH) :: res
  public :: ucase
  public :: lcase
  
  contains

    function ucase(text) result(res)
      !<Mmake the whole string to upper case
      CHARACTER(len=*), intent(in)        :: text
      !< Input string of any case
      character(len=STRING_BUFFER_LENGTH) :: res
      !< Output string of upper case
      integer ::  I,C
  
      res=text
      DO I = 1,LEN(TEXT)
        C = INDEX("abcdefghijklmnopqrstuvwxyz",TEXT(I:I))
        IF (C.GT.0) res(I:I) = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"(C:C)
      END DO
  
    end function ucase

    function lcase(text) result(res)
      !< Make the whole string to lower case
      CHARACTER(len=*), intent(in)         :: text
      !< Input string of any case
      character(len=STRING_BUFFER_LENGTH) :: res
      !< Output string of lower case
      integer ::  I,C
  
      res=text
      DO I = 1,LEN(TEXT)
        C = INDEX("ABCDEFGHIJKLMNOPQRSTUVWXYZ",TEXT(I:I))
        IF (C.GT.0) res(I:I) = "abcdefghijklmnopqrstuvwxyz"(C:C)
      END DO
  
    end function lcase

end module str_case
