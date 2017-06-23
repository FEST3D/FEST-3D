module str_case
  use global, only: STRING_BUFFER_LENGTH
  implicit none
  private
  
  character(len=STRING_BUFFER_LENGTH) :: res
  public :: ucase
  public :: lcase
  
  contains

    function ucase(text) result(res)
      CHARACTER(len=*), intent(in)        :: text
      character(len=STRING_BUFFER_LENGTH) :: res
      integer ::  I,C
  
      res=text
      DO I = 1,LEN(TEXT)
        C = INDEX("abcdefghijklmnopqrstuvwxyz",TEXT(I:I))
        IF (C.GT.0) res(I:I) = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"(C:C)
      END DO
  
    end function ucase

    function lcase(text) result(res)
      CHARACTER(len=*), intent(in)         :: text
      character(len=STRING_BUFFER_LENGTH) :: res
      integer ::  I,C
  
      res=text
      DO I = 1,LEN(TEXT)
        C = INDEX("ABCDEFGHIJKLMNOPQRSTUVWXYZ",TEXT(I:I))
        IF (C.GT.0) res(I:I) = "abcdefghijklmnopqrstuvwxyz"(C:C)
      END DO
  
    end function lcase

end module str_case
