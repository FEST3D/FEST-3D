#define Fatal_error write(*,'(A,A,A,I0,A,I0)')"* ERROR *"//Achar(10)//"In file: ", __FILE__,Achar(10)//"At line number: ",__LINE__, Achar(10)//"For process number: ", process_id; STOP
#define Warning write(*,'(A,A,A,I0,A,I0)')"* WARNING *"//Achar(10)//"In file: ", __FILE__,Achar(10)//"At line number: ",__LINE__, Achar(10)//"For process number: ", process_id
#define Error_msg "* ERROR *"//Achar(10)//"In file: "// __FILE__//Achar(10)//"At line number: "+__LINE__
