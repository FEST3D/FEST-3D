#define Release

#ifdef Debug
#define DebugCall(arg) if(process_id==0) write(*,'(A,I0,A)') "Call_one (In "//__FILE__//", at ",__LINE__,", "//arg//" )"
#define DebugInfo(arg) if(process_id==0) write(*,'(A,I0,A)') "Info_one (In "//__FILE__//", at ",__LINE__,", "//arg//" )"
#endif

#ifdef Debug_All
#define DebugCall(arg) write(*,'(A,I0,A)') "Call_all (In "//__FILE__//", at ",__LINE__,", "//arg//" )"
#define DebugInfo(arg) write(*,'(A,I0,A)') "Info_all (In "//__FILE__//", at ",__LINE__,", "//arg//" )"
#endif

#ifdef Release
#define DebugCall(arg) !should be nothing
#define DebugInfo(arg) !should be nothing
#endif
