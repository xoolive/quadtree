set (ENV{LIBGL_ALWAYS_INDIRECT} "1")

execute_process (
    COMMAND ./test_simu
    DEPENDS test_simu
    RESULT_VARIABLE rv
    )

message ("${rv}")

