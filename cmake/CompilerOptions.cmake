include_guard(GLOBAL)

function(gmd_apply_default_warnings target_name)
    if(MSVC)
        target_compile_options(${target_name} INTERFACE /W4 /permissive-)
    else()
        target_compile_options(${target_name} INTERFACE -Wall -Wextra -Wpedantic)
    endif()
endfunction()
