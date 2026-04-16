include_guard(GLOBAL)

function(gmd_configure_cuda_target target_name)
    target_compile_features(${target_name} PUBLIC cuda_std_17)
endfunction()
