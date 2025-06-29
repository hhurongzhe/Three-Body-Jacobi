add_rules("mode.release")
add_requires("openmp")


target("apwd3-c1")
    set_kind("binary")
    set_languages("c17", "c++17")
    set_optimize("fastest")
    add_packages("openmp")
    add_files("src/numlib/*.cpp")
    add_files("src/aPWD3_part_c1.cpp")
    add_files("src/main_part_c1.cpp")
    if is_host("windows") then
        set_filename("apwd3-c1.exe")
    elseif is_plat("linux") then
        set_filename("apwd3-c1.x")
    elseif is_plat("macosx") then
        add_includedirs("/opt/homebrew/Cellar/libomp/17.0.6/include")
        add_linkdirs("/opt/homebrew/Cellar/libomp/17.0.6/lib")
        add_links("omp")
        set_filename("apwd3-c1")
    end


target("apwd3-c3")
    set_kind("binary")
    set_languages("c17", "c++17")
    set_optimize("fastest")
    add_packages("openmp")
    add_files("src/numlib/*.cpp")
    add_files("src/aPWD3_part_c3.cpp")
    add_files("src/main_part_c3.cpp")
    if is_host("windows") then
        set_filename("apwd3-c3.exe")
    elseif is_plat("linux") then
        set_filename("apwd3-c3.x")
    elseif is_plat("macosx") then
        add_includedirs("/opt/homebrew/Cellar/libomp/17.0.6/include")
        add_linkdirs("/opt/homebrew/Cellar/libomp/17.0.6/lib")
        add_links("omp")
        set_filename("apwd3-c3")
    end


target("apwd3-c4")
    set_kind("binary")
    set_languages("c17", "c++17")
    set_optimize("fastest")
    add_packages("openmp")
    add_files("src/numlib/*.cpp")
    add_files("src/aPWD3_part_c4.cpp")
    add_files("src/main_part_c4.cpp")
    if is_host("windows") then
        set_filename("apwd3-c4.exe")
    elseif is_plat("linux") then
        set_filename("apwd3-c4.x")
    elseif is_plat("macosx") then
        add_includedirs("/opt/homebrew/Cellar/libomp/17.0.6/include")
        add_linkdirs("/opt/homebrew/Cellar/libomp/17.0.6/lib")
        add_links("omp")
        set_filename("apwd3-c4")
    end


target("apwd3-cD")
    set_kind("binary")
    set_languages("c17", "c++17")
    set_optimize("fastest")
    add_packages("openmp")
    add_files("src/numlib/*.cpp")
    add_files("src/aPWD3_part_cD.cpp")
    add_files("src/main_part_cD.cpp")
    if is_host("windows") then
        set_filename("apwd3-cD.exe")
    elseif is_plat("linux") then
        set_filename("apwd3-cD.x")
    elseif is_plat("macosx") then
        add_includedirs("/opt/homebrew/Cellar/libomp/17.0.6/include")
        add_linkdirs("/opt/homebrew/Cellar/libomp/17.0.6/lib")
        add_links("omp")
        set_filename("apwd3-cD")
    end


target("apwd3-cE")
    set_kind("binary")
    set_languages("c17", "c++17")
    set_optimize("fastest")
    add_packages("openmp")
    add_files("src/numlib/*.cpp")
    add_files("src/aPWD3_part_cE.cpp")
    add_files("src/main_part_cE.cpp")
    if is_host("windows") then
        set_filename("apwd3-cE.exe")
    elseif is_plat("linux") then
        set_filename("apwd3-cE.x")
    elseif is_plat("macosx") then
        add_includedirs("/opt/homebrew/Cellar/libomp/17.0.6/include")
        add_linkdirs("/opt/homebrew/Cellar/libomp/17.0.6/lib")
        add_links("omp")
        set_filename("apwd3-cE")
    end


target("benchmark")
    set_kind("binary")
    set_languages("c17", "c++17")
    set_optimize("fastest")
    add_packages("openmp")
    add_files("src/numlib/*.cpp")
    add_files("src/aPWD3_part_c1.cpp","src/aPWD3_part_c3.cpp","src/aPWD3_part_c4.cpp")
    add_files("src/main_benchmark.cpp")
    if is_host("windows") then
        set_filename("apwd3-benchmark.exe")
    elseif is_plat("linux") then
        set_filename("apwd3-benchmark.x")
    elseif is_plat("macosx") then
        add_includedirs("/opt/homebrew/Cellar/libomp/17.0.6/include")
        add_linkdirs("/opt/homebrew/Cellar/libomp/17.0.6/lib")
        add_links("omp")
        set_filename("apwd3-benchmark")
    end
