# 显式指定 Clang 编译器
set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_C_COMPILER "clang")

cmake_minimum_required(VERSION 3.10)
project(rank-fhe)

# 启用调试模式
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找OpenFHE库
set(OpenFHE_DIR "/home/delld/spz/openfhe-development/build")  # 替换为实际构建目录
find_package(OpenFHE)
if(NOT OpenFHE_FOUND)
    message(FATAL_ERROR "OpenFHE library not found. Please install OpenFHE or set OpenFHE_DIR.")
endif()

# 包含OpenFHE头文件目录
include_directories(
    ${OpenFHE_INCLUDE}
    ${OpenFHE_INCLUDE}/third-party/include
    ${OpenFHE_INCLUDE}/core
    ${OpenFHE_INCLUDE}/pke
    ${OpenFHE_INCLUDE}/binfhe
)

# 添加 OpenFHE 编译标志
set(CMAKE_CXX_FLAGS "${OpenFHE_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

# 强制替换不兼容标志（针对Clang）
string(REPLACE "-fopenmp=libomp" "-fopenmp" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# 处理 OpenMP 依赖
find_package(OpenMP)
if(OpenMP_FOUND)
    set(OPENMP_LINK_FLAGS OpenMP::OpenMP_CXX)
else()
    # Clang 回退方案（如果未找到 OpenMP）
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(OPENMP_LINK_FLAGS "-lgomp")
endif()

# 包含项目头文件
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# 列出所有源文件
set(SOURCES
    main.cpp
    encryption.cpp
    rank.cpp
    sign.cpp
)

# 先定义可执行目标
add_executable(rank-fhe ${SOURCES})

# 后链接库（修复顺序错误）
target_link_libraries(rank-fhe PRIVATE
    ${OpenFHE_SHARED_LIBRARIES}
    ${OPENMP_LINK_FLAGS}
)

