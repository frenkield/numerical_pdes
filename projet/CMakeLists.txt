cmake_minimum_required(VERSION 3.15)
project(projet)

set(CMAKE_CXX_STANDARD 14)

# pour umfpack
# mettez le chemin à jour si umfpack.h se trouve ailleurs sur votre système
include_directories(/usr/local/include)

include_directories(include)
link_directories(/usr/local/lib)

file(GLOB hecht_src "src/hecht/*.cpp")
file(GLOB solver_src "src/solver/*.cpp")

set(libraries umfpack cholmod ccolamd colamd camd amd suitesparseconfig)
add_link_options(-framework Accelerate)

# ============================================================================

add_executable(SolveOven src/apps/SolveOven.cpp ${hecht_src} ${solver_src})
target_link_libraries(SolveOven ${libraries})
