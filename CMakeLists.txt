cmake_minimum_required(VERSION 3.21)
project(entzerr)

set(CMAKE_CXX_STANDARD 17)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GTK3 REQUIRED gtk+-3.0)

file(GLOB_RECURSE LIB_SOURCES src/entzerr/*.cpp)

if(LIB_SOURCES)
  add_library(entzerr_lib ${LIB_SOURCES})
  target_include_directories(entzerr_lib PUBLIC include ${GTK3_INCLUDE_DIRS})
  target_link_libraries(entzerr_lib PUBLIC ${GTK3_LIBRARIES})
  target_compile_options(entzerr_lib PRIVATE ${GTK3_CFLAGS_OTHER})
endif()

add_executable(entzerr src/main.cpp)
target_include_directories(entzerr PRIVATE include ${GTK3_INCLUDE_DIRS})
target_link_libraries(entzerr PRIVATE ${GTK3_LIBRARIES})
target_compile_options(entzerr PRIVATE ${GTK3_CFLAGS_OTHER})

if(LIB_SOURCES)
  target_link_libraries(entzerr PRIVATE entzerr_lib)
endif()

include(GNUInstallDirs)

install(TARGETS entzerr
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(FILES LICENSE
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
)

install(FILES res/entzerr.desktop
    DESTINATION ${CMAKE_INSTALL_DATADIR}/applications
)

install(FILES res/entzerr.png
    DESTINATION ${CMAKE_INSTALL_DATADIR}/pixmaps
)
