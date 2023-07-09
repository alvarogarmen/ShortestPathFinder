# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "C:/Users/alvar/CLionProjects/praktikum/cmake-build-debug/sparsehash/src/sparsehash_external"
  "C:/Users/alvar/CLionProjects/praktikum/cmake-build-debug/sparsehash/src/sparsehash_external-build"
  "C:/Users/alvar/CLionProjects/praktikum/cmake-build-debug/sparsehash"
  "C:/Users/alvar/CLionProjects/praktikum/cmake-build-debug/sparsehash/tmp"
  "C:/Users/alvar/CLionProjects/praktikum/cmake-build-debug/sparsehash/src/sparsehash_external-stamp"
  "C:/Users/alvar/CLionProjects/praktikum/cmake-build-debug/sparsehash/src"
  "C:/Users/alvar/CLionProjects/praktikum/cmake-build-debug/sparsehash/src/sparsehash_external-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "C:/Users/alvar/CLionProjects/praktikum/cmake-build-debug/sparsehash/src/sparsehash_external-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "C:/Users/alvar/CLionProjects/praktikum/cmake-build-debug/sparsehash/src/sparsehash_external-stamp${cfgdir}") # cfgdir has leading slash
endif()
