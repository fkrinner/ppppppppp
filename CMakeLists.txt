
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
PROJECT(ppppppppp)
SET(CMAKE_CXX_COMPILER "/opt/rh/devtoolset-7/root/usr/bin/g++")
SET(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR} )

FIND_PACKAGE(Nlopt)
FIND_PACKAGE(ROOT 6.14.202 REQUIRED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Werror -O3 -fopenmp  -Wno-ignored-attributes -Wno-deprecated-declarations -Wno-int-in-bool-context -Wno-misleading-indentation")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Werror -O3")

SET( INCLUDE_DIRS 
    ${INCLUDE_DIRS}
    ${Nlopt_INCLUDE_DIR} 	
    ${PROJECT_BINARY_DIR}  
    ${PROJECT_SOURCE_DIR}  
    $ENV{FABILI_INCLUDE_DIR}
    $ENV{EIGEN_INCLUDE_DIR}
)

INCLUDE_DIRECTORIES( ${INCLUDE_DIRS} )

#set libraries
ADD_LIBRARY(massShape           SHARED massShape.cc           )
ADD_LIBRARY(angularDependence   SHARED angularDependence.cc   )
ADD_LIBRARY(amplitude           SHARED amplitude.cc           )
ADD_LIBRARY(integrator          SHARED integrator.cc          )
ADD_LIBRARY(generator           SHARED generator.cc           )
ADD_LIBRARY(kinematicSignature  SHARED kinematicSignature.cc  )
ADD_LIBRARY(modelAmplitude      SHARED modelAmplitude.cc      )
ADD_LIBRARY(modelGenerator      SHARED modelGenerator.cc      )
ADD_LIBRARY(logLikelihood       SHARED logLikelihood.cc       )
ADD_LIBRARY(efficiencyFunction  SHARED efficiencyFunction.cc  )
ADD_LIBRARY(BELLE_efficiency    SHARED BELLE_efficiency.cc    )
ADD_LIBRARY(getBELLEdata        SHARED getBELLEdata.cc        )
ADD_LIBRARY(fabiliLL            SHARED fabiliLL.cc            )
ADD_LIBRARY(fabiliLL_openMP     SHARED fabiliLL_openMP.cc     )
ADD_LIBRARY(phaseSpaceFunctions SHARED phaseSpaceFunctions.cc )

#set executables
ADD_EXECUTABLE(p9.exe ppppppppp.cxx )
ADD_EXECUTABLE(K9.exe KKKKKKKKK.cxx )
ADD_EXECUTABLE(e9.exe eeeeeeeee.cxx )
ADD_EXECUTABLE(D9.exe DDDDDDDDD.cxx )
ADD_EXECUTABLE(Dnine.exe Dnine.cxx )
ADD_EXECUTABLE(D9_bg.exe DDDDDDDDD_bgFitter.cxx )
ADD_EXECUTABLE(D6.exe DDDDDD.cxx )

ADD_EXECUTABLE(testBelleEfficiecny.exe testBelleEfficiecny.cxx)

# link libraries
TARGET_LINK_LIBRARIES(p9.exe massShape          )
TARGET_LINK_LIBRARIES(p9.exe angularDependence  )
TARGET_LINK_LIBRARIES(p9.exe amplitude          )
TARGET_LINK_LIBRARIES(p9.exe generator          )
TARGET_LINK_LIBRARIES(p9.exe integrator         )
TARGET_LINK_LIBRARIES(p9.exe modelAmplitude     )
TARGET_LINK_LIBRARIES(p9.exe modelGenerator     )
TARGET_LINK_LIBRARIES(p9.exe logLikelihood      )
TARGET_LINK_LIBRARIES(p9.exe efficiencyFunction )
TARGET_LINK_LIBRARIES(p9.exe ${ROOT_LIBS} )

TARGET_LINK_LIBRARIES(K9.exe massShape          )
TARGET_LINK_LIBRARIES(K9.exe angularDependence  )
TARGET_LINK_LIBRARIES(K9.exe amplitude          )
TARGET_LINK_LIBRARIES(K9.exe generator          )
TARGET_LINK_LIBRARIES(K9.exe integrator         )
TARGET_LINK_LIBRARIES(K9.exe modelAmplitude     )
TARGET_LINK_LIBRARIES(K9.exe modelGenerator     )
TARGET_LINK_LIBRARIES(K9.exe logLikelihood      )
TARGET_LINK_LIBRARIES(K9.exe efficiencyFunction )

TARGET_LINK_LIBRARIES(e9.exe massShape          )
TARGET_LINK_LIBRARIES(e9.exe angularDependence  )
TARGET_LINK_LIBRARIES(e9.exe amplitude          )
TARGET_LINK_LIBRARIES(e9.exe generator          )
TARGET_LINK_LIBRARIES(e9.exe integrator         )
TARGET_LINK_LIBRARIES(e9.exe modelAmplitude     )
TARGET_LINK_LIBRARIES(e9.exe modelGenerator     )
TARGET_LINK_LIBRARIES(e9.exe logLikelihood      )
TARGET_LINK_LIBRARIES(e9.exe efficiencyFunction )

TARGET_LINK_LIBRARIES(D9.exe massShape          )
TARGET_LINK_LIBRARIES(D9.exe angularDependence  )
TARGET_LINK_LIBRARIES(D9.exe amplitude          )
TARGET_LINK_LIBRARIES(D9.exe generator          )
TARGET_LINK_LIBRARIES(D9.exe integrator         )
TARGET_LINK_LIBRARIES(D9.exe modelAmplitude     )
TARGET_LINK_LIBRARIES(D9.exe modelGenerator     )
TARGET_LINK_LIBRARIES(D9.exe logLikelihood      )
TARGET_LINK_LIBRARIES(D9.exe efficiencyFunction )
TARGET_LINK_LIBRARIES(D9.exe getBELLEdata       )
TARGET_LINK_LIBRARIES(D9.exe fabiliLL           )
TARGET_LINK_LIBRARIES(D9.exe fabiliLL_openMP    )
TARGET_LINK_LIBRARIES(D9.exe /nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Fabili/build/libfabili.so )

TARGET_LINK_LIBRARIES(Dnine.exe massShape          )
TARGET_LINK_LIBRARIES(Dnine.exe angularDependence  )
TARGET_LINK_LIBRARIES(Dnine.exe amplitude          )
TARGET_LINK_LIBRARIES(Dnine.exe generator          )
TARGET_LINK_LIBRARIES(Dnine.exe integrator         )
TARGET_LINK_LIBRARIES(Dnine.exe modelAmplitude     )
TARGET_LINK_LIBRARIES(Dnine.exe modelGenerator     )
TARGET_LINK_LIBRARIES(Dnine.exe logLikelihood      )
TARGET_LINK_LIBRARIES(Dnine.exe efficiencyFunction )
TARGET_LINK_LIBRARIES(Dnine.exe getBELLEdata       )
TARGET_LINK_LIBRARIES(Dnine.exe fabiliLL           )
TARGET_LINK_LIBRARIES(Dnine.exe fabiliLL_openMP    )
TARGET_LINK_LIBRARIES(Dnine.exe /nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Fabili/build/libfabili.so )

TARGET_LINK_LIBRARIES(D9_bg.exe massShape          )
TARGET_LINK_LIBRARIES(D9_bg.exe angularDependence  )
TARGET_LINK_LIBRARIES(D9_bg.exe amplitude          )
TARGET_LINK_LIBRARIES(D9_bg.exe generator          )
TARGET_LINK_LIBRARIES(D9_bg.exe integrator         )
TARGET_LINK_LIBRARIES(D9_bg.exe modelAmplitude     )
TARGET_LINK_LIBRARIES(D9_bg.exe modelGenerator     )
TARGET_LINK_LIBRARIES(D9_bg.exe logLikelihood      )
TARGET_LINK_LIBRARIES(D9_bg.exe efficiencyFunction )
TARGET_LINK_LIBRARIES(D9_bg.exe getBELLEdata       )
TARGET_LINK_LIBRARIES(D9_bg.exe fabiliLL           )
TARGET_LINK_LIBRARIES(D9_bg.exe fabiliLL_openMP    )
TARGET_LINK_LIBRARIES(D9_bg.exe /nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Fabili/build/libfabili.so )

TARGET_LINK_LIBRARIES(D6.exe massShape          )
TARGET_LINK_LIBRARIES(D6.exe angularDependence  )
TARGET_LINK_LIBRARIES(D6.exe amplitude          )
TARGET_LINK_LIBRARIES(D6.exe generator          )
TARGET_LINK_LIBRARIES(D6.exe integrator         )
TARGET_LINK_LIBRARIES(D6.exe modelAmplitude     )
TARGET_LINK_LIBRARIES(D6.exe modelGenerator     )
TARGET_LINK_LIBRARIES(D6.exe logLikelihood      )
TARGET_LINK_LIBRARIES(D6.exe efficiencyFunction )
TARGET_LINK_LIBRARIES(D6.exe getBELLEdata       )
TARGET_LINK_LIBRARIES(D6.exe fabiliLL           )
TARGET_LINK_LIBRARIES(D6.exe fabiliLL_openMP    )
TARGET_LINK_LIBRARIES(D6.exe /nfs/freenas/tuph/e18/project/compass/analysis/fkrinner/Fabili/build/libfabili.so )

TARGET_LINK_LIBRARIES(efficiencyFunction BELLE_efficiency)

TARGET_LINK_LIBRARIES(amplitude massShape          )
TARGET_LINK_LIBRARIES(amplitude angularDependence  )
TARGET_LINK_LIBRARIES(amplitude kinematicSignature )
TARGET_LINK_LIBRARIES(amplitude kinematicSignature )

TARGET_LINK_LIBRARIES(generator kinematicSignature )

TARGET_LINK_LIBRARIES(angularDependence kinematicSignature )

TARGET_LINK_LIBRARIES(modelAmplitude amplitude )

TARGET_LINK_LIBRARIES(modelGenerator modelAmplitude )
TARGET_LINK_LIBRARIES(modelGenerator generator      )

TARGET_LINK_LIBRARIES(logLikelihood amplitude     )
TARGET_LINK_LIBRARIES(logLikelihood integrator    )
TARGET_LINK_LIBRARIES(logLikelihood ${Nlopt_LIBS} )
#TARGET_LINK_LIBRARIES(logLikelihood ${ROOT_LIBS}  )

TARGET_LINK_LIBRARIES(integrator amplitude          )
TARGET_LINK_LIBRARIES(integrator generator          )
TARGET_LINK_LIBRARIES(integrator efficiencyFunction )

TARGET_LINK_LIBRARIES(getBELLEdata ${ROOT_LIBS}  )

TARGET_LINK_LIBRArIES(testBelleEfficiecny.exe BELLE_efficiency    )
TARGET_LINK_LIBRArIES(testBelleEfficiecny.exe ${ROOT_LIBS}        )
TARGET_LINK_LIBRArIES(testBelleEfficiecny.exe getBELLEdata        )
TARGET_LINK_LIBRArIES(testBelleEfficiecny.exe phaseSpaceFunctions )
TARGET_LINK_LIBRArIES(testBelleEfficiecny.exe massShape           )

TARGET_LINK_LIBRARIES(fabiliLL logLikelihood )

TARGET_LINK_LIBRARIES(fabiliLL_openMP logLikelihood )
