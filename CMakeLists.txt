cmake_minimum_required(VERSION 3.5)
project(version_2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lpthread -pthread")

set(SOURCE_FILES
    CMakeLists.txt
    coordinaryTransform.cpp
    fileDisposal.cpp
    initialization.cpp
    main.cpp
    matrixFunction.cpp
    MOPSO.h
    MOPSOAidFunction.cpp
    MOPSOFunction.cpp
    parameterSettings.cpp
        setTime.cpp
        preDispose.cpp
        multiplyBetterInput.cpp
        disposedCodes.cpp
        seiveRep.cpp
        inputParticles.cpp
        createRep.cpp
        checkSimilarity.cpp
        Debug.cpp
        freeUpSpace.cpp
        applyVariables.cpp
        Matrix.cpp
        randFixedSum.cpp disposeTM_align.cpp)

add_executable(version_2 ${SOURCE_FILES} setTime.cpp preDispose.cpp multiplyBetterInput.cpp disposedCodes.cpp seiveRep.cpp inputParticles.cpp Matrix.cpp randFixedSum.cpp)