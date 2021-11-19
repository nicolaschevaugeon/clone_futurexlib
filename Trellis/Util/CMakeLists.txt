# Pull up source list to the Util target defined in the upper directory
target_sources(Util
    PRIVATE
# Headers   
include/mPoint.h
include/mTensor2.h
include/mTensor4.h
include/mVector.h
# Sources
src/mTensor2.cc
src/mVector.cc
)

# Pull up include directories to the Util target defined in the upper directory
target_include_directories(Util PUBLIC include)
