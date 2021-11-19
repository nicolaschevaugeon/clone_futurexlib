# Pull up source list to the AOMD target defined in the upper directory
# Note: needs Cmake 3.13
target_sources(AOMD
    PRIVATE
# Headers   
include/AOMD_LoadBalancer.h
include/AOMD_SharedInfo.h
include/AOMD_ZoltanLoadBalancer.h
include/AOMDfwd.h
include/MeshModel.h
include/mAttachableDataContainer.h
include/mBuildAdj.h
include/mCommonBdry.h
include/mDebugUtil.h
include/mEdge.h
include/mEntityContainer.h
include/mException.h
include/mExchangeData.h
include/mFace.h
include/mGeomSearch.h
include/mHex.h
include/mIterator.h
include/mIteratorNew.h
include/mMirrorEntity.h
include/mPrism.h
include/mPyramid.h
include/mRegion.h
include/mTet.h
include/mUnits.h
include/mVertex.h
include/pmEntity.h
include/AOMD_METIS_Partitioner.h
include/AOMD_OwnerManager.h
include/ParUtil.h
include/mAOMD.h
include/mEntity.h
include/mExchangeDataGeneric.h
include/mIdGenerator.h
include/mMesh.h
include/mMigrateUtil.h
include/mParallelUtility.h
include/pmGraphs.h
cint/AOMD.h
cint/AOMDInternals.h
cint/AOMDList.h
cint/AOMD_Defs.h
cint/AOMD_Internals.h
cint/MSops.h
cint/AOMD_cint.h
# Sources
src/MeshModel.cc
src/mBuildGTopology.cc
src/mCommonBdry.cc
src/mEdge.cc
src/mEntityContainer.cc
src/mException.cc
src/mFace.cc
src/mGambit.cc
src/mHex.cc
src/mIterator.cc
src/mIteratorNew.cc
src/mMirrorEntity.cc
src/mOneLevel.cc
src/mPrism.cc
src/mRegion.cc
src/mSMS_oneLevel.cc
src/mTet.cc
src/mVertex.cc
src/mSMS.cc
src/mAMR.cc
src/mAOMD.cc
src/mDebugUtil.cc
src/mEntity.cc
src/mGraphs.cc
src/mIdGenerator.cc
src/mMesh.cc
cint/AOMD_Exist.cc
cint/AOMD_Mod.cc
cint/AOMD_FlexDB.cc
cint/AOMD.cc
cint/AOMD_Entity.cc
cint/AOMD_Internals.cc
cint/AOMD_Mesh.cc
cint/PAOMD_Internals.cc
cint/PARALLEL.cc
pmodel/pmEntity.cc
pmodel/ParUtil.cc
pmodel/pmMeshEntity.cc
pmodel/pmModelUtil.cc
parallel/AOMDBdryLinkSetup.cc
parallel/AOMD_OwnerManager.cc
parallel/mExchangeData.cc
parallel/mMigrateUtil.cc
)

# Pull up include directories to the AOMD target defined in the upper directory
target_include_directories(AOMD PUBLIC include cint)



