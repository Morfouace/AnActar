add_custom_command(OUTPUT MRecoDict.cxx COMMAND ../scripts/build_dict.sh "MTrack.h MHough.h MSimpleRansac.h" MRecoDict.cxx MReco.rootmap libACTARcore.so ACTARRecoLinkDef.h DEPENDS MTrack.h MHough.h MSimpleRansac.h)

add_library(ACTARreco SHARED MTrack.cc MHough.cc MSimpleRansac.cc MRecoDict.cxx)
target_link_libraries(ACTARreco ${ROOT_LIBRARIES})

install(FILES MTrack.h MHough.h MSimpleRansac.h DESTINATION ${CMAKE_BINARY_DIR}/include)
