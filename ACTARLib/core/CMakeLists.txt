add_custom_command(OUTPUT MEventReducedDict.cxx COMMAND ../scripts/build_dict.sh "MEventReduced.h MDetectorConfig.h" MEventReducedDict.cxx MEventReduced.rootmap libACTARcore.so MTreeStructureLinkDef.h DEPENDS MEventReduced.h MDetectorConfig.h)

add_library(ACTARcore SHARED MEventReducedDict.cxx MDetectorConfig.cc)
target_link_libraries(ACTARcore ${ROOT_LIBRARIES})

#install(FILES MChannel.h MCoboAsad.h MEventReduced.h MUnpacker.h DESTINATION ${CMAKE_BINARY_DIR}/include)
install(FILES MEventReduced.h MDetectorConfig.h DESTINATION ${CMAKE_BINARY_DIR}/include)


