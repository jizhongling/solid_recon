cmake_minimum_required(VERSION 3.16)

add_subdirectory(interfaces)
add_subdirectory(calorimetry)
add_subdirectory(tracking)
add_subdirectory(digi)
add_subdirectory(reco)

if(USE_ONNX)
  add_subdirectory(onnx)
endif()
