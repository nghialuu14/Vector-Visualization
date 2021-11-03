cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(XMLStructuredGridWriter)

find_package(VTK COMPONENTS 
  vtkCommonCore
  vtkCommonDataModel
  vtkIOXML
  vtkIOLegacy
  QUIET
)

if (NOT VTK_FOUND)
  message("Skipping XMLStructuredGridWriter: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(XMLStructuredGridWriter MACOSX_BUNDLE XMLStructuredGridWriter.cxx )
  target_link_libraries(XMLStructuredGridWriter PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(XMLStructuredGridWriter MACOSX_BUNDLE XMLStructuredGridWriter.cxx )
  target_link_libraries(XMLStructuredGridWriter PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS XMLStructuredGridWriter
    MODULES ${VTK_LIBRARIES}
    )
endif ()