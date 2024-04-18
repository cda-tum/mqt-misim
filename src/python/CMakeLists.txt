if(NOT TARGET pymisim)
  pybind11_add_module(
    _misim
    # Prefer thin LTO if available
    THIN_LTO
    # Optimize the bindings for size
    OPT_SIZE
    # Source code goes here
    bindings.cpp)
  target_link_libraries(_misim PRIVATE MQT::MiSiM MQT::ProjectOptions MQT::ProjectWarnings)

  # Install directive for scikit-build-core
  install(
    TARGETS _misim
    DESTINATION .
    COMPONENT mqt-misim_PythonModule)
endif()
