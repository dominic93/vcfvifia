add_executable (tests
  tests.cpp
  FilterTest.cpp
  UtilsTest.cpp
  BasicStatisticTest.cpp
  ValidatorTest.cpp
  ../src/Plotter.cpp
  ../src/BasicStatistic.cpp
  ../src/Filter.cpp
  ../src/utils.cpp
  ../src/Validator.cpp
  )
  target_link_libraries(tests ${VCFLIB_PATH}/libvcflib.a ${VCFLIB_PATH}/tabixpp/htslib/libhts.a ${Boost_LIBRARIES} pthread z m)
  add_test(NAME all-unittests COMMAND tests)
