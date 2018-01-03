add_executable(gh-issue-1 gh-issue-1.cpp)
target_link_libraries(gh-issue-1 ${PROJECT_NAME})
add_test(NAME gh-issue-1 COMMAND gh-issue-1)

add_executable(gh-issue-2 gh-issue-2.cpp)
target_link_libraries(gh-issue-2 ${PROJECT_NAME})
add_test(NAME gh-issue-2 COMMAND gh-issue-2)

find_package(Armadillo)
include_directories(${ARMADILLO_INCLUDE_DIRS})
add_executable(testWigner testWigner.cpp)
target_link_libraries(testWigner ${PROJECT_NAME} ${ARMADILLO_LIBRARIES})
add_test(NAME testWigner COMMAND testWigner 10)
