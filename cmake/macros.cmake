FUNCTION(ADD_NO_COMMENT_TEST TEST_NAME)
    ADD_EXECUTABLE(${TEST_NAME}.exe ${PROJECT_SOURCE_DIR}/test/${TEST_NAME}.cc)
    TARGET_LINK_LIBRARIES(${TEST_NAME}.exe no_comment)
    TARGET_LINK_LIBRARIES(${TEST_NAME}.exe deal_II)
    TARGET_LINK_LIBRARIES(${TEST_NAME}.exe ${DataTransferKit_LIBRARIES})
    TARGET_LINK_LIBRARIES(${TEST_NAME}.exe ${Boost_LIBRARIES})
    ADD_TEST(NAME ${TEST_NAME} COMMAND ${TEST_NAME}.exe)
ENDFUNCTION()
