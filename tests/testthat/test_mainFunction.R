########## Test That ##########
## main function ##


context("Description")  # Describe the file, grouping multiple tests



test_that("str_length is number of characters", { # A test, grouping multiple expectations
  expect_equal(length("a"), 1)
})



# Expectations

#expect_equal()
#expect_identical()


#expect_match()
string="TestinG"
expect_match(string, "testing", ignore.case = TRUE)

#expect_output()
#expect_message()
#expect_warning()
#expect_error()
# With expect_message(), expect_warning(), expect_error() you can leave the second argument blank if you just want to see if a message, warning or error is created.
# However, it’s normally better to be explicit, and provide some text from the message.

#expect_is()


#expect_true()
#expect_false()
