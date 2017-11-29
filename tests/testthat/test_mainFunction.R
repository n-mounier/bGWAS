########## Test That ##########
## main function ##


context("Test of main function, bGWAS()")  # Describe the file, grouping multiple tests



test_that("Handling bad parameters", { # A test, grouping multiple expectations
  # Require mandatory parameters

  # Name should be a string
  expect_error(bGWAS(Name=123, GWAS="Test"), "Name : non-character argument", fixed=T)
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
