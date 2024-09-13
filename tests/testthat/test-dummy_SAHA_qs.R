test_that("ann data can be created.", {
   data("x")
   data("y")
   ann=Create_SAHA_object(query = x,db = y,data_type = "Markers")
   expect_equal(typeof(ann),"S4")
})
