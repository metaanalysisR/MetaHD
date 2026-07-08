test_that("NULL Y colnames assigned from Slist", {
  d <- make_test_data()
  Y_nonames <- d$Y
  colnames(Y_nonames) <- NULL
  
  expect_no_error(
    model <- MetaHD(Y = Y_nonames,
                    Slist = d$Slist,
                    method = "multi")
  )
})

test_that("NULL Slist colnames assigned from Y", {
  d <- make_test_data()
  Slist_nonames <- lapply(d$Slist, function(Sk) {
    colnames(Sk) <- rownames(Sk) <- NULL
    Sk
  })
  
  expect_no_error(
    MetaHD(Y = d$Y,
           Slist = Slist_nonames,
           method = "multi")
  )
})

test_that("both NULL assigned generic names", {
  d <- make_test_data()
  Y_nonames <- d$Y
  colnames(Y_nonames) <- NULL
  Slist_nonames <- lapply(d$Slist, function(Sk) {
    colnames(Sk) <- rownames(Sk) <- NULL
    Sk
  })
  
  expect_no_error(
    MetaHD(Y = Y_nonames,
           Slist = Slist_nonames,
           method = "multi")
  )
})

test_that("mismatched Y and Slist names throw error", {
  d <- make_test_data()
  Y_wrongnames <- d$Y
  colnames(Y_wrongnames) <- paste0("X", seq_len(d$N))  # different names
  
  expect_error(
    MetaHD(Y = Y_wrongnames,
           Slist = d$Slist,
           method = "multi"),
    "Column names of Y and Slist do not match"
  )
})