context("Test that the safedata network comms fail gracefully")
library(safedata)

test_that("no internet fails gracefully", {
    Sys.setenv(NETWORK_DOWN = TRUE)

    success <- safedata:::try_to_download("https://httpbin.org")

    expect_false(success)
    expect_match(attr(success, "fail_msg"), regexp = "No internet connection")

    Sys.unsetenv("NETWORK_DOWN")
})

# All of the tests below rely on having the internet available to generate
# the various error messages from safedata:::try_to_download(), so use a
# check function to skip if there is a network outage.

internet_unavailable <- function() {
    return(!curl::has_internet())
}

test_that("bad host fails gracefully", {
    if (internet_unavailable()) {
        skip("No internet - skipping test")
    }

    success <- safedata:::try_to_download("https://httpbinzzzzz.org")

    expect_false(success)
    expect_match(
        attr(success, "fail_msg"),
        regexp = "URL not found|URL timed out"
    )
})

test_that("timeout fails gracefully", {
    if (internet_unavailable()) {
        skip("No internet - skipping test")
    }

    success <- safedata:::try_to_download(
        "https://httpbin.org/delay/2",
        timeout = 1
    )

    # If httpbin genuinely times out, the end result is the same.
    expect_false(success)
    expect_match(attr(success, "fail_msg"), regexp = "URL timed out")
})

test_that("URL errors fails gracefully", {
    if (internet_unavailable()) {
        skip("No internet - skipping test")
    }

    success <- safedata:::try_to_download("https://httpbin.org/status/404")
    expect_false(success)

    # If httpbin issues a Gateway timeout 504, then "URL error" is still
    # reported, but if it _properly_ times out then we get "URL timed out"
    # and that will cause a test failure unless trapped.
    if (!grepl("URL timed out", attr(success, "fail_msg"))) {
        expect_match(attr(success, "fail_msg"), regexp = "URL error")
    }
})

test_that("Good URL works and returns object to memory", {
    if (internet_unavailable()) {
        skip("No internet - skipping test")
    }

    success <- safedata:::try_to_download("https://httpbin.org/base64/c2FmZWRhdGE=")

    # Screen for failure of the download due to timeout
    if (is.logical(success) && isFALSE(success)) {
        skip("Download failed")
    }

    # Check we get content not logical
    expect_false(is.logical(success) && isTRUE(success))
    # Check the content is as expected
    expect_equal(rawToChar(success$content), "safedata")
})

test_that("Bad path fails gracefully", {
    if (internet_unavailable()) {
        skip("No internet - skipping test")
    }

    success <- safedata:::try_to_download(
        "https://httpbin.org/base64/c2FmZWRhdGE=",
        local_path = "/does/not/exist"
    )

    expect_false(success)
    expect_match(attr(success, "fail_msg"), regexp = "Failed to open file")
})

test_that("no internet and safedata dir creation fails gracefully", {
    Sys.setenv(NETWORK_DOWN = TRUE)
    temp_safe_dir <- file.path(tempdir(), "test_safe_data")

    success <- expect_message(
        safedata::set_safe_dir(temp_safe_dir, create = TRUE),
        regexp = "Could not download required files: SAFE data directory not created",
    )
    expect_false(success)

    # Check it tidied up
    expect_false(file.exists(temp_safe_dir))
    Sys.unsetenv("NETWORK_DOWN")
})

test_that("API down and safedata dir creation fails gracefully", {
    Sys.setenv(URL_DOWN = TRUE)
    temp_safe_dir <- file.path(tempdir(), "test_safe_data")

    success <- expect_message(
        safedata::set_safe_dir(temp_safe_dir, create = TRUE),
        regexp = "Could not download required files: SAFE data directory not created",
    )

    expect_false(success)
    # Check it tidied up
    expect_false(file.exists(temp_safe_dir))
    Sys.unsetenv("URL_DOWN")
})
