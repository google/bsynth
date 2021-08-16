## Copyright 2021 Google LLC
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     https://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
##
test_that("the factory works", {
  # Create a new factory
  germany_synth <- bayesianSynth$new(
    data = germany,
    time = year,
    id = country,
    treated = D,
    outcome = gdp,
    ci_width = 0.95
  )
  expect_is(germany_synth, "bayesianSynth")
})
