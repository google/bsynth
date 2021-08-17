#!/bin/bash

docker build --tag bsynth .

docker run --rm \
  -v $(pwd):/tmp/working_dir \
  -v $(pwd)/../..:/tmp/bsynth \
  -w /tmp/working_dir \
  bsynth \
  R -e "devtools::build(pkg = '/tmp/bsynth', path = '/tmp/working_dir, binary = TRUE')"


