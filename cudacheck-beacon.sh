#!/bin/bash

module load cuda/10.0
exec ${CUDA_DIR}/extras/demo_suite/deviceQuery
