#!/usr/bin/env bash
mkdir Build
(cd Build && cmake -DMATHLINK_BUILD=ON ../ && make)