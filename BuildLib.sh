#!/usr/bin/env bash
mkdir Build
(cd Build && cmake -DLIB_BUILD=ON ../ && make)