cmake_minimum_required(VERSION 3.0)
project(libpyquante2 C)

add_library(pyquante2 STATIC chgp.c cints.c)

install(TARGETS pyquante2 DESTINATION .)
