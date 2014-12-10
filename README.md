# pedigree - pedigree representation and extractors

[![Build Status](https://travis-ci.org/Rpedigree/pedigree.jl.svg?branch=master)](https://travis-ci.org/Rpedigree/pedigree.jl)

## Pedigree representation

A pedigree is represented as a [Julia](http://www.julialang.org) type with two integer
members, `sire` and `dam`.  A parent not in the pedigree is coded as `0`.

The constructor checks for an ordered pedigree (parents occur before children) and re-orders
the members if they are not already ordered.  The longest ancestral path for each animal is
evaluated and stored in the pedigree.
