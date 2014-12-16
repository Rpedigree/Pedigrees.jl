# Pedigrees - pedigree representation and extractors

[![Build Status](https://travis-ci.org/Rpedigree/pedigree.jl.svg?branch=master)](https://travis-ci.org/Rpedigree/pedigree.jl)

## Pedigree representation

A pedigree is represented as a [Julia](http://www.julialang.org) type
constructed from two integer vectors, `sire` and `dam`.  A parent not
in the pedigree is coded as `0`.

By default the `Pedigree`constructor reorders the pedigree according
to the "longest ancestral path" for each animal.  This guarantees that
parents in the pedigree occur before their offspring, because the
longest ancestral path of the offspring is greater than that of either
parent.  The permutation is stored in the `Pedigree` object.

## Additive relationship matrix

An important tool in fitting models that take into account genetic
relationships is the additive relationship matrix, usually written
`A`, for a pedigree.  (See R.A. Mrode, _Linear Models for the
Prediction of Animal Breeding Values, 2nd edition_, 2004, chapter 2).

In practice, the Cholesky factor, `L`, which is a lower triangular
matrix with positive diagonal elements and satisfying $A = LL'$ is
more useful.  `L` or its transpose, often written `R` because $A=R'R$
so `R` is the one on the right and `L` is the one on the left, can be
evaluated directly.




