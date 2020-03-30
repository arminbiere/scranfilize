Scramble CNFs in DIMACS format.
==================================================================

To build `scranfilize` use

  `./configure && make`

or include testing

  `./configure && make test`

which produces scrambled versions of the CNFs in 'cnfs' in 'log'.

To understand what `scranfilize` can do run

  `./scranfilize -h`

after building it.

This tool is described in our [POS'18 paper](http://fmv.jku.at/papers/BiereHeule-POS18.pdf)  ([bibtex](http://fmv.jku.at/papers/BiereHeule-POS18.bib)):

  Armin Biere, Marijn Heule.
  [The Effect of Scrambling CNFs](http://fmv.jku.at/papers/BiereHeule-POS18.pdf).
  In Proceedings 9th Workshop on Pragmatics of SAT 2015 and 2018,
  EPiC Series in Computing, vol. 59, pages 111-126, EasyChair, 2019.

The set of [experiments](http://fmv.jku.at/scranfilize) are described at <http://fmv.jku.at/scranfilize>.
