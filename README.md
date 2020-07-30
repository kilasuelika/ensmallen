This is a fork of Ensmallen using Eigen as backend.

3 main goals:

1. Remove checks.
2. Rewrite traits.
3. Rewrite API calls of matrices.
4. Rewrite tests.

### Requirements
1. Template library Eigen.

### Important Notes
1. Use **VectorXd** or **VectorXf** as type of arguments and gradients. They are column vectors.

### Finished
 - [x] DE
 - [x] L_BFGS
 - [x] GradientDescent
 
### Developers
1. Zhou Yao

### Important changes
1. Traits and SFINAE utilities are removed. Now you are forced to write at most 3 methods for function type if necessary. 
2. Callbacks are removed. If you really want to use it, modify source code directly.

With these changes, Ensmallen-eigen are simply a collection of unrelated optimization algorithms. In fact, you can just copy the single folder of one optimizer without install whole library.


#### Optimizers that only require Evaluate()
1. DE

#### Optimizers that requires Evaluate(), Gradient(), EvaluateWithGradient()
1. L_BFGS
2. GradientDescent


-----


<h2 align="center">
  <a href="http://ensmallen.org/"><img src="http://ensmallen.org/img/ensmallen_text.svg" style="background-color:rgba(0,0,0,0);" height=230 alt="ensmallen: a C++ header-only library for numerical optimization"></a>
</h2>

**ensmallen** is a C++ header-only library for numerical optimization.

Documentation and downloads: http://ensmallen.org

ensmallen provides a simple set of abstractions for writing an objective
function to optimize. It also provides a large set of standard and cutting-edge
optimizers that can be used for virtually any numerical optimization task.
These include full-batch gradient descent techniques, small-batch techniques,
gradient-free optimizers, and constrained optimization.


### Requirements

* C++ compiler with C++11 support
* Armadillo: http://arma.sourceforge.net
* OpenBLAS or Intel MKL or LAPACK (see Armadillo site for details)


### Installation

ensmallen can be installed with CMake 3.3 or later.
If CMake is not already available on your system, it can be obtained from https://cmake.org

If you are using an older system such as RHEL 7 or CentOS 7,
an updated version of CMake is also available via the EPEL repository via the `cmake3` package.


### License

Unless stated otherwise, the source code for **ensmallen** is licensed under the
3-clause BSD license (the "License").  A copy of the License is included in the
"LICENSE.txt" file.  You may also obtain a copy of the License at
http://opensource.org/licenses/BSD-3-Clause .

### Citation

Please cite the following paper if you use ensmallen in your research and/or
software. Citations are useful for the continued development and maintenance of
the library.

* S. Bhardwaj, R. Curtin, M. Edel, Y. Mentekidis, C. Sanderson.
  [ensmallen: a flexible C++ library for efficient function optimization](http://www.ensmallen.org/files/ensmallen_2018.pdf).
  Workshop on Systems for ML and Open Source Software at NIPS 2018.

```
@article{DBLP:journals/corr/abs-1810-09361,
  author    = {Shikhar Bhardwaj and
               Ryan R. Curtin and
               Marcus Edel and
               Yannis Mentekidis and
               Conrad Sanderson},
  title     = {ensmallen: a flexible {C++} library for efficient function optimization},
  journal   = {CoRR},
  volume    = {abs/1810.09361},
  doi       = {10.5281/zenodo.2008650},
  year      = {2018},
  url       = {http://arxiv.org/abs/1810.09361},
  archivePrefix = {arXiv},
  eprint    = {1810.09361},
  timestamp = {Wed, 31 Oct 2018 14:24:29 +0100},
  biburl    = {https://dblp.org/rec/bib/journals/corr/abs-1810-09361},
  bibsource = {dblp computer science bibliography, https://dblp.org}
}
```

### Developers and Contributors

* Ryan Curtin
* Dongryeol Lee
* Marcus Edel
* Sumedh Ghaisas
* Siddharth Agrawal
* Stephen Tu
* Shikhar Bhardwaj
* Vivek Pal
* Sourabh Varshney
* Chenzhe Diao
* Abhinav Moudgil
* Konstantin Sidorov
* Kirill Mishchenko
* Kartik Nighania
* Haritha Nair
* Moksh Jain
* Abhishek Laddha
* Arun Reddy
* Nishant Mehta
* Trironk Kiatkungwanglai
* Vasanth Kalingeri
* Zhihao Lou
* Conrad Sanderson
* Dan Timson
* N Rajiv Vaidyanathan
* Roberto Hueso
* Sayan Goswami
