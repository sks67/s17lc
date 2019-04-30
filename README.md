# s17lc
This is a repository for the radio light curve program used in [Sarbadhicary et al 2017](http://adsabs.harvard.edu/abs/2017MNRAS.464.2326S). **Please cite the paper if you use the code for your project**. Here is the BibteX citation:
```@ARTICLE{2017MNRAS.464.2326S,
   author = {{Sarbadhicary}, S.~K. and {Badenes}, C. and {Chomiuk}, L. and 
	{Caprioli}, D. and {Huizenga}, D.},
    title = "{Supernova remnants in the Local Group - I. A model for the radio luminosity function and visibility times of supernova remnants}",
  journal = {\mnras},
archivePrefix = "arXiv",
   eprint = {1605.04923},
 primaryClass = "astro-ph.HE",
 keywords = {acceleration of particles, ISM: supernova remnants, Local Group, radio continuum: ISM},
     year = 2017,
    month = jan,
   volume = 464,
    pages = {2326-2340},
      doi = {10.1093/mnras/stw2566},
   adsurl = {http://adsabs.harvard.edu/abs/2017MNRAS.464.2326S},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

The code is a simply Python script written with Python 2.7, although it should work with higher versions. 

### Installation
You can do the following :-
1. If you have [git](https://git-scm.com/) installed, you can directly clone this repository using `git clone https://github.com/sks67/s17lc.git` into your computer
2. Copy paste the `s17lc.py` into a blank `.py` file directly in your computer and use it accordingly.

### Description
The program `s17lc.py` is the coded up version of the equations A1-A11 in the paper. 

### Usage
I have included a Jupyter notebook, `test_s17lc.ipynb` to show how to generate radius, velocity, and light curves with the code. 
