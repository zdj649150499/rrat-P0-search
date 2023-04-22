The software is to search the rotation period of RRAT using a single pulse with the method of Keane et al.(2010) (doi: 10.1111/j.1365-2966.2009.15693.x)

This software is used in the FAST GPPS survey: http://zmtt.bao.ac.cn/GPPS/.



## Install

```shell
gcc rratbestP0.c  -o rratbestP0  -lgsl -lgslcblas -lcpgplot -lpgplot -lX11 -lgfortran -lpng -lz -O3 -lm -Wall
```

'-lpng -lz'  could not included if you PGPLOT not installed with png.





## Usage

```shell
rratbestP0  -bp 0.1,10,0.00001 TOA.list
```

1.  The value '0.1,10,0.00001' is the min, max, and step of period for searching period.
2.  The TOA.list file is a file with only TOAs (in MJD) in one column.



A 'TOA-P0-sigma.ps' image file will out put, and the P0 and sigma will output in the terminal.



If you used this software, plese cite:

>   @ARTICLE{2023arXiv230317279Z,
>          author = {{Zhou}, D.~J. and {Han}, J.~L. and {Xu}, Jun and {Wang}, Chen and {Wang}, P.~F. and {Wang}, Tao and {Jing}, Wei-Cong and {Chen}, Xue and {Yan}, Yi and {Su}, Wei-Qi. and {Gan}, Heng-Qian and {Jiang}, Peng and {Sun}, Jing-Hai and {Wang}, Hong-Guang and {Wang}, Na and {Wang}, Shuang-Qiang and {Xu}, Ren-Xin and {You}, Xiao-Peng},
>           title = "{The FAST Galactic Plane Pulsar Snapshot Survey: II. Discovery of 76 Galactic rotating radio transients and their enigma}",
>         journal = {arXiv e-prints},
>        keywords = {Astrophysics - High Energy Astrophysical Phenomena},
>            year = 2023,
>           month = mar,
>             eid = {arXiv:2303.17279},
>           pages = {arXiv:2303.17279},
>             doi = {10.48550/arXiv.2303.17279},
>   archivePrefix = {arXiv},
>          eprint = {2303.17279},
>    primaryClass = {astro-ph.HE},
>          adsurl = {https://ui.adsabs.harvard.edu/abs/2023arXiv230317279Z},
>         adsnote = {Provided by the SAO/NASA Astrophysics Data System}
>   }
