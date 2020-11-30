# AMUSIC
AMUSIC is a tool for approximate counting of minimal unsatisfiable subsets of a given Boolean formula in CNF. 

## Installation
AMUSIC employs several third-party tools that you have to install:

#### PySAT
AMUSIC uses the python API of PysAT [3]. You can install PySAT via `pip3 install python-sat`. In case that does not work, follow instructions at https://github.com/pysathq/pysat 

#### CAQE
Download CAQE [4] from https://github.com/ltentrup/caqe and build it. Copy the binary of caqe (located at caqe/target/release/caqe) to the folder amusic/tools/.

#### QRATPre+
Download QRATPre+ [5] from https://github.com/lonsing/qratpreplus and build it.  Copy the binary of qratpre+ (located at qratpreplus/qratpre+) to the folder amusic/tools/. 

#### CADET
Download CADET [6] from https://github.com/MarkusRabe/cadet and build it. Copy the binary of cadet to the folder amusic/tools/.

Note that we use also a single MUS extractor muser2 [2]; we distribute its binary with our tool (you do not have to build it).
After installing all the necessary tools, the folder amusic/tools/ should contain four binaries: caqe, muser2-para, qratpre+, cadet.

## Running our tool
You can run the tool via "python3 counter.py <input_file>", e.g.:
```
python3 counter.py examples/generated_16.gcnf
```
To run the tool with a time limit use e.g.:
```
timeout 300 python3 counter.py examples/generated_16.gcnf
```
To use a 2QBF solver as a backend instead of a 3QBF solver, use the flag "--qbf2", e.g.:
```
python3 counter.py --qbf2 examples/generated_16.gcnf
```
To see all the available parameters, run:
```
python3 counter.py -h
```

## Related Tools
AMUSIC is the first (and as far as we know the only one) probabilistic approximate counter of Minimal Unsatisfiable Subsets. If you want to find the exact MUS count, you should use a complete MUS enumeration algorithm to enumerate all the MUSes. In particular, for the exact MUS count, give a try to our MUS enumeration tool UNIMUS [1]: https://github.com/jar-ben/unimus

Also, note that AMUSIC is suitable especially in the case where the MUS count is too large (e.g., millions of MUSes) to be able to be determined by a complete MUS enumeration algorithm in a reasonable time. However, if the MUS count is relatively low, we recommend to use a complete MUS enumeration algorithm (https://github.com/jar-ben/unimus) as it will be probably much faster than AMUSIC. 

## Citation
If you use AMUSIC in your research, please, cite our paper that presented AMUSIC:
```
@inproceedings{DBLP:conf/cav/BendikM20,
  author    = {Jaroslav Bend{\'{\i}}k and
               Kuldeep S. Meel},
  title     = {Approximate Counting of Minimal Unsatisfiable Subsets},
  booktitle = {{CAV} {(1)}},
  series    = {Lecture Notes in Computer Science},
  volume    = {12224},
  pages     = {439--462},
  publisher = {Springer},
  year      = {2020}
}
```


## References

* [1] https://github.com/jar-ben/unimus
* [2] https://bitbucket.org/anton_belov/muser2
* [3] https://pysathq.github.io/
* [4] https://github.com/ltentrup/caqe
* [5] https://github.com/lonsing/qratpreplus
* [6] https://github.com/MarkusRabe/cadet

## Contact
In case of any troubles, do not hesitate to contact me, Jaroslav Bendik, the developer of the tool, at xbendik=at=gmail.com.
