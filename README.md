# AMUSIC
AMUSIC is a tool for approximate counting of minimal unsatisfiable subsets of a given Boolean formula in CNF. 

## Installation
AMUSIC employs several third-party tools that you have to install:

#### PySAT
AMUSIC uses the python API of PysAT [3]. You can install PySAT via `pip3 install python-sat`. In case that does not work, follow instructions at https://github.com/pysathq/pysat 

#### CAQE
Download CAQE [4] from https://github.com/ltentrup/caqe and build it. Copy the binary of caqe (located at caqe/target/release/caqe) to the folder amusic/tools/.

#### QRATPre+
Download QRATPre+ [5] from https://github.com/ltentrup/caqe and build it.  Copy the binary of qratpre+ (located at qratpreplus/qratpre+) to the folder amusic/tools/. 

Note that we use also a single MUS extractor muser2 [2]; we distribute its binary with our tool (you do not have to build it).
After installing all the necessary tools, the folder amusic/tools/ should contain three binaries: caqe, muser2-para, qratpre+.

## Running our tool
You can run the tool via "python3 counter.py <input_file>", e.g.:
```
python3 counter.py examples/generated_16.gcnf
```
To run the tool with a time limit use e.g.:
```
timeout 300 python3 counter.py examples/generated_16.gcnf
```
To see all the available parameters, run:
```
python3 counter.py -h
```

## Related Tools
AMUSIC is the first (and as far as we know the only one) probabilistic approximate counter of Minimal Unsatisfiable Subsets. If you want to find the exact MUS count, you should use a comple MUS enumeration algorithm to enumerate all the MUSes. In particular, for the exact MUS count, give a try to our MUS enumeration tool: https://github.com/jar-ben/mustool

Also, note that AMUSIC is suitable especially in the case where the MUS count is too large (e.g., millions of MUSes) to be able to be determined by a complete MUS enumeration algorithm in a reasonable time. However, if the MUS count is relatively low, we propose to use a complete MUS enumeration algorithm (https://github.com/jar-ben/mustool) as it will be probably much faster than AMUSIC. 

## References

* [1] Jaroslav Bendík and Kudeep S. Meel: Approximate Counting of Minimal Unsatisfiable Subsets. CAV 2020.
* [2] https://bitbucket.org/anton_belov/muser2
* [3] https://pysathq.github.io/
* [4] https://github.com/ltentrup/caqe
* [5] https://github.com/lonsing/qratpreplus


## Contact
In case of any troubles, do not hesitate to contact me, Jaroslav Bendik, the developer of the tool, at xbendik=at=gmail.com.
