# amusic
AMUSIC: a tool for approximate counting of minimal unsatisfiable subsets of a given Boolean formula in CNF

## Installation
TBA

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

## Contact
In case of any troubles, do not hesitate to contact me, Jaroslav Bendik, the developer of the tool, at xbendik=at=gmail.com.
