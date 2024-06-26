# tlbo-omp

## 1º Phase. Implementation

Implementation in C++ of evolution algorithm *TLBO*. 

The pseudocode of the algorithm  is the follow

![pseudocode tlbo](https://github.com/Pmcb04/tlbo-omp/blob/main/pseudocode.png?raw=true)

You can found more information how algorithm works in the paper[*Teaching–learning-based optimization (TLBO): A novel method for constrained mechanical design optimization problems*](https://www.sciencedirect.com/science/article/abs/pii/S0010448510002484?via%3Dihub)


### How to execute

```bash
g++ teaching-learning.cpp -o tlbo 
./tlbo
```


## 2º Phase. Parallelization

In this phase I use *OpenOMP* to obtain the algorithm TLBO can execute in multiple threads. More information of [OpenOMP](https://www.openmp.org/spec-html/5.0/openmp.html).

### How to execute

```bash
g++ teaching-learning-omp.cpp -o tlbo-omp -fopenmp
./tlbo-omp
```