# Su-Schrieffer-Heeger-Hubbard Model
The repository for the research on Su-Schrieffer-Heeger-Hubbard model, including codes, notes, plot subscripts, and manuscript.
The Hamiltonian of the model reads
$$
\begin{aligned} 
{H}=&-t \sum_{\langle i, j\rangle, \sigma}\left({c}_{i, \sigma}^{\dagger} {c}_{j, \sigma}+\text { H.c. }\right) \\ &+\frac{U}{2} \sum_{i} {n}_{i}^{2}+\alpha \sum_{i} {n}_{i} \hat{X}_{i}+\sum_{i}\left[\frac{\hat{P}_{i}^{2}}{2 m}+\frac{k \hat{X}_{i}^{2}}{2}\right]. 
\end{aligned}
$$

## Reference
The research on the SSHH model was summarized in [arXiv:2211.09143](https://arxiv.org/abs/2211.09143).

## Author
Hao-Xin Wang  <wanghx18@mails.tsinghua.edu.cn>

## Dependence
- [QuantumLiquids/tensor](https://github.com/QuantumLiquids/tensor), branch dev-v0.2-MagicChangeExpansion
- [QuantumLiquids/MPS2](https://github.com/QuantumLiquids/MPS2), branch dev-v0.2-OptimizeLanczos
- Boost::mpi, Boost::serialization, Boost headers
- MPI
- Intel MKL
- [springer13/hptt](https://github.com/springer13/hptt)


## Usage
First install above dependencies. Note the Boost::mpi should keep the same version with bard MPI.
Then configure the project in a directory `build` by CMakeLists.txt. 
Remember to hint the dependencies path to Cmake. After automatically finding all the dependencies and 
configuration, run
```bash
make -j16
```
to compile the binary files.

To get the ground state MPS of the SSHH model, firstly run
```bash
./mpogen params.json
```
and  it will generate the MPO files of the SSHH model. Then run
```bash
mpirun -n $NUMPROC ./vmps_ssh_pbc params.json
```
to utilize variational MPS algorithm (people always call it DMRG) to optimize ground state.
You should run it once and once again and adjust params.json at the same time.

Other programs in CMake targets are used to measure correlations, or calculate for open boundary conditions, or calculate for the Holstein-Hubbard model.

## Calculation time scale
For a $4\times 40$ lattice with phonon pseudo-site number 3, if you're lucky so that the cluster is stable,
and if you're working hard enough so that the program is always running under appropriate parameters,
using about 240 CPU cores at the same time, it requests about 1 year to push the bond dimension to $D\sim 18000$, 
due to the complexity induced by the phonon hilbert space.