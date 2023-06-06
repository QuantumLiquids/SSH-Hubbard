# SSH-Hubbard
The repository for the research on Su-Schrieffer-Heeger-Hubbard model, including codes, notes, plot subscripts, and manuscript.

## Reference
The research on the SSHH model was summarized in [arXiv:2211.09143](https://arxiv.org/abs/2211.09143).

## Author
Hao-Xin Wang  <wanghx18@mails.tsinghua.edu.cn>

## Dependence
- LinuxDaFaHao/tensor, branch dev-v0.2-MagicChangeExpansion
- LinuxDaFaHao/MPS2, branch dev-v0.2-OptimizeLanczos
- Boost::mpi, Boost::serialization, Boost headers
- MPI
- Intel MKL
- HPTT


## Usage
Compile by CMakeLists.txt. RTFM of CMake.

To get the ground state MPS, firstly run
```bash
./mpogen params.json
```
to get the MPO files of the model. Then run
```bash
mpirun -n $NUMPROC ./vmps_ssh_pbc params.json
```
to utilize variational MPS algorithm (generally it is called DMRG) to optimize ground state.
You should run it once and once again and adjust params.json at the same time.

Other programs in CMake targets are used to measure correlations, or calculate for open boundary conditions, or calculate for the Holstein-Hubbard model.

## Calculation time scale
For a $4\times 40$ lattice with phonon pseudo-site number 3, if you're lucky so that the cluster is stable,
and if you're working hard enough so that the program is always running under appropriate parameters,
using about 240 CPU cores at the same time, it requests about 1 year to push the bond dimension to $D\sim 18000$, 
due to the complexity induced by the phonon hilbert space.