hartree-fock
============

A simple one dimensional restricted closed-shell Hartree-Fock implementation 
based upon the examples in "Modern Quantum Chemistry" by Szabo and Ostlund.

The program is written in Julia and when run it will calculate the energy
for H2 and HeH+ and verify that the energies are correct. It will also
calculate the energy of HeH+ for a range of bond lengths and write this
information to the file `heh_pes.dat`. The potential energy surface
can be visualized with the supplied gnuplot script `heh_pes.gp`.
