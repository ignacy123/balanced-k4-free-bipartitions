This repository contains files for the computer-assisted proof of Theorem 8.1 in the paper "Balanced partitions of K_4-free graphs".

# References
We make use of the tools written by Bernard Lidický, and distributed with the paper "10 Problems for Partitions of Triangle-free Graphs" (authored by József Balogh, Felix Christian Clemen and Bernard Lidický).
You can view the source here: https://lidicky.name/pub/10problems/.

# Dependencies
To calculate the proofs, you need the following dependencies:
 - `g++`
 - `sage`
 - `csdp`.

# Flags format
We encode flags in the following format:
```
[num vertices] [num labeled vertices] [adjacency list of the first vertex]  [adjacency list of second vertex] [adjacency list of third vertex]  ...,
```
where 1 denotes no edge, and 2 denotes an edge.

# Assumptions
Input assumptions from the paper are in the `assump_X.txt` files, and the objective function is in `objective.txt`.
They are written in the form of vectors, (we have performed flag multiplication beforehand).
We multiplied some of them by positive constants greater than 1 to ensure the underlying SDP consists of integers.

# Usage
The `run_all.sh` script performs the following steps:
- encodes Theorem 8.1 as an SDP instance, together with the input assumptions (20), (22), and (24) from the paper.
- obtains an approximate solution to the SDP problem with CSDP.
- rounds the output matrix to a rational, positive semidefinite matrix.
- calculates the final, formal bound obtained from the rounded matrix.
