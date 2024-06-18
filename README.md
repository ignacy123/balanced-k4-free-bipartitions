This repository contains files for computer-assisted steps in the paper "Balanced partitions of K_4-free graphs".

# References
We make use of the tooling written by Bernard Lidický, and distributed with the paper "10 Problems for Partitions of Triangle-free Graphs" (authored by József Balogh, Felix Christian Clemen and Bernard Lidický).
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
