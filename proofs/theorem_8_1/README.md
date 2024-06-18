In this directory, you can find the proof of Theorem 8.1.

The `run_all.sh` script performs the following steps:
- encodes Theorem 8.1 as an SDP instance, together with assumptions (20), (22), and (24) from the paper.
- obtains an approximate solution to the SDP problem with CSDP.
- rounds the output matrix to a rational, positive semidefinite matrix.
- calculates the final, formal bound obtained from the rounded matrix.

Input assumptions from the paper are in the `assump_X.txt` files, and the objective function is in `objective.txt`.
They are written in the form of vectors, (we have performed flag multiplication beforehand).
We multiplied some of them by positive constants larger than 1 to ensure the underlying SDP instance is integer.
