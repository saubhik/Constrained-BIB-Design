# Constrained Partially Balanced Incomplete Block (PBIB) Design

Suppose that we have `N` treatments (or alternatives) and `B` blocks, with `k` treatments per block (or alternatives per block). We have:

1. `k < N`, or in other words an *incomplete* design
2. A treatment can occur only once in a block `b`
3. Each treatment occur equal number of times (equal frequencies)
4. Each treatment pair occur equal number of times (we try to maintain this), or in other words a *partially balanced* design

We call a matrix of dimension `(N, N)`, whose `(i, j)`th entry shows the pairwise frequency of `i` & `j` occurring together, a `co-occurrency matrix`. Call this `M`.

`M` is a symmetric matrix. Ideally, we would like to have the variance of the lower triangular entries to be `0`, i.e. all the lower triangular entries should be same (equal pairwise frequencies) and equal diagonal entries (equal treatment frequencies).

In addition to the above 4, we can have some additional constraints which can arise in any application:

5. `p` treatment pairs `[(i1, j1), (i2, j2), ..., (ip, jp)]` cannot occur together. Call these `forbidden pairs`.

I use `optBlock` from the `R` package `AlgDesign` for getting a PBIB design.

And then I use a swapping-based algorithm to satisfy the constraint 5, which I explain below:

### Algorithm
* **Step 0** : For every pair `(i, j)` in the `p` treatment pairs, repeat step 1.
* **Step 1** : For every block `b` containing `(i, j)`, repeat steps 2 through 4.
* **Step 2** : We will be swapping either `i` or `j` with some other treatment. Find the set of all blocks that doesn't contain `(i, j)`. Call this set `Brep`.
* **Step 3** : Choose a treatment from a block `brep` in `Brep` which doesn't contain any of the treatments in `b`, and also when combined with a treatment from `b` doesn't become a forbidden pair.
* **Step 4** : Choose the treatment in Step 3, and `i` or `j` which minimizes the variance of the lower triangular entries of the co-occurrence matrix. Make sure that when we are considering `i`, `j` should not be in the block `brep` and vice versa.

So we see that a swap of `x` from block `bx` with `y` from block `by` is `valid` if:
1. `y` is not in `bx`.
2. `x` is not in `by`.
3. `y` does not form a `forbidden pair` with any element in `bx`.
4. `x` does not form a `forbidden pair` with any element in `by`.

This ensures that the occurrences of `forbidden pairs` decreases with every iteration of the algorithm, in a way that keeps the design as close to BIB as possible.

How to make the output design invariant of the order of the constraints?