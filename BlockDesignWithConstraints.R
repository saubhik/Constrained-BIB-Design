# This code runs a Partially Balanced Incomplete Block design,
# where `constraints` such as some pairs of treatments cannot appear together in any block are taken care of.
# 
# `Incomplete` means that number of treatments per block is less than the total number of treatments
# `Balanced` means that all pairwise treatments occur equal number of times in the block design
# Also, a treatment cannot occur more than once for a single block and the frequencies of each treatment
# should be equal.
#
# `optBlock` from `AlgDesign` gives us a PBIB design.
# This code modifies the design output from `optBlock` so as to accommodate for the constraints in a manner
# that the design deviates the least from what we mean from a BIB design
# 
# The algorithm relies on a `swapping` technique to swap a treatment of the `forbidden pairs` with another
# treatment from another block in a way which minimizes the variance of the lower triangular co-occurrence matrix.

rm(list=ls())
library(AlgDesign)

number.alternatives = 123
number.blocks = 123
alternatives.per.block = 20
constraint.pairs = list(c(3, 4), c(11, 12), c(19, 20), c(27, 28), c(35, 36), c(43, 44), c(51, 52),
                        c(59, 60), c(68, 69), c(76, 77), c(84, 85), c(92, 93), c(100, 101), c(108, 109),
                        c(117, 118))
n.repeats=5


# run the BIB design
set.seed(123)
results <- optBlock(~.,
                    withinData=factor(1:number.alternatives),
                    blocksizes=rep(alternatives.per.block, number.blocks), 
                    nRepeats=n.repeats)


# get the co-occurrence matrix for a design
get.cooccurence.matrix = function(design) {

  co.occurrence.matrix = crossprod(table(c(rep(1:number.blocks, 
                                          rep(alternatives.per.block, number.blocks))), 
                                          design))
  return (co.occurrence.matrix)
}


# find the set of blocks not containing any of the `forbidden` pairs
# this will be used later in the algorithm, in the next code block
get.legit.blocks = function(results, constraint.pairs) {
  block.set.rep = c()
  block.id = 0
  for (block in results$Block) {
    block.id = block.id + 1
    items.in.current.block = as.numeric(as.character(unlist(block)))
    
    # find if the block is legit or not
    legit.block = TRUE
    for (pair in constraint.pairs) {
      if (all(pair %in% items.in.current.block)) {
        legit.block = FALSE
        break
      }
    }
    
    # save the block id if it's a legit block
    if (legit.block) {
      block.set.rep = c(block.set.rep, block.id)
    }
  }
  
  return (block.set.rep)
}


# run the function above to get the legit blocks that can be used
# for replacement
block.set.rep = get.legit.blocks(results = results, constraint.pairs = constraint.pairs)


# get the design from the results
design = as.numeric(results$design[,1])


# get table from design
get.design.table = function(design) {
  return (matrix(design, nrow = number.blocks, ncol = alternatives.per.block, byrow = TRUE))
}


# core function for swapping based on minimum variance of lower
# diagonal entries of co-occurrence matrix
swapper = function(design, pair, block.id, block.set.rep) {
  
  design.table = get.design.table(design)
  
  # function to calculate the variance of the lower diag entries
  # of co-occurrence matrix
  get.variance = function(design) {
    coocmat = get.cooccurence.matrix(design)
    return (var(coocmat[lower.tri(coocmat)]))
  }
  
  # get the numbers which cannot be used as candidates
  forbidden = design.table[block.id, ]
  for (p in constraint.pairs) {
    if (any(p %in% design.table[block.id, ])) {
      forbidden = c(forbidden, p[1], p[2])
    }
  }
  
  # SOURCE BLOCK CONDITION CHECK
  # get the candidate alternatives which can be used for swapping
  # make sure that these candidate alternatives doesn't contain:
  # 1. any number from original block.id
  # 2. any number which when combined with a number from the original
  # block.id would result in a forbidden pair
  # we have already computed these 2 types of numbers in `forbidden`
  design.rep.set = design.table[block.set.rep,]
  candidates = split(design.rep.set, block.set.rep)
  candidates = lapply(candidates, function(c) base::setdiff(c, forbidden))
  
  min.variance = Inf
  best.design.table = NULL
  for (elem1 in pair) {
    
    # uncomment next line for debugging
    # print(elem1)
    
    for (ind in 1:length(candidates)) {
      
      block.rep.id = as.numeric(names(candidates[ind]))
      
      # DESTINATION BLOCK CONDITION CHECK
      # if swap element is already present in destination block or
      # the other element of the pair of swap element is present in
      # destination block then don't consider that destination block
      if (any(pair %in% design.table[block.rep.id, ])) {
        next
      }
      
      for (elem2 in candidates[[ind]]) {
        
        # swap elem1 and elem2
        design.table.copy = design.table
        design.table.copy[block.id, match(elem1, design.table.copy[block.id,])] = elem2
        design.table.copy[block.rep.id, match(elem2, design.table.copy[block.rep.id,])] = elem1
        
        # check the variance
        if (get.variance(design.table.copy) < min.variance) {
          
          # uncomment following lines for debugging
          # print(paste0("swapped ", block.id, " and ", block.rep.id))
          # print(paste0("swapped elements ", elem1, " and ", elem2))
          
          min.variance = get.variance(design.table.copy)
          best.design.table = design.table.copy
        }
      }
    }
  }
  
  return (matrix(t(best.design.table), nrow = number.blocks * alternatives.per.block, ncol = 1)[,1])
}


# get the blocks containing a `forbidden` pair
illegit.blocks = base::setdiff(seq(1, number.blocks), block.set.rep)


st = Sys.time()
# algorithm for each constraint pair
for (pair in constraint.pairs) {
  # for each block containing a `forbidden` pair
  for (block.id in illegit.blocks) {
    items.in.current.block = as.numeric(as.character(unlist(results$Blocks[[block.id]])))
    if (all(pair %in% items.in.current.block)) {
      # update the design after solving each pair
      design = swapper(design, pair, block.id, block.set.rep)
    }
  }
}
en = Sys.time()
print(en - st)


check.constraints = function(design) {
  coocmat = get.cooccurence.matrix(design)
  for (pair in constraint.pairs) {
    if (coocmat[pair[1], pair[2]] != 0) {
      return (FALSE)
    }
  }
  return (TRUE)
}

# check whether all the constraints are satisfied: no forbidden pair occurs together
check.constraints(design)

write.csv(get.design.table(design), file = "output_design.csv")
