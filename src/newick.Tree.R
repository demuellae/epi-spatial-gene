"newick.Tree" <-
    function(P, numLeaf, N, Z) {
        newick <- array(0, dim=(199))
        newick[1:numLeaf] <- Z[1:numLeaf]
        for (parent in seq(numLeaf+1,N)) {
            children <- which(P==parent)
            newick[parent] <- paste( "(", toString(newick[children[1]]), ",",  toString(newick[children[2]]) ,  ")",  toString(Z[parent]), sep="")
        }
        paste(newick[N], ";")
    }
