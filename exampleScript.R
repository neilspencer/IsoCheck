# example script demonstrating the isomorphism check for spreads
# and stars using the IsoCheck package
library(IsoCheck)

############
# Notation #
############

# Spreads and stars are formatted as 3-dimensional arrays. For a spread or star
# spr, the entry spr[i,j,k] indicates the presence of the ith basic factor
# in the jth effect of the kth flat of spr. The following provides examples
# in which we build the representation of a 1-spread of PG(3,2) directly
# and using shorthand functions.

# direct instantiation of the spread
spread <- spread <- array(NA, c(4,3,5))

spread[,1,1] <- c(0, 0, 0, 1)
spread[,2,1] <- c(0, 1, 1, 0)
spread[,3,1] <- c(0, 1, 1, 1)
spread[,1,2] <- c(0, 0, 1, 0)
spread[,2,2] <- c(1, 1, 0, 0)
spread[,3,2] <- c(1, 1, 1, 0)
spread[,1,3] <- c(0, 1, 0, 0)
spread[,2,3] <- c(1, 0, 1, 1)
spread[,3,3] <- c(1, 1, 1, 1)
spread[,1,4] <- c(1, 0, 0, 0)
spread[,2,4] <- c(0, 1, 0, 1)
spread[,3,4] <- c(1, 1, 0, 1)
spread[,1,5] <- c(0, 0, 1, 1)
spread[,2,5] <- c(1, 0, 1, 0)
spread[,3,5] <- c(1, 0, 0, 1)


# using shorthand
stringtovector <- function(string){
  ret <- rep(0, 4)
  ret[1] <- (grepl("A", string)==T)
  ret[2] <- (grepl("B", string)==T)
  ret[3] <- (grepl("C", string)==T)
  ret[4] <- (grepl("D", string)==T)
  return(ret)
}

spread <- array(NA, c(4,3,5)) #

spread[,1,1] <- stringtovector("D")
spread[,2,1] <- stringtovector("BC")
spread[,3,1] <- stringtovector("BCD")
spread[,1,2] <- stringtovector("C")
spread[,2,2] <- stringtovector("AB")
spread[,3,2] <- stringtovector("ABC")
spread[,1,3] <- stringtovector("B")
spread[,2,3] <- stringtovector("ACD")
spread[,3,3] <- stringtovector("ABCD")
spread[,1,4] <- stringtovector("A")
spread[,2,4] <- stringtovector("BD")
spread[,3,4] <- stringtovector("ABD")
spread[,1,5] <- stringtovector("CD")
spread[,2,5] <- stringtovector("AC")
spread[,3,5] <- stringtovector("AD")

spread

############
#Example 1 #
############

# load two example 1-spreads of PG(3,2) from the Isocheck package
data(spreadn4t2a)
data(spreadn4t2b)

# test their isomorphism
test1 <- checkSpreadIsomorphism(spreadn4t2a, spreadn4t2b)

test1$result # the test indicates that they are isomorphic
(IEC1 <- (test1$IECs)[[1]])
# we store the first isomorphism establishing collineation as IEC1

# Let us verify that the function is working properly

# apply IEC1 to spread4t2a
(Cspreadn4t2a <- applyCollineation(IEC1, spreadn4t2a))
spreadn4t2b # By visual inspection it appears to be the same as spreadn4t2b

# bitstring representation verifies this fact
prod(getBitstrings(Cspreadn4t2a) == getBitstrings(spreadn4t2b))

# we could also use the checkSpreadEquivalence function
checkSpreadEquivalence(Cspreadn4t2a, spreadn4t2b)


############
#Example 2 #
############

# For a bigger problem, consider two 2-spreads of PG(5,2)
data(spreadn6t3a)
data(spreadn6t3b)

# we are using the option to return just the first IEC to cut down on runtime
test2 <- checkSpreadIsomorphism(spreadn6t3a, spreadn6t3b, returnfirstIEC = T)

test2$result # the test indicates that they are isomorphic
(IEC2 <- (test2$IECs)[[1]])
# we can again verify the isomorphism
Cspreadn6t3a <- applyCollineation(IEC2, spreadn6t3a) # apply IEC2 to spread6t3a
checkSpreadEquivalence(Cspreadn6t3a, spreadn6t3b) # verified!


############
#Example 3 #
############

# For a problem involving two non-isomorphic spreads,
#consider two 1-spreads of PG(5,2)
data(spreadn6t2a)
data(spreadn6t2c)
test3 <- checkSpreadIsomorphism(spreadn6t2a, spreadn6t2c, returnfirstIEC = T)
# The above test takes a bit of time

test3$result # the test indicates that they are not isomorphic

############
#Example 4 #
############

# Consider the two 1-spreads of PG(5,2)
data(spreadn6t2a)
data(spreadn6t2b)
test4 <- checkSpreadIsomorphism(spreadn6t2a, spreadn6t2b, returnfirstIEC = T)
test4$result # the test indicates that they are isomorphic

############
#Example 5 #
############

# Now, we consider stars.
# First, we consider two stars of PG(4,2) consisting
# of 4-flats.
data(starn5t3a)
data(starn5t3b)
test5 <- checkStarIsomorphism(starn5t3a, starn5t3b, returnfirstIEC = T)
test5$result # the test indicates that they are isomorphic
(IECstar <- test5$IECs[[1]]) # the first IEC
checkStarEquivalence(applyCollineation(IECstar, starn5t3a), starn5t3b) # verify


############
#Example 6 #
############

# As another example , we consider two stars of PG(7,2) consisting
# of 6-flats.
data(starn8t5a)
data(starn8t5b)
test6 <- checkStarIsomorphism(starn8t5a, starn8t5b, returnfirstIEC = T)
test6$result # the test indicates that they are also isomorphic
