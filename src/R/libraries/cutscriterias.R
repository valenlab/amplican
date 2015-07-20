{
# TOO MANY CHANGES, DOCUMENTION NO LONGER VALID!
#
# Defining Aligning Position (AP):
# 
#   An alignment position is defined in the genome/amplicon string with three uppercase letters. For example:
#   
#   aaacccTTTgggggg
#   123456789012345
#   
#   In this case, the alignment position is in position 7. This can vary with respect the alignment. For example:
#     
#     ---aaacccTTTgggggg
#   tttaaacccTTTggg---
#     123456789012345678
#   
#   In here, the aligment position is in position 10.
#   
#   The alignment position is always defined by 3 letters, and could be more than one. But never one with more or less
#   than three letters:
#     
#   nnnNNNnnnNNNnnn: We have TWO alignment positions
#   nnnNNNNNNnnnnnn: We have ZERO alignment positions
#   nnNNnnNNnnNNnnn: We have ZERO alignment positions
#   nnNNNnnnnnNnnNN: We have ONE alignment position.
#   
#   The alignement position has an sphere of influence of certain range. This range is defined by the user and by default
#   is 5 nucleotides distance. In the following example, all 'a' nucleotides are inside the sphere of influence of the AP.
#   
#   ..._____xxx_____...
#   cccaaaaaNNNaaaaaccc
#   8765432100012345678
#   
#   Note that the right side of the AP counts from the last nucleotide of the AP. This means that the 3 letters which
#   compose the AP counts to the right but not to the left, since we defined the AP to be in the position of the first
#   letter.
#   
#   In an alignment, the AP sequence will always be the subject. Hence the cuts are produced when gaps appears in the
#   patterns, and insertions, when appears in the subject.
#   
#   In the alignment, those gaps in the pattern could be within range, or outside range, of the sphere of influence of
#   the AP. For example:
#     
#   aaaacccgg-------ttggagattagggatagagaggaarttaag-------agatagaagatag--------atagata----cgat------------- PRIMER
#   aaaacccgggggggggtTGGagattagggatagagaggaarttaaggccgctaagataGAAgatagggggaaaaaTAGataacaccgatcggttTTTgggtg GENOME
#   
#   _____ A _____                            _____ B _____    _____ C _____      _____ D _____
#   
#   In this example:
#     A has one deletion inside his sphere.
#     B has nothing in his sphere.
#     C has two deletions, one to the left and one to the right
#     D has two deletions, starts from the left and reach the right side overlapping D. Also, one start from the
#       right, even though we can't see the end of it.
#   
# 
# Based on different AP configurations, we have three main cut criterias.
# 
# ALL INSIDE RANGE:
#   
#   This criteria ensure that the cuts/insertions sites are near the alignment positions.
#   
#   Lets define a set of deletions in an alignment called D, and a set of APs called A.
# 
#   We have six posibilities:
# 
#       1.- All members of D are within range of at least one member of A.
#       2.- All members of D are within range of all members of A.
#       3.- At least one member of D is within range of at least one member of A.
#       4.- At least one member of D is within range of all members of A.
#       5.- We don't care if member of D are with in range of members of A.
#       6.- No member of D is within range of any member of A.  
# 
#   This criteria deals with cases 1, 3 and 5.
# 
#   Cases 2, 4 and 6 are over zealous and they are not (yet) implemented.
# 
#   When you analize an alignment, you choose either 1 xor 3 and will tells you which members of D comply with the
#   criteria. In the particular case of 1, either all of them comply, or none of them does. In the case of 5, all of
#   them are always valid deletions.
#
#   Note that this criteria also apply to insertions in the genome.
#
#   This criteria DO NOT APPLY to insertions and deletions at the same time. Menaning that if we have all valid
#   deletions and no valid insertions, the deletions will still be valid, and the insertions won't nullify everything.
# 
# SAME TOTAL CUTS :
#   
#   This criteria deals with the number of valid cuts in a pair of alignments.
# 
#   When we do alignments we do a forward reads with a genome, and a reverse read with the same genome. These alignments
#   do not guaranty to have the same cuts even thought the forward and reverse reads come from the same sequencer.
# 
#   This criteria simply counts how many valid cut are in the forward alignment , and compare them with the number
#   of cuts in the reverse alignment. If both are the same, all cuts are valid. But if they differs, none of the cuts
#   are valid anymore.
#
#   This criteria also applies to the number of valid insertions in the genome for each alignment.
# 
# SAME CUTS INTERVALS :
#   
#   This criteria watch that the cuts in the forwards and reverse alignments starts and end in the same position with
#   respect the genome.
# 
#   First of all, is necessary condition that the number of cuts in the forward and reverse alignment are the same.
#   So this particular criteria is more restrictive than the SAME CUT NUMBER criteria. But you can of course choose
#   to select both, none, or only one of them.
# 
#   Lets see an example and how this criteria apply to it:
#       
#   aaaacccgg-------ttggagattagggatagagaggaarttaag-------agatagaagatag--------atagata---ccgatc------------ FWD
#   -----ccgg-------ttggagattagggatagagaggaarttaag-------agatagaagatag--------atagataa---cgatc-------gggtg RVS
#   aaaacccgggggggggtTGGagattagggatagagaggaarttaaggccgctaagataGAAgatagggggaaaaaTAGataacaccgatcggttTTTgggtg GENOME
#   1234567890          1234567890          1234567890          1234567890          1234567890          12
# 
#               _____ A _____                            _____ B _____    _____ C _____      _____ D _____
# 
#   Note that for this criteria the AP are irrelevant, we are just keeping them for the sake of explaining this
#   criteria when combined with the ALL INSIDE RANGE criteria. Also we are keeping the GENOME free of insertions
#   since this criteria allways apply with the cuts respect genome coordinates. In this case we got this results:
#   
#       A: In this range, both FWD and RVS have a deletion that starts in 10 and ends in 16. This is a good cut.
# 
#       B: In this case, nothing is in range of B so there is no valid cuts. Note that, if you apply this criteria
#          alone, the cut that are in [37,43] would be a valid cut.
# 
#       C: The cut to the left of C is a valid one that starts at 47 and ends at 54. However, the cut to the right
#          is not. Both FWD and RVS have a valid cut because both falls inside the range of C. But the FWD cut is at
#          [52,54] and the RVS cut is at [53,55]. Because both of them do not match, is an invalid cut for this
#          criteria.
# 
#       D: Here, an special case occurs. Is valid cut that is within range of D. However the FWD reads stops before
#          the end appears. This means that we do not have sufficient information to justify this as a valid or invalid
#          cut. For this cases we follow the following approach:
#              - Do we have the 4 interval numbers? Act normal.
#              - Do we have 3 interval numbers?     If only the start or the ends agree, then is a valid cut.
#              - Do we have 2 or less intervals?    Then is an ivalid cut
# 
# FUTURE CRITERIAS:
#     - Cuts must have an specific length.
#     - ...
}

{
# This function deals with the ALL INSIDE RANGE criteria. Please read documentation in the header of this file.
# 
# The function takes the following paramenters:
#   
#   array<int>  alignmentPositions: Represent the positions of the APs respect the alignment
#               eventsStarts: Where the insertions start in the genome
#               eventEnds: Where the insertions ends in the genome
# 
#   int         mode: Represent which mode do you want to apply, can be set from 1 to 5.
#               sphereRange: How many nucleotides away do you want extend the sphere of influence. Minimum 1.
#
# 
# The function returns the following variables:
#   
#   array<bool> , Each position in the array represent a valid event TRUE, or an invalid one FALSE
#   
# Invariant:
#   
#   The legnth of eventStarts and eventEnds and the return array is the same.
#   eventStarts and eventEnds can be Empty, NA or NULL, in this case they have length 0 and NULL will be returned.
#   alignmentPositions can be Empty, NA or NULL, in this case, none of the events are valid.
}

allInsideRange <- function(alignmentPositions, eventsStarts, eventsEnds, mode, sphereRange, positionWide, trailing, alignmentLength){

  # Get how many events we have. In some cases we can have a NA vector or a NULL vector;
  # this is interpreted as a vector of length 0, so no indel in that case.
  if(is.null(eventsStarts[1])){ # You need to check for NULL before NA or else you get a warning. Also you need to
                                # to check only the first element or warnings too.
  
    totalEvents <- 0
  
  } 
  else{
  
    if(is.na(eventsStarts[1])){
      totalEvents <- 0
    }
    else{
      totalEvents <- length(eventsStarts)
    }
    
  
  }
  
  # We can try to call this function with no events, in that case we return NULL
  validEvents <- NULL
  
  if(totalEvents>0){
  
    # Mode 5 means that we don't care what going on in here, so everything that we return is valid, hence TRUE.
    if(mode !=5){
    
      # Prepare the vectors where we are going to track the correct indels
      validEvents <- rep(FALSE,totalEvents)
      
      # For every valid AP position (usually only one)
      if(alignmentPositions>=1){
      for(l in 1:length(alignmentPositions)){

        # For every event that is in the alignment        
        for(m in 1:totalEvents){
                
          # The event can be to the left, to the right, or overlapping the cut point.
                  
          # If the end is smaller than the cut point, then it is to the left
          if(eventsEnds[m]<alignmentPositions[l]){
                    
            # If the end plus the error distance falls to the right or the center of the cut, then is a potential valid
            if(eventsEnds[m] + sphereRange >= alignmentPositions[l]){
              
              # If we don't care if it is at the end of the alignment, then is valid
              if(trailing==TRUE){
                validEvents[m] <- TRUE
              }
              # If we do care about that, lets check out the start
              else{
                
                # If the start is zero, then is a bad cut; otherwise, is good.
                if(eventsStarts[m]!=1){
                  validEvents[m] <- TRUE
                }
              
              }
                      
            }
            
          }
          # Otherwise the events can be to the right or overlapping
          else{
                    
            # If the start is to the right of the cut site, then the event is to the right of the cut site.
            if(eventsStarts[m]>alignmentPositions[l]){
                    
              # If the start plus the error distance falls to the left or the center of the cut, then is valid
              if(eventsStarts[m] - sphereRange <= alignmentPositions[l] + positionWide){
                        
                # If we don't care if it is at the end of the alignment, then is valid
                if(trailing==TRUE){
                  validEvents[m] <- TRUE
                }
                # If we do care about that, lets check out the end
                else{
                  
                  # If the end is the limit, then is a bad cut; otherwise, is good.
                  if(eventsEnds[m]!=alignmentLength){
                    validEvents[m] <- TRUE
                  }
                  
                }
                
                        
              }                
            }
            # Otherwise, is overlapping the cut side, and thus, is valid.
            else{
                      
              # If we don't care if it is at the end of the alignment, then is valid
              if(trailing==TRUE){
                validEvents[m] <- TRUE
              }
              # If we do, we might have a bad start or bad ending
              else{
                
                # Check out both, start and end, to be different from the extremes
                if(eventsStarts[m]!=1 && eventsEnds[m]!=alignmentLength){
                  validEvents[m] <- TRUE  
                }
                
              } 

            }
                    
          }
                
        }
            
      }}
        
      # Now that we are here, we know which events are within range of the AP sites.
      totalValidEvents <- sum(validEvents)
      
      # If mode is set to 1, then we need to have all the events with in range of at least one cut site
      # (We have all, or we have nothing)
      if(mode == 1){
          
          if(totalValidEvents != length(validEvents)){
            
            # If we fail, reset the validEvents because we have nothing.
            validEvents <- rep(FALSE,totalEvents)
           
          }
      }
        
      # If mode is set to 3, then valid events are valid indel, no matter the rest. 
      # We can still have no valid indels thought.
      # In any case, do nothing because mode 2,4 and 6 are not cover (yet).

    
    }
    # If mode 5 is TRUE
    else{
    
      # If we don't care, set everything to true
      validEvents <- rep(TRUE,totalEvents)
      
    }
  
  }
  
  return(validEvents)
  
  
}


{
# This function deals with the SAME TOTAL CUTS criteria. Please read documentation in the header of this file.
#  
# The function simply checks that the total of valid events are the same in both the forward and the reverse alignment.
#
# The function takes the following paramenters:
#  
#     int totalEventForward, totalEventsReverse: How many events are in the forward and the reverse alignment
#
#
# The function returns the following variables:
#
#     bool, TRUE if they are the same, FALSE in another case. Even both of them are 0, the function still return TRUE.
#
#
# Invariant:
#
#     Both totals must be greater or equal to zero.  
}
sameValidsTotal <- function(totalEventsForward, totalEventsReverse){

  sameAmount <- TRUE
  
  # If the number of valid insertions doesn't match in the foward and reverse
  if( totalEventsForward != totalEventsReverse ){
            
    sameAmount <- FALSE
                         
  }
  
  return (sameAmount)

}

{
# This function deals with the SAME CUTS INTERVAL criteria. Please read documentation in the header of this file.
# The function takes the following paramenters:
#
#    array<int> relativeForwardStarts, relativeForwardEnds: Two arrays with the starts and ends of the deletions in
#                                                           the forward alignment which are relative to the genome.
#    
#               relativeReverseStarts, relativeReverseEnds: Two arrays with the starts and ends of the deletions in
#                                                           the reverse alignment which are relative to the genome.
# 
#    array<bool> currentForwardValids:    In here we have an array of booleans values that represent which of those
#                                         events are valid or invalid for the forward alignment.
#
#                currentReverseValids:    In here we have an array of booleans values that represent which of those
#                                         events are valid or invalid for the reverse alignment.
#
#    int genomeLength: This tells us how long is the genome. We need this to know if we are in the special case
#                      were we know the start of the last forward cut but not the end (the end is set up to the
#                      end of the genome, hence the need for this data).
#
# The function returns the following variables:
#
#   A list of two components:
#   [1]  array<bool> , Each position represent a valid event TRUE, or an invalid one FALSE in the forward alignment
#   [2]  array<bool> , Each position represent a valid event TRUE, or an invalid one FALSE in the reverse alignment
#
#   If any of the arrays is Empty, NA, or NULL, we return NULL (see invariant).
#
#
# Invariant:
#     The length of relativeForwardStarts, relativeForwardEnds,  currentForwardValids, and the return list[1]
#     must be the same and greater than 0.
#  
#     The length of relativeReverseStarts, relativeForwardEnds,  currentReverseValids, and the return list[2]
#     must be the same and greater than 0.
#
#     The genomeLength must be greater than 1.  
}
sameCutIntervals <- function(relativeForwardStarts, relativeForwardEnds, relativeReverseStarts, relativeReverseEnds,
                             genomeLength){

  #Variable to test if that the input is valid
  validInput <- TRUE
  
  # Final Return variable, the list will be here.
  toReturn <- list(NULL,NULL)
  
  # Get how many events we have. In some cases we can have a NA vector or a NULL vector;
  # this is interpreted as a vector of length 0, so no indel in that case.
  if(is.null(relativeForwardStarts[1]) || is.null(relativeForwardEnds[1]) ||
     is.null(relativeReverseStarts[1]) || is.null(relativeReverseEnds[1])){ # You need to check for NULL before NA

    validInput <- FALSE
    
  } 
  else{
    
    if(is.na(relativeForwardStarts[1]) || is.na(relativeForwardEnds[1]) ||
       is.na(relativeReverseStarts[1]) || is.na(relativeReverseEnds[1])){
      
      validInput <- FALSE
      
    }
    
  }
  
  if(validInput == TRUE){
  
    # Count how many events we have in the forward and reverse
    totalForwardCuts <- length(relativeForwardStarts)
    totalReverseCuts <- length(relativeReverseStarts)
    
    # Prepare the arrays that you are going to return.
    returningForwards <- rep(FALSE,totalForwardCuts)
    returningReverses <- rep(FALSE,totalReverseCuts) # We set these to all FALSE and mark TRUE when we found a good one.
      
    # For each valid event in the forward we must find a valid event in the reverse with the same coordinates.
    # If we found one, we stop looking and pass to the next one. If we don't, we leave it as invalid (FALSE)
      
    # We need to do the same thing with the reverse. For each valid reverse must be a valid forward. However we can
    # speed things up. We can start assuming that all reverse are invalid. One we find a valid forward, the reverse
    # is also set as valid. In that way, when we finish with the valid forward, the reverse are also done.
      
    # In order to do this we have two pointers, one for the forwards and the other one for the reverse.
    # When we found a valid pair, each pointer set the current array to TRUE, and the forward continues.
  
    # Giving the nature of the problem, a valid forward only have ONE EXCLUSIVE valid reverse. So there is a bijective
    # relationship. Furthermore, if we draw a line between each valid forward/reverse pair, those lines will never
    # cross because they corresponde to the same genome position and we have a bijective application. Because all of
    # this, we have a O(N) and not a O(N²) search here.
      
    # TODO: Find out how the bijective application with total order are call.
      
    forwardIndex <- 1 #?
    reverseIndex <- 1
      
    # Search all indels in the foward
    for(i in 1:totalForwardCuts){
            
      # Start looking from the last position, so this part is O(2N) instead O(N²)
      for(j in reverseIndex:totalReverseCuts){
              
    
        # We have now two posibilities, we have a normal 4 points cut check, or we are in the case where
        # we are at the start/end of one of the reads.
              
        # If this is TRUE, we are on a far away indel.
        if(relativeReverseStarts[j] == 1 || relativeForwardEnds[i] == genomeLength){
        
#           print("3 points case")
          
          # If this condition is TRUE, the boundaries agree and we have a good link.
          # In this case we just need to check that the ends agree or that the starts agree
          if( (relativeReverseStarts == 1          && relativeReverseEnds[j]   == relativeForwardEnds[i]) ||
              (relativeForwardEnds == genomeLength && relativeReverseStarts[j] == relativeForwardStarts[i])  ){
                
            returningForwards[i] <- TRUE
            returningReverses[j] <- TRUE
            reverseIndex <- j
                
          }
        }
        # In the other case, we have a 4 point check
        else{
          
#           print("4 points case")

          
          # If this condition is TRUE, the boundaries agree and we have a good link.
          if(relativeReverseStarts[j] == relativeForwardStarts[i] &&
             relativeReverseEnds[j]   == relativeForwardEnds[i]){
                  
            returningForwards[i] <- TRUE
            returningReverses[j] <- TRUE
            reverseIndex <- j
                  
          }
        }
            
      } # j index
    } # i index

    # Prepare the return results for a valid case. A couple of array full of booleans
    toReturn = list(returningForwards, returningReverses)
  
  }
  
  return (toReturn)

}







sameCutHybrid <- function(relativeForwardStarts, relativeForwardEnds, relativeReverseStarts, relativeReverseEnds,
                             genomeLength){
  
  debug <- FALSE
  intervalsResults <- sameCutIntervals(relativeForwardStarts, relativeForwardEnds,
                                       relativeReverseStarts, relativeReverseEnds,
                                       genomeLength)
  
#   if( sum(intervalsResults[[1]]) > 0){
#   
#     debug <- TRUE
#   
#   }
  
  if(debug == TRUE){
  
    print("Debugging Hybrid")
    print("Forwards")
    print(relativeForwardStarts)
    print(relativeForwardEnds)
    print("Reverses")
    print(relativeReverseStarts)
    print(relativeReverseEnds)
    print("Length")
    print(genomeLength)
  
  }
  
  
  #Variable to test if that the input is valid
  validInput <- TRUE
  
  # Final Return variable, the list will be here.
  toReturn <- list(NULL,NULL)
  
  # Get how many events we have. In some cases we can have a NA vector or a NULL vector;
  # this is interpreted as a vector of length 0, so no indel in that case.
  if(is.null(relativeForwardStarts[1]) || is.null(relativeForwardEnds[1]) ||
       is.null(relativeReverseStarts[1]) || is.null(relativeReverseEnds[1])){ # You need to check for NULL before NA
    
    validInput <- FALSE
    
  } 
  else{
    
    if(is.na(relativeForwardStarts[1]) || is.na(relativeForwardEnds[1]) ||
         is.na(relativeReverseStarts[1]) || is.na(relativeReverseEnds[1])){
      
      validInput <- FALSE
      
    }
    
  }
  
  if(validInput == TRUE){
    
    # Count how many events we have in the forward and reverse
    totalForwardCuts <- length(relativeForwardStarts)
    totalReverseCuts <- length(relativeReverseStarts)
    
    # Prepare the arrays that you are going to return.
    returningForwards <- rep(FALSE,totalForwardCuts)
    returningReverses <- rep(FALSE,totalReverseCuts) # We set these to all FALSE and mark TRUE when we found a good one.
    
    reverseIndex <- 1
    
    # Search all indels in the foward
    for(i in 1:totalForwardCuts){
      
      # Start looking from the last position, so this part is O(2N) instead O(N²)
      if(totalReverseCuts-reverseIndex >= 0){
      for(j in reverseIndex:totalReverseCuts){
        
        # TODO:
        # VERY VERY IMPORTANT
        # The Gotoh has a bug where it returns an interval of [start,end-1] instead of [start,ends] if ends is the
        # end of the amplicon. This piece of code correct that BUT THAT MUST BE CORRECTED IN GOTOH.CPP!!!!
        if((relativeForwardEnds[i]+1) == genomeLength){
          relativeForwardEnds[i] <- relativeForwardEnds[i] + 1
        }
        if((relativeReverseEnds[j]+1) == genomeLength){
          relativeReverseEnds[j] <- relativeReverseEnds[j] + 1
        }
        
        
        # Lets find out the differences between the ends and the starts
        distanceStarts <- relativeForwardStarts[i] - relativeReverseStarts[j]
        distanceEnds   <- relativeForwardEnds[i] - relativeReverseEnds[j]
        
        if(debug == TRUE){
        
          print("Distances")
          print(distanceStarts)
          print(distanceEnds)
          print("Candidates")
          print(i)
          print(j)
          
        }
        
        # First, we look for the start of the I forward indel.
        # If it is equal to 0 we need to do something special about it
        if(relativeForwardStarts[i]==1){
          
          if(debug == TRUE){
          
            print("Case 1")
            
          }
          
          # Any case which is not this, is a valid one
          if( !(distanceEnds < 0 && relativeReverseStarts!=1) ){
          
            returningForwards[i] <- TRUE
            returningReverses[j] <- TRUE
            reverseIndex <- j+1
            
          }
          
        }
        # If neither the start not the end is at the extreme, we have two points to compare
        else if(relativeForwardStarts[i]!=1 && relativeForwardEnds[i]!=genomeLength ){
        
          if(debug == TRUE){
            
            print("Case 2")
            
          }
          
          
          # In this case, the end in reverse match
          # nnnn--------nnnnn
          # nnnn----------???
          if(distanceStarts==0){
            
            # In this case, the end is at the end of the alignment, so doesn't matter what happen with the reverse
            # nnnn-----nnnnn
            # nnnn--------- VALID!
            if(relativeReverseEnds[j]==genomeLength){
            
              returningForwards[i] <- TRUE
              returningReverses[j] <- TRUE
              reverseIndex <- j+1
              
            }
            
            # In this case, the end in both match
            # nnnn------nnn
            # nnnn------nnn VALID!
            else if(distanceEnds==0){  
              
              returningForwards[i] <- TRUE
              returningReverses[j] <- TRUE
              reverseIndex <- j+1
              
            }
          }
          # In this case, the end in reverse match
          # nnnn--------nnnnn
          # nn---------------
          else if(distanceStarts>0 && relativeReverseEnds[j] == genomeLength){
          
            returningForwards[i] <- TRUE
            returningReverses[j] <- TRUE
            reverseIndex <- j+1
            
          }
          # In this case, the end in reverse match
          # nnnn--------nnnnn
          # ------------nnnnn
          
          # nnnn--------nnnnn
          # ---------------nn
          else if(distanceEnds<=0 && relativeReverseStarts[j]==1){
          
            returningForwards[i] <- TRUE
            returningReverses[j] <- TRUE
            reverseIndex <- j+1
          
          }
        }
        # If the end is the actuall end of the alignment
        else if(relativeForwardEnds[i]==genomeLength){
        
          if(debug == TRUE){
            
            print("Case 3")
            
          }
          
          
          # Any other case which IS NOT this one, is good.
          # nnnnnnn----------
          # nnnn---------nnn
          if(  !(distanceStarts>0 && relativeReverseEnds[j] != genomeLength && relativeReverseStarts[j]!=1) ){
          
            returningForwards[i] <- TRUE
            returningReverses[j] <- TRUE
            reverseIndex <- j+1
            
          }
        }
        
        if(debug == TRUE){
          
          print("Indexes")
          print(reverseIndex)
          print(i)
          print(j)
          
        }
        
      }} # j index + if
    } # i index
    
    # Prepare the return results for a valid case. A couple of array full of booleans
    toReturn = list(returningForwards, returningReverses)
    
  }
  
  return (toReturn)
  
}







# This function test the primer dimer criteria.
# For a given list of events, we need to check each one and if it meet the condition that the size of the event
# is smaller than a certain range.
# The range is given by the size of the amplicon minus the sumation of the lengths of the primers plus a constant.
# The constant can be choose by the user, the default is 10.
primerDimerTest <- function(eventsStarts, eventsEnds, forwardPrimerLength, reversePrimerLength, ampliconLength, base = 10){

  # Get how many events we have. In some cases we can have a NA vector or a NULL vector;
  # this is interpreted as a vector of length 0, so no indel in that case.
  if(is.null(eventsStarts[1])){ # You need to check for NULL before NA or else you get a warning. Also you need to
    # to check only the first element or warnings too.
    
    totalEvents <- 0
    
  } 
  else{
    
    if(is.na(eventsStarts[1])){
      totalEvents <- 0
    }
    else{
      totalEvents <- length(eventsStarts)
    }
    
    
  }
  
  # We can try to call this function with no events, in that case we return NULL
  validEvents <- NULL

  if(totalEvents>0){
  
    # Prepare the vectors where we are going to track the correct indels
    validEvents <- rep(FALSE,totalEvents)
    
    limitLength <- ampliconLength - (reversePrimerLength + forwardPrimerLength + base)
    
    # For every event that is in the alignment        
    for(i in 1:totalEvents){
      
      cutLength <- eventsEnds[i] - eventsStarts[i] + 1
      
      if(cutLength < limitLength){
      
        validEvents[i] <- TRUE
        
      }
      
    }
    
  }
  
  return (validEvents)

}

# Get the deletions in the pattern and tells if the sequence is frameshifted
frameShift <- function(eventsStarts, eventsEnds, alignmentLength){
  
  # We have the following return scenarios
  # 0 -- Something when wrong, we must never get this
  # 1 -- No deletions
  # 2 -- One deletion in frame shifted
  # 3 -- One deletion out frame shifted
  # 4 -- Multiple deletions
  frameClass <- 0
  
  totalEvents <- 0
  totalDeletions <- 0
  allDeletionsLength <- 0
  
  # Get how many events we have. In some cases we can have a NA vector or a NULL vector;
  # this is interpreted as a vector of length 0, so no indel in that case.
  if(is.null(eventsStarts[1])){ # You need to check for NULL before NA or else you get a warning. Also you need to
    # to check only the first element or warnings too.
    totalEvents <- 0
  } 
  else{
    if(is.na(eventsStarts[1])){
      totalEvents <- 0
    }
    else{
      totalEvents <- length(eventsStarts)
    }
  }
  
  # Fill the deletions stats
  if(totalEvents>0){
   
    for(m in 1:length(eventsStarts)){
    
      if(eventsStarts[m] != 1 && eventsEnds[m] != alignmentLength){
      
        totalDeletions <- totalDeletions + 1
        allDeletionsLength <- allDeletionsLength + (eventsEnds[m] - eventsStarts[m])
        
      }
    }
  }
  
  # Find out the class where we are now
  # -- If we have no deletion, we are class 1
  if(totalDeletions == 0){
    frameClass <- 1
  }
  else{
  
    # -- If we have more than 1, we are in class 4
    if(totalDeletions > 1){
      frameClass <- 4
    }
    else{
      # -- If the deletion is modulus 3 we are in frame shift
      if(allDeletionsLength %% 3 == 0){
      
        frameClass <- 2
      }
      # -- If not, we are out of frame
      else{
      
        frameClass <- 3
      }
      
    }
    
  }

  return (frameClass)

}