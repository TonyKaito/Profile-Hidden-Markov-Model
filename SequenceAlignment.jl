function DiagOneMatrix(n)
  #return diag of 1 matrix
  matrix = zeros((n, n))
  for j in range(1, length=n)
    for i in range(1, length=n)
      if i == j
        matrix[i, j] = 1
      end
    end
  end
  return matrix
end

function FindPairAlignment(s, w, indel=1, find="min")
  s_new = ""
  w_new = ""
  DistanceMatrix = [[]]
  
  #get the matrix
  for j in range(1, length=length(w)+1)
    for i in range(1, length=length(s)+1)
      if i == 1
        push!(DistanceMatrix[i], j-1)
      elseif j == 1
        push!(DistanceMatrix,[i-1])
      else
        above = DistanceMatrix[i-1][j]
        above_delta = indel
        left = DistanceMatrix[i][j-1]
        left_delta = indel
        # to be done, get scoring matrix, but default is of edit distance
        diag = DistanceMatrix[i-1][j-1]
        diag_delta = 1
        
        if s[i-1] == w[j-1]
          diag_delta = 0
        end
            
        value = min(above + above_delta, 
                    left + left_delta, 
                    diag + diag_delta)
        push!(DistanceMatrix[i], value)
      end
    end
  end
  dist = last(last(DistanceMatrix))
  
  # get the alignment
  row = length(s)+1
  col = length(w)+1
  while col > 1 || row > 1
    # border case of searching for values outside the matrix
    if row-1 > 0
      above = DistanceMatrix[row-1][col]
      if col-1 > 0
        left = DistanceMatrix[row][col-1]
        diag = DistanceMatrix[row-1][col-1]
      else
        left = Inf
        diag = Inf
      end
    else
      above = Inf
      diag = Inf
      if col-1 > 0
        left = DistanceMatrix[row][col-1]
      else
        left = Inf
      end
    end
      
    if diag <= left
      if diag <= above
        s_new = s[row-1] * s_new
        w_new = w[col-1] * w_new
        row -= 1
        col -= 1
      else
        s_new = s[row-1] * s_new
        w_new = '_' * w_new
        row -= 1
      end
    else
      if left <= above
        s_new = '_' * s_new
        w_new = w[col-1] * w_new
        col -= 1
      else
        s_new = s[row-1] * s_new
        w_new = '_' * w_new
        row -= 1
      end
    end
  end
  return DistanceMatrix, dist, s_new, w_new
end
