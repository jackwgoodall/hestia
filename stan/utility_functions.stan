
functions {
  
  // Calculate diagonal element i of matrix m
  // Diagonal element = 1 - sum(non-diagonal elements)
  real get_diagonal_element(matrix m, int i){
    real out;
    out = 1;
    for(j in 1:rows(m)) {
      if(j != i) {
        out = out - m[j, i];  
      }
    }
    return out;
  }
  
  // function that is comparable to R's %in% function
  // pos is the value to test for matches
  // array pos_var is the 1-dimensional array of possible matches
  // returns pos_match=1 if pos matches at least one element of pos_var and pos_match=0 otherwise
  // example code: 
  // if(r_in(3,{1,2,3,4})) will evaluate as TRUE
  // from: https://discourse.mc-stan.org/t/stan-equivalent-of-rs-in-function/3849
  int is_in(int pos,int[] pos_var) {
    int pos_match;
    int all_matches[size(pos_var)];
    
    for (p in 1:(size(pos_var))) {
      all_matches[p] = (pos_var[p]==pos);
    }
    
    if(sum(all_matches)>0) {
      pos_match = 1;
      return pos_match;
    } else {
      pos_match = 0;
      return pos_match;
    }
  }
  
  // Normailze the columns of a matrix to sum to 1
matrix normalize_cols(matrix m) {
  matrix[rows(m), cols(m)] out;
  
  for(i in 1:cols(m)) {
    out[,i] = m[,i]/sum(m[,1]);
  }
  return out;
}
  
}
