# This function collapses the cells of a contingency table according to a column 
collapse_ftable <- function(x, col, min.collapse=1) {
  
  tmp <- x
  
  # check the locations of the cells in the selected column that have frequencies less than the minimum criterion
  loc_less <- which(tmp[, col] < min.collapse)
  
  # collapse cells
  while(length(loc_less) > 0) {
    
    # check the last row number in the contingency table
    last.num <- nrow(tmp)
    
    # check the center point 
    center <- round(mean(1:last.num) + 0.01, 0)
    
    # relocate the location of rows that have frequencies less than minimum criterion
    # "loc_less_low" is the rows in which location is less than or equal to the center point 
    loc_less_low <- loc_less[loc_less <= center]
    
    # "loc_less_high" is the rows in which location is greater than the center point 
    loc_less_high <- rev(loc_less[loc_less > center])
    
    # relocation of the rows
    loc_less <- c(loc_less_low, loc_less_high)
    
    # check the location of the first selected cell
    start.num <- loc_less[1]
    
    # when the location of the first selected cell is 1
    if(start.num == 1) {
      row.1 <- tmp[start.num,]
      row.2 <- tmp[start.num+1,]
      row.sum <- row.1 + row.2
      row.other <- tmp[-c(start.num, start.num+1), ]
      tmp <- rbind(row.sum, row.other)
    }
    
    # when the location of the first selected cell is between 1 and last row number of the contingency table
    if(start.num >=2 & start.num < last.num) {
      dif.num <- start.num - center
      sel <- ifelse(dif.num > 0, 1, 2)
      row.1 <- tmp[start.num, ]
      row.2 <- rbind(tmp[start.num-1,], tmp[start.num+1,])[sel, ]
      row.sum <- row.1 + row.2
      collapsed.nums <- sort(c(start.num, c(start.num-1, start.num+1)[sel]))
      row.other1 <- tmp[-c(collapsed.nums[1]:last.num), ]
      row.other2 <- tmp[-c(1:collapsed.nums[2]), ]
      tmp <- rbind(row.other1, row.sum, row.other2)
      
    }
    
    # when the location of the first selected cell is the last row number of the contingency table
    if(start.num == last.num) {
      row.1 <- tmp[start.num, ]
      row.2 <- tmp[start.num-1, ]
      row.sum <- row.1 + row.2
      row.other <- tmp[-c(start.num, start.num-1), ]
      tmp <- rbind(row.other, row.sum)
    }
    
    # after collapsing cells, 
    # recheck the locations of the cells in the selected column that have frequencies
    # less than the minimum criterion
    loc_less <- which(round(tmp[, col], 15) < min.collapse)
  }
  
  tmp
  
}
