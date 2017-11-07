library(data.tree)
library(rstack)

str_to_char <- function(string){
  chars <- strsplit(string, "")[[1]]
  remove <- which(chars == " ")
  if (length(remove) > 0){
    chars <- chars[-remove]
  }
  return(chars)
}

char_to_str <- function(chars){
  string <- paste(chars,sep = '', collapse = '')
  return(string)
}

extract_logical_from_base <- function(chars){
  and <- which(chars == "&")
  or <- which(chars == "|")
  
  if ((length(and) > 0) & (length(or) > 0)){
    print('not base')
    return()
  }
  
  if ((length(and) > 0)){
    return('&')
  }
  
  if ((length(or) > 0)){
    return('|')
  }
}

extract_elems_from_base <- function(chars){
  str <- char_to_str(chars) #paste(str,sep = '', collapse = '')
  elem <- strsplit(str, split = '\\||&')[[1]]
  remove <- which(elem == "")
  if (length(remove) > 0){
    elem <- elem[-remove]
  }
  return(elem)
}

create_base_node <- function(chars){
  # and <- which(chars == "&")
  # or <- which(chars == "|")
  # 
  # if ((length(and) > 0) & (length(or) > 0)){
  #   print('not base')
  #   return()
  # }
  # 
  logic_node <- Node$new('_')
  
  logical <- extract_logical_from_base(chars)

  logic_node$name <- logical
  elems <- extract_elems_from_base(chars) #strsplit(char_to_str(chars), logical)[[1]]
  for (i in elems){
    logic_node$AddChild(i)
  }
  
    
  # if ((length(and) > 0)){
  #   logic_node$name <- "&"
  #   elems <- strsplit(char_to_str(chars), '&')[[1]]
  #   for (i in elems){
  #     logic_node$AddChild(i)
  #   }
  # }
  # if ((length(or) > 0)){
  #   logic_node$name <- "|"
  #   elems <- strsplit(char_to_str(chars), '|')[[1]]
  #   for (i in elems){
  #     logic_node$AddChild(i)
  #   }
  # }
  
  return(logic_node)
}

create_node <- function(string){
  chars <- strsplit(string, "")[[1]]
}

string <- "(k &((d & i & j) | a) & g & ((b & (e | h)) | c) & f)" # "((a&b)|(c&d))" # yeast_model_mod@gprRules[913]
chars <- strsplit(string, "")[[1]]
# print(paste(chars,sep = '', collapse = ''))
chars <- chars[-which(chars == " ")]
# print(paste(chars,sep = '', collapse = ''))
len <- length(chars)
open <- which(chars == "(")
close <- which(chars == ")")

# for (i in 1:max(max(open), max(close))){
#   print(i)
# }
#print(which(chars == "("))
#print(which(chars == ")"))

root_node <- Node$new('root_node')
# parent <- root_node$AddChild('1')
# parent2 <- Node$new('2')
# parent$AddChildNode(parent2)
# parent2 <- Node$new('3')
# parent$AddChildNode(parent2)
# print(root_node)

# s <- stack$new()
# 
# n <- FALSE
# for (i in chars){
#   s$push(i)
#   if (i == ")"){
#     j <- s$pop()
#     j <- s$pop()
#     char <- c()
#     while (j != "("){
#       char <- c(j, char)
#       j <- s$pop()
#     }
#     #str <- c(str, j)
#     num_logical <- length(which(char == '&' | char == '|'))
#     # str <- char_to_str(char) #paste(str,sep = '', collapse = '')
#     # elem <- strsplit(str, split = '\\||&')[[1]]
#     elem <- extract_elems_from_base(char) #elem[-which(elem == "")]
#     num_elem <- length(elem)
#     print(paste(num_elem, num_logical))
#     print(char)
# 
#     # new_node <- Node$new("_")
#     if (num_logical >= num_elem){
#       temp_new <- Node$new(extract_logical_from_base(char))
#       for (i in elem){
#         temp_new$AddChild(i)
#       }
#       temp_new$AddChildNode(new_node)
#       new_node <- temp_new
#       print(temp_new)
#       print(new_node)
#       n <- TRUE
#     }
#     else {
#       # print(char)
#       new_node <- create_base_node(char)
#       if (n){
#         temp_new$AddChildNode(new_node)
#         print(temp_new)
#       }
#       else{
#         print(new_node)
#       }
#       #n <- TRUE
#       # print('smth')
#     }
#      # Node$new(str)
#     # print(new_node)
#   }
# }
i <- 443
a <- capture.output(call_tree(parse(text = yeast_model_mod@gprRules[i])))
print(yeast_model_mod@gprRules[i])
# for (i in 1:length(a)){
#   print(paste(i, ":", (which(strsplit(a[i], "")[[1]] == "\\")-1)/2, a[i]))
# }

# root_node <- Node$new('root node')
# new_node <- Node$new('_')
# for (i in 1:length(a)){
#   chars <- strsplit(a[i], "")[[1]]
#   depth <- (which(strsplit(a[i], "")[[1]] == "\\")-1)/2
#   id <- chars[length(chars)] 
#   if (id == '|' | id == '&'){
#     new_node <- Node$new(id)
#   }
#   print(j)
# }

tree_building <- function(list, index){
  
  chars <- strsplit(list[index], "")[[1]]
  og_depth <- (which(strsplit(list[index], "")[[1]] == "\\")-1)/2
  og_id <- chars[length(chars)]
  id1 <- "_" #og_id
  id2 <- "_"
  
  # print(paste('start', og_id))
  new_node <- Node$new(og_id)
  
  start_1 <- index + 1
  chars <- strsplit(list[start_1], "")[[1]]
  depth1 <- (which(chars == "\\")-1)/2
  id1 <- chars[length(chars)]
  
  start_2 <- start_1 + 1
  chars <- strsplit(list[start_2], "")[[1]]
  depth2 <- (which(chars == "\\")-1)/2
  id2 <- chars[length(chars)]
  while (depth2 != depth1 | id2 != id1){
    start_2 <- start_2 + 1
    chars <- strsplit(list[start_2], "")[[1]]
    depth2 <- (which(chars == "\\")-1)/2
    id2 <- chars[length(chars)]
  }
  
  # print(paste(start_1, start_2))
  
  # depth1 <- og_depth
  # child_ct <- 0
  index1 <- start_1
  # chars <- strsplit(list[index1], "")[[1]]
  # depth1 <- (which(chars == "\\")-1)/2
  # id1 <- chars[length(chars)]
  # print('starting id1')
  while (!(id1 %in% c( '|', '&', 'x'))){
    # print(paste(index1, id1, depth1))
    index1 <- index1 + 1
    chars <- strsplit(list[index1], "")[[1]]
    depth1 <- (which(chars == "\\")-1)/2
    id1 <- chars[length(chars)]
  }
  
  index2 <- start_2
  while (!(id2 %in% c( '|', '&', 'x'))){
    # print(paste(index1, id1, depth1))
    index2 <- index2 + 1
    chars <- strsplit(list[index2], "")[[1]]
    depth2 <- (which(chars == "\\")-1)/2
    id2 <- chars[length(chars)]
  }
  
  # print(paste('indexes',index1, index2))
  # print(paste('ids',id1, id2))
  # print('ending')
  # 
  left_node <- Node$new("_")
  right_node <- Node$new("_")
  if (id1 == 'x'){
    chars <- strsplit(list[index1+1], "")[[1]]

    len <- length(chars)
    id1 <- chars[len]
    while (id1 == " "){
      len <- len - 1
      id1 <- chars[len]
    }

    # print(paste('1 x', id1, id2))
    # new_node$AddChild(id1)
    left_node <- Node$new(id1)
  }
  else{
    # print(paste('1 node', id1, id2))
    left_node <- tree_building(list, index1)
    # print('extra node')
    # print(next_node)
    # new_node$AddChildNode(next_node)
    # print(new_node)
  }
  # 
  if (id2 == 'x'){
    chars <- strsplit(list[index2+1], "")[[1]]

    len <- length(chars)
    id2 <- chars[len]
    while (id2 == " "){
      len <- len - 1
      id2 <- chars[len]
    }

    # print(paste('2 x', id1, id2))
    right_node <- Node$new(id2)
    # new_node$AddChild(id2)
  }
  else {
    # print(paste('2 node', id1, id2))
    right_node <- tree_building(list, index2)
    # print(next_node)
    # new_node$AddChildNode(next_node)
    # print(new_node)
  }
  
  if (right_node$name == left_node$name){
    right_node$name <- paste(right_node$name, 2, sep = '')
  }
  
  new_node$AddChildNode(left_node)
  new_node$AddChildNode(right_node)
  # left_node$AddSiblingNode(right_node)
  # print('check nodes')
  # print('left')
  # print(left_node)
  # print('right')
  # print(right_node)
  # print('final')
  # print(new_node)
  return(new_node)
}

a_ <- tree_building(a, 4)
print(a_)