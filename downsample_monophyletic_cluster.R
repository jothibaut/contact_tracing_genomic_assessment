# @Author: Samuel L. Hong

library(ape)
library(glue)
library(dplyr)
library(readr)

# Gets loc from name of strain.
get_loc <- function(string,delimiter,index,reverse=FALSE){
    if(reverse==TRUE){
        splitted = strsplit(string,split=delimiter,fixed=T)
        reved = lapply(splitted,rev)
        return(lapply(reved,"[[",index))
    }
    else{
    return(lapply(strsplit(string,split=delimiter,fixed=T),"[[",index))
    }
}


# This function returns the location of sequence n.
get_loc_from_meta <- function(n,metadf){
    loc<- (metadf %>% filter(name==n))$loc
    return(as.character(loc))
}

# This function gets the location trait for each sequence in name_array.
# The metadata file must have a 'loc' and a 'name' column.
get_loc_array_meta <- function(name_array,meta){
    locs<-sapply(name_array,get_loc_from_meta,metadf=meta)
    return(as.character(locs))
}

# Returns the two direct descendant nodes.
getDescendants<-function(tree,node){
    if(node<=length(tree$tip.label)){
        return(-1)
    } else{
        return( tree$edge[which(tree$edge[,1]==node),2] )
    }
}

getAllDescendants<-function(tree,node,curr=NULL){
    if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
    if(is.null(curr)) curr<-vector()
    daughters<-tree$edge[which(tree$edge[,1]==node),2]
    curr<-c(curr,daughters)
    if(length(curr)==0&&node<=Ntip(tree)) curr<-node
    w<-which(daughters>Ntip(tree))
    if(length(w)>0) for(i in 1:length(w)) 
        curr<-getAllDescendants(tree,daughters[w[i]],curr)
    return(curr)
}

getTipsFromNode<-function(tree,node){
    desc <- getAllDescendants(tree,node,curr=NULL)
    return(desc[desc<=length(tree$tip.label)]) # MISTAKE: CHANGE FROM INITIAL CODE: < to <=. Since, the last tip has the number length(tree$tip.label). (Iteration starts at 1, not 0).
}


getNodeAnnots <- function(tree,tip_annot){
    tree<-reorder(tree,"postorder") # For each internal node, we will first iterate descendants.
    # prev: "pruningwise". According to reorder documentation, "postorder" is more efficient.
    num_nodes <- (length(tree$tip.label)+tree$Nnode)
    node_annot_list <- rep(0, length = num_nodes)
    rootnode <- length(tree$tip.label)+1
    for(node in c(as.array(tree$edge[,2]),rootnode)){ #Iterate over all nodes
        
        if(node<=length(tree$tip.label)){ # If node is tip, annot is loc.
            node_annot_list[node] <- tip_annot[node]
        }
        
        else{
            children_nodes <- getDescendants(tree,node)
            locs_in_children <- c()
            for (c in children_nodes){
                locs_in_children <- c(locs_in_children, node_annot_list[c])
            }

            locs_in_children <- unique(unlist(locs_in_children))
            
            # If all descendant tips have the same location, label node with this location, MIXED otherwise.
            if( ("MIXED" %in% locs_in_children) | (length(locs_in_children)>1) ){
                node_annot_list[node] = "MIXED"
            }
            else{
                node_annot_list[node] = node_annot_list[c]
            }
        }
    
    }
    return(node_annot_list)
}


# Returns the biggest non-mixed clades.
getClusters <- function(tree, node_annot, loc_to_keep){
    non_mixed <- which(node_annot!="MIXED" & node_annot!=loc_to_keep)
    children_nodes <- c() # This vector will contain each non-mixed node, which is child of another non-mixed node.
    for( internal_node in non_mixed[non_mixed>length(tree$tip.label)] ){
        children_nodes <- c(children_nodes, getDescendants(tree,internal_node))
    }
    return(setdiff(non_mixed[non_mixed>length(tree$tip.label)], children_nodes)) 
    # Only return non-mixed node which are not children --> Returns the big non-mixed clade and not the corresponding descending subclades.
}

# This function creates a df: [taxa, cluster, remove]
# taxa: tip's name
# cluster: Parent node of the biggest homogeneous cluster the tip belongs to.
# remove [bool]: Do we remove this tip? Based on a random selection.
pruneClusters<-function(tree,cluster_nodes){
    toremove_df <- as.data.frame(tree$tip.label)
    colnames(toremove_df) <- c("taxa")
    toremove_df$cluster <- 0
    toremove_df$remove <- 0
    for (c in cluster_nodes){
        
        subtree_tips <- tree$tip.label[getTipsFromNode(tree, c)] # Here again, we already go through all tips when going backward --> There's probably a way to save them in advance.
        filtered <- sample(subtree_tips, length(subtree_tips)-1) # Randomly select all tips but one.
        toremove_df$cluster[which(toremove_df$taxa %in% subtree_tips )] = c
        toremove_df$remove[which(toremove_df$taxa %in% filtered )] = 1
    }
    return(toremove_df)
}


# keep only N per cluster
update_metadata <- function(meta_df, clust_df){
  df <- left_join(meta_df, select(clust_df, taxa, remove), by = c('name'='taxa'))
  df <- df[df$remove == 0, ] %>% select(-remove)
  return(df)
}


# This function takes a large tree, with its corresponding metadata as an input.
# It then collapses all maximum clades that only contain 'loc_to_keep' sequences and reduced them to one sequence.
# It outputs the selected sequences within a new metadata file.

down_sample_tree <- function(tree_file, metadata_file, loc_to_keep, metadata_subsampled_file){
  start_time <- Sys.time()
  tree <- read.tree(tree_file)
  metadata <- read.table(file = metadata_file, sep = '\t', header=TRUE)
  
  num_taxa <- length(tree$tip.label)
  tip_annot <- get_loc_array_meta(tree$tip.label, metadata)
  annotated_nodes<-getNodeAnnots(tree, tip_annot)
  clusters<-getClusters(tree, annotated_nodes, loc_to_keep)
  
  print(glue("Starting loop after {Sys.time() - start_time}"))
  
  clustered_df <- pruneClusters(tree, clusters)
  
  num_removed <- length(clustered_df[which(clustered_df$remove==1),'taxa'])
  print(glue("Tree pruned from {num_taxa} taxa to {num_taxa-num_removed} after collapsing {length(clusters)} clusters"))
  
  keep<-(clustered_df %>% filter(remove==0))$taxa
  ktre<-keep.tip(tree,keep)
  print(glue("The tree now contains {length(ktre$tip.label)} tips."))
  
  new_metadata <- update_metadata(metadata, clustered_df)
  write.table(new_metadata, metadata_subsampled_file, sep = '\t', row.names=FALSE, quote=FALSE)
  print(Sys.time() - start_time)
}



# Down-sampling BA.1 and BA.2 Fasttrees
down_sample_tree(tree_file=paste0( "analyses/", "trees/", "ba1_fasttree.nwk"), 
                 metadata_file=paste0("data/", "genomics/", "ba1_all.tsv"), 
                 loc_to_keep="leuven",
                 metadata_subsampled_file=paste0("data/", "genomics/", "ba1_subsampled.tsv"))

down_sample_tree(tree_file=paste0( "analyses/", "trees/", "ba2_fasttree.nwk"), 
                 metadata_file=paste0("data/", "genomics/", "ba2_all.tsv"), 
                 loc_to_keep="leuven",
                 metadata_subsampled_file=paste0("data/", "genomics/", "ba2_subsampled.tsv"))
