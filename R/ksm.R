# clean-up
rm(list = ls())

# Libraries
library(dplyr)
library(stringr)
library(seqinr)
library(data.table)
library(igraph)

# Main functions

#' Load kinase-subtrate mapping library
#'
#' Description
#'
#' @param path character
#' @export

loadLibrary <- function(path = ".") {
  setwd(path)
  .makeFolders(path)
  if (file.exists("./_library/substrateLibrary.rda")) {
    load("./_library/substrateLibrary.rda", envir = .GlobalEnv)
  } else {
    cat("WARNING: Library file does not exists.\n")
  }
  
}

#' Generate kinase-subtrate mapping library
#'
#' Description
#'
#' @param path character
#' @import data.table
#' @importFrom dplyr filter %.% rbind_all mutate group_by ungroup arrange
#' @export

generateLibrary <- function(path = ".") {
  # load uniprot sequence annotation file
  setwd(path)
  .makeFolders(path)
  if (file.exists("./_library/uniTable.rda")) {
    load("./_library/uniTable.rda")
  } else {
    uniTable <- generateUniprotMapping(path)
  }
  
  # read full data set
  files <- dir(recursive = TRUE, pattern = ".+_.+_.+_.+_.+_.+_.+_.+.txt")
  fullData <- rbind_all(lapply(files, .readMSFile))
  
  # generate data set filtered by specific phosphopeptides
  .dataFiltered <- fullData %.%
    # filter phosphopeptides
    filter(grepl("Phospho", Protein.Mods)) %.%
    # make type consistent
    mutate(type = ifelse(toupper(type) %in% c("GK", "AS"), "AS", toupper(type))) %.%
    # filter specific
    group_by(experiment, DB.Peptide) %.%
    mutate(specific = ifelse(any(type != "AS"), "unspecific", "specific")) %.%
    ungroup() %.%
    filter(specific == "specific") %.%
    # extract phosphosites
    mutate(sites = .extractPattern(Protein.Mods, "Phospho@([0-9|&]+)=*[0-9]*")) %.%
    # extract SLIP score
    mutate(slip = .extractPattern(Protein.Mods, "Phospho@[0-9|&]+=*([0-9]*)")) %.%
    # extract independent IDs
    group_by(UniProt.ID, sites, experiment) %.%
    mutate(count = min(length(sites), length(unique(sample)))) %.%
    filter(as.numeric(P.Value) == min(as.numeric(P.Value))) %.%
    ungroup()
  
  # generate data set aggregated by kinase
  .dataAggregate <- .dataFiltered %.%  
    group_by(UniProt.ID, sites, kinase) %.%
    mutate(count = sum(count)) %.%
    filter(as.numeric(P.Value) == min(as.numeric(P.Value))) %.%
    ungroup() %.%
    # get number of times peptide was found
    group_by(kinase) %.%
    mutate(tries = length(unique(paste(experiment, sample, sep = "_"))), 
           found = round(count / tries, 3)) %.%
    ungroup() %.%
    # sort by kinase and protein
    arrange(kinase, Gene)
  
  # add consensus sequence to data set
  .dataAggregate["consensus"] <- apply(.dataAggregate, 1, .extractConsensus, uniTable)
  
  # save tables
  save(.dataAggregate, .dataFiltered, file = "./_library/substrateLibrary.rda")
  
  # load tables
  loadLibrary(path)
}

#' Generate mapping of Uniprot accession to protein sequence
#'
#' Description
#'
#' @param path character
#' @import data.table
#' @importFrom seqinr read.fasta
#' @importFrom dplyr rbind_list
#' @return data.table

generateUniprotMapping <- function(path = ".") {
  setwd(path)
  files <- dir("./_uniprot/", pattern = ".fasta")
  uniTable <- data.frame()
  for (i in seq(files)) {
    file <- paste0("./_uniprot/", files[i])
    seqs <- read.fasta(file, seqtype = "AA", set.attributes = FALSE)
    accessions <- sapply(names(seqs), 
                         function(name) unlist(strsplit(name, "|", fixed = TRUE))[2])
    seqs <- sapply(seqs, paste0, collapse = "")
    d <- data.frame("accession" = accessions, "sequence" = seqs, stringsAsFactors = FALSE)
    uniTable <- rbind_list(uniTable, d)
  }
  uniTable <- data.table(uniTable)
  setkey(uniTable, accession)
  save(uniTable, file = "./_library/uniTable.rda")
  return(uniTable)
}

#' Generate summary of kinase-substrate relationships for a kinase.
#'
#' Description
#'
#' @param kin character
#' @param toFile boolean
#' @param path character
#' @param data data.table
#' @import data.table
#' @importFrom dplyr filter %.%
#' @export

getKinaseSummary <- function(kin, toFile=TRUE, path=".", data = .dataAggregate) {
  ks <- data %.% filter(kinase == kin)
  if (toFile == TRUE) {
    fn <- paste0("./_results/", kin, "_summary.txt")
    write.table(ks, file = fn, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

#' Generate summary of kinase-substrate relationships for an experiment.
#'
#' Description
#'
#' @param exp character
#' @param toFile boolean
#' @param path character
#' @param data data.table
#' @import data.table
#' @importFrom dplyr filter %.%
#' @export

getExperimentSummary <- function(exp, toFile=TRUE, path=".", data=.dataFiltered) {
  es <- data %.% filter(exp == experiment)
  if (toFile == TRUE) {
    fn <- paste0("./_results/", exp, "_summary.txt")
    write.table(es, file = fn, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

#' Generate list of local sequences around phosphorylation site. 
#'
#' Description
#'
#' @param kin character
#' @param ids integer
#' @param toFile boolean
#' @param path character
#' @param data data.table
#' @import data.table
#' @importFrom dplyr filter %.% select
#' @export

getConsensus <- function(kin, ids=1, toFile=TRUE, path=".", data=.dataAggregate) {
  cs <- data %.% filter(kinase == kin & count >= ids) %.% select(consensus)
  cs <- unlist(lapply(cs, strsplit, ";"))
  # filter out too short sequences
  cs <- cs[nchar(cs) == 13]
  # filter out non-unique sequences
  cs <- unique(cs)
  # TODO: filter out sequences not centered around phosphorylateable residue
  if (toFile == TRUE) {
    fn <- paste0("./_results/", kin, "_summary.txt")
    write.table(cs, file = fn, sep = "\n", row.names = FALSE, quote = FALSE)
  }
}

#' Generate graph of substrates
#'
#' Description
#'
#' @param kin character
#' @param org character
#' @param interaction character
#' @param go boolean
#' @param data data.table
#' @import igraph
#' @importFrom dplyr filter %.% select group_by summarise arrange
#' @export

getGraph <- function(kin, org, interaction = "physical", go = TRUE, data = .dataAggregate) {
  # map organism to full name
  organisms <- c("sp" = "Schizosaccharomyces_pombe", 
                 "hs" = "Homo_sapiens", 
                 "mm" = "Mus_musculus", 
                 "sc" = "Saccharomyces_cerevisiae")
  org <- organisms[org]
  # select substrate data
  gd <- data %.% 
    filter(kinase == kin) %.% 
    select(Gene, found) %.%
    group_by(Gene) %.% 
    summarise(weight = max(found))
  # eliminate nodes with empty gene names
  gd <- gd[nchar(gd[, 1]) > 0, ]
  substrates <- gd[, 1]
  weights <- ifelse(gd[, 2] == 1, 0.99, gd[, 2])
  # add kinase to node list if not present, assign minimum weight
  if (!kin %in% substrates) {
    substrates <- c(substrates, kin)
    weights <- c(weights, min(weights))
  }
  # TODO: replace known erroneous names with correct ones (ded1 -> sum3)
  substrates[grep("ded1", substrates)] <- "sum3"
  
  # make phosphorylation edges table
  edgesPhospho <- data.frame(A = kin, B = substrates, type = "phosphorylation")
  
  # make interaction edges table
  bgFiles <- dir("./_biogrid/", pattern = "BIOGRID-ORGANISM-")
  bgFile <- bgFiles[grep(org, bgFiles)]
  if (length(bgFile) == 1) {
    bgFile <- paste0("./_biogrid/", bgFile)
    bg <- read.table(bgFile, sep = "\t", quote = "", fill = TRUE)
    bg <- select(bg, V8, V9, V13)
    names(bg) <- c("A", "B", "type")
    edgesInteraction <- filter(bg, A %in% substrates & B %in% substrates)
    if (interaction != "all") {edgesInteraction <- filter(edgesInteraction, type == interaction)}
    edgesInteraction <- unique(edgesInteraction)
    edges <- rbind(edgesInteraction, edgesPhospho)
  } else {
    cat("WARNING: No BioGrid interactions found.\n")
    edges <- edgesPhospho
  }
  # get GO annotation attributes for substrates
  goFiles <- dir("./_go/", pattern = "DAVID")
  classFiles <- dir("./_go/", pattern = "CLASS")
  goFile <- goFiles[grep(kin, goFiles)]
  classFile <- classFiles[grep(kin, classFiles)]
  if (go == TRUE & length(goFile) == 1 & length(classFile) == 1) {
    goFile <- paste0("./_go/", goFile)
    classFile <- paste0("./_go/", classFile)
    ge <- read.table(goFile, sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)
    ge <- ge %.% filter(grepl("^GOTERM", Category)) %.%
      select(Term, Benjamini, Genes) %.%
      arrange(as.numeric(Benjamini))
    idsByGo <- lapply(split(ge, ge$Term), "[", 3)
    idsByGo <- lapply(idsByGo, as.character)
    idsByGo <- lapply(idsByGo, strsplit, ", ")
    idsByGo <- lapply(idsByGo, "[[", 1)
    # map uniprot ids to gene names
    uni2gene <- data %.% filter(kinase == kin) %.% select(Acc.., Gene) %.% unique()
    uni <- uni2gene[, "Acc.."]
    uni2gene <- as.character(uni2gene[, "Gene"])
    names(uni2gene) <- uni
    idsByGo <- lapply(idsByGo, sapply, function(ID) uni2gene[ID])
    # create GO membership matrix
    goMat <- matrix(NA, nrow = length(substrates), ncol = length(idsByGo), 
                    dimnames = list(substrates, names(idsByGo)))
    for (i in substrates) {
      for (j in seq(idsByGo)) { goMat[i, j] <- ifelse(i %in% idsByGo[[j]], 1, 0) }
    }
    classes <- read.table(classFile, sep = "\t", stringsAsFactors = FALSE)
    groups <- list()
    for (i in seq(nrow(classes))) {
      groups[[classes[i, 1]]] <- grep(classes[i, 2], names(idsByGo))
    }
    grouping <- data.frame(node = substrates, group = "Other", stringsAsFactors = FALSE)
    for (i in seq(substrates)) {
      for (j in seq(groups)) {
        if (sum(goMat[i, groups[[j]]]) > 0 ) {
          grouping[i, "group"] <- names(groups)[j]
          break
        }
      }
    }
    go <- select(grouping, group)
  } else {
    cat("WARNING: No GO annotations found.\n")
    go <- "Other"
  }
  
  # make node attribute table
  nodeAttr <- data.frame(node = substrates, 
                         class = "substrate", 
                         weight = weights, 
                         go = go,
                         stringsAsFactors = FALSE)
  nodeAttr[grep(kin, nodeAttr$node), "class"] <- "kinase"
  
  # make gml graph
  g <- graph.data.frame(edges, directed = TRUE, vertices = nodeAttr)
  g <- delete.edges(g, which(is.loop(g)))
  g <- delete.edges(g, which(is.multiple(g)))
  g <- delete.vertices(g, which(igraph::degree(g) < 2 & nodeAttr$group == "Other"))
  write.graph(g, file = paste0("./_results/", kin, "_", interaction, ".gml"), format = "gml")
}

#' Get list of experiments
#'
#' Description
#'
#' @param data data.table
#' @export

experiments <- function(data = .dataAggregate) {
  unique(data$experiment)
}


#' Get list of kinases
#'
#' Description
#'
#' @param data data.table
#' @export

kinases <- function(data = .dataAggregate) {
  unique(data$kinase)
}

#' Summarize all experiments
#'
#' Description
#'
#' @param data data.table
#' @export

summarizeExperiments <- function() {
  for (e in experiments()) { getExperimentSummary(exp = e) }
}

#' Summarize all experiments
#'
#' Description
#'
#' @param data data.table
#' @export

summarizeKinases <- function() {
  for (k in kinases()) { getKinaseSummary(kin = k) }
}

# Utility functions

#' Read mass spectrometry file.
#'
#' Description
#'
#' @param path character
#' @import stringr

.readMSFile <- function(path) {
  rawSubstrates <- read.delim(path, skip = 2, fill = TRUE, stringsAsFactors = FALSE,
                              colClasses = "character")
  expInfo <- str_replace(basename(path), ".txt", "")
  expInfo <- unlist(str_split(expInfo, "_"))
  expNames <- c("date", "experiment", "kinase", "condition", "type", 
                "method", "sample", "scientist")
  rawSubstrates[, expNames] <- rep(expInfo, each = nrow(rawSubstrates))
  rawSubstrates
}

#' Extract pattern from raw modification string.
#'
#' Description
#'
#' @param column character
#' @param pattern character

.extractPattern <- function(column, pattern) {
  matches <- str_match_all(column, pattern)
  matches <- sapply(matches, function(m) paste(m[, 2], collapse = ";"))
}

#' Extract consensus sequence for modification from protein sequence.
#'
#' Description
#'
#' @param row character
#' @param uniTable data.table

.extractConsensus <- function(row, uniTable) {
  sites <- row["sites"]
  species <- row["Species"]
  accession <- row["Acc.."]
  seq <- uniTable[accession, "sequence", with = FALSE]
  seq <- as.character(seq)
  sites <- unlist(strsplit(sites, ";", fixed = TRUE))
  ambiguous <- grepl("|", sites, fixed = TRUE)
  sites[ambiguous == TRUE] <- NA
  sites <- as.integer(sites)
  consensus <- str_sub(seq, (sites - 6), (sites + 6))
  badConsensus <- is.na(consensus) | !(str_sub(consensus, 7, 7) %in% c("S", "T", "Y"))
  consensus[badConsensus] <- ""
  consensus <- paste(consensus, collapse = ";")
}

#' Create neccessary folders.
#'
#' Description
#'
#' @param path character

.makeFolders <- function(path = ".") {
  if (!file.exists("./_library/")) {dir.create("./_library/", showWarnings = FALSE)}
  if (!file.exists("./_results/")) {dir.create("./_results/", showWarnings = FALSE)}
  if (!file.exists("./_go/")) {dir.create("./_go/", showWarnings = FALSE)}
  if (!file.exists("./_biogrid/")) {dir.create("./_biogrid/", showWarnings = FALSE)}
  if (!file.exists("./_uniprot/")) {dir.create("./_uniprot/", showWarnings = FALSE)}
}