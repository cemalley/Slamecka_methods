
# Project: Establishment of trophoblast stem cells from primed human pluripotent stem cells
# Analysis by: Jaroslav Slamecka, Claire Malley


# enrichR ----

# function to run enrichR
# draw plots
# save them to file
enrichr.draw.plots = function(cell.type = "TE",
                                   enrichr.databases = enrichr.DB,
                                   n.terms = 5,
                                   plot.width = 3.8
                                   ) {
  library(enrichR)
  
  # run enrichR
  enrichr.list = enrichr(markers[[cell.type]], databases=enrichr.databases)
  
  # function to plot the result of enrichR run on a single database
  single.enrichr.plot = function(enrichr.list = enrichr.list,
                               database
                               ) {
    library(extrafont)
    library(ggplot2)
    df.plot = enrichr.list[[database]] # choose a database result from the list
    df.plot = df.plot[1:n.terms, ] # keep only the top n terms
    df.plot = df.plot[seq(nrow(df.plot),1), ] # reverse order of rows
    df.plot$Adjusted.P.value = log(df.plot$Adjusted.P.value) # log transform P values
    # add gene counts to the Terms
    # the variable Overlap contains gene counts of those that belong to the term
    # the second number of the Overlap variable is total counts of genes belonging to the term
    # keep only gene counts out of submitted counts, discard overall gene counts, they are irrelevant
    # strsplit splits a string matching "split" argument
    # however, it will return a list
    # to extract only the first element of each list element vector, use sapply
    # at the moment, I don't know why this works
    gene.count.list = strsplit(df.plot$Overlap, split="/")
    df.plot$GeneCount = sapply(gene.count.list, `[[`, 1)
    df.plot$GeneCount = paste0(df.plot$GeneCount,
                               "/",
                               length(markers[[cell.type]])-length(grep(markers[[cell.type]], pattern="ENSG"))
                               )
    df.plot$Term = paste0(df.plot$Term," (",df.plot$GeneCount, ")")
    df.plot$Term = factor(df.plot$Term, levels=unique(df.plot$Term)) # convert Term to factor
    database = gsub(database, pattern="_", replacement=" ")
    database = gsub(database, pattern="-", replacement=" ")
    plot = ggplot(df.plot, aes(x=Term,
                               y=Adjusted.P.value, # Adjusted.P.value or Combined.Score, if Combined.Score is used, comment scale_y_reverse
                               fill=Adjusted.P.value)) + # the same variable has to be assigned to a fill variable as well in order for gradient scaling of the fill to work
      geom_col(show.legend=FALSE) + # if fill="darkcyan" is declared then gradient-scaling in the next line doesn't override it
      scale_fill_gradient(low="coral2", high="firebrick") + # darkcyan # flip if using Combined.Score as a metric
      xlab(paste(database, "Term", sep="\n")) +
      ylab("Adjusted P Value") +
      coord_flip() + # flip coordinate system to make x vertical
      scale_y_reverse() + # reverse y axis # comment this if Combined.Score is used as a metric
      theme_minimal(base_size=12, base_family="Noto Sans Cond") +
      theme(axis.title.y=element_text(size=10), axis.title.x=element_blank())
    return(plot)
  }
  
  # apply the function to all databases
  enrichr.plots = lapply(enrichr.databases, FUN=function(x) {single.enrichr.plot(enrichr.list, x) })
  
  # combine the results together
  library(patchwork)
  wrap_plots(enrichr.plots, ncol=1)
  ggsave(filename=paste0("enrichr.terms.",cell.type,".top",n.terms,"terms.",length(markers[[cell.type]]),"features.log",CO$logFC[[cell.type]],".p",CO$adjPval$global,".pdf"), # .png
         device=cairo_pdf, # to properly embed font
         width=plot.width,
         height=length(enrichr.plots)*(n.terms/4.4)
         ) 
}

# requirements for the function to work:
# markers - a list of markers in categories: markers[[cell.type]]
# enrichr.DB - declared below

# choose databases
enrichr.DB = c("ARCHS4_Tissues","Human_Gene_Atlas","Jensen_TISSUES","ARCHS4_Cell-lines")

# execute
enrichr.draw.plots(cell.type="TE", enrichr.DB, n.terms=5, plot.width=3.6)




