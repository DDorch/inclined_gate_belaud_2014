trim.trailing <- function (x) sub("\\s+$", "", x)

trace_test_results <- function (sTest, lVersions) {

  versions = names(lVersions)
  varsToDraw = c("Q", "Cd", "Cc", "h1", "nbIter")
  
  bFirst = T
  
  
  l = list()
  bFirst = T
  for (version in versions) {
    for (modele in lVersions[[version]]$modeles) {
      dfi = trace_test_results.read.data(version, sTest, modele, lVersions[[version]]$widths, lVersions[[version]]$colnames)
      for(varToDraw in varsToDraw) {
        if(bFirst) {
          l[[varToDraw]] = data.frame(matrix(ncol = nrow(dfi) + 1, nrow = 0))
          colnames(l[[varToDraw]]) <- c("Version Modele", dfi$expe)
        }
        df = l[[varToDraw]]
        df[nrow(df) + 1, "Version Modele"] = trim.trailing(paste(version, modele))
        for (expe in dfi$expe) {
            df[nrow(df), expe] = dfi[dfi$expe == expe, varToDraw]
        }
        l[[varToDraw]] = df
      }
      bFirst = F
    }
  }
  
  for(varToDraw in varsToDraw) {
    df = l[[varToDraw]]
    df$"Version Modele" = NULL
    m = as.matrix(df)
    
    png(paste0(paste(sTest, varToDraw, sep = "_"), ".png"), width = 800, height = 600)
    par(mar=c(5.1, 14 ,4.1 ,2.1))
    barplot(m,
            main = paste("Test", sTest, ": Versions comparison on", varToDraw), 
            xlab = paste("Value of", varToDraw),
            # ylab = "Experiment",
            horiz = TRUE, 
            beside = TRUE, 
            las=2,
            col = rainbow(nrow(m))
    )
    legend("topright", l[[varToDraw]]$"Version Modele", fill= rainbow(nrow(m)))
    dev.off()
  }
  return(l)
}

trace_test_results.read.data <- function (version, sTest, modele, widths, colnames) {
  sFileName = paste("results_calcul_Q", version, sTest, sep = "_")
  if(modele != "") {sFileName = paste(sFileName, modele, sep = "_")}
  sFileName = paste0(sFileName, ".txt")
  df = read.fwf(sFileName, widths = widths, col.names = colnames, header = FALSE, skip = 2)
  df$expe=  paste0("h0=", df$h0, " h2=", df$h2, " W=", df$W, " angle=", df$angle)
  return (df)
}

lVersions = list(
  v20170601 = list(
    widths = c(6, 7, 7, 6, 9, 2, 7, 7, 7, 4), 
    colnames = c("h0", "h2", "W", "angle", "Q", "R", "Cd", "Cc", "h1", "nbIter"),
    modeles = c(""), 
    colors = c("")
  ), 
  v20181031 = list(
    widths = c(6, 7, 7, 6, 9, 2, 7, 7, 7, 7, 7, 7, 4), 
    colnames = c("h0", "h2", "W", "angle", "Q", "R", "Cd", "Cc", "h1", "dQdW", "dQdh0", "dQdh2","nbIter"),
    modeles = c("formule", "tabul")
  )
)

l1 = trace_test_results("data_test", lVersions)
l2 = trace_test_results("test45deg", lVersions)
