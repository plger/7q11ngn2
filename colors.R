WBScols = c(WBS="#ff2600", AtWBS="darkorange1", CTRL="#c0c0c0", DUP="#0433ff", CTL="#c0c0c0", "7Dup"="#0433ff")
layercols = c(RNA="#FF9896", RPF="#98DF8A", TE="#BCBD22", Protein="#9467BD")
systemcols <- c(isogenic="#4477AA", patients="#CC6677", patientDerived="#CC6677", "patient-derived"="#CC6677", "patient-\nderived"="#CC6677")

# defaults for SEtools
annocolors <- list( genotype=WBScols, layer=layercols,
            system=systemcols,
            SFARI=c("3"="lightblue", "2"="blue", "1"="darkblue"))
options("SEtools_def_anno_colors"=annocolors)
setSechmOption("anno_colors", value=annocolors)

options("SEtools_def_hmcols"=c("purple", "black", "orange"))
setSechmOption("hmcols", value=c("purple", "black", "orange"))

ggheatcol <- ggplot2::scale_color_gradient2(low="purple", mid="black", high="orange")

qualcols <- c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#C7C7C7", "#BCBD22", "#17BECF")

theme_set(theme_cowplot(12))
theme_replace(text=element_text(family="Arial", face="bold", colour="black"))
suppressMessages(extrafont::loadfonts())