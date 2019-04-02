# GGPLOT functions #############################################################

# Multiplot for ggplot2 by winston@stdout.org from Cookbook for R
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    # library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}

# Boxplot shortcut function (ggplot2)
ggBoxplot <- function(matrix, title = "", xlab = "x", ylab = "y", outLCol = NA){
    ggplot(data=melt(as.data.frame(matrix)), aes(variable, value)) + 
        geom_boxplot(outlier.colour= outLCol, outlier.size = 1) + xlab(xlab) + 
        ylab(ylab) + ggtitle(title) + theme_classic() + 
        stat_n_text(size = 3, angle = 90) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Pie Chart Nucleotide frequency
ggPieFreq <- function(Freq, labSize = 5){
    palette(brewer.pal(9, "Set1")) #RColorBrewer
    tmpDF <- data.frame(nucs = names(Freq), Percent = Freq, stringsAsFactors = F)
    tmpDF <- tmpDF[order(tmpDF[,2], decreasing = T),]
    tmpDF <- data.frame(rbind(tmpDF[1:5,], c("All Others", sum(tmpDF[-1:-5,2]))))
    tmpDF[,2] <- as.numeric(tmpDF[,2]) / sum(as.numeric(tmpDF[,2]))
    tmpDF[,1] <- factor(tmpDF[,1], levels = tmpDF[,1])
    ggPie <- ggplot(tmpDF, aes(x="", y=Percent, fill=nucs)) +
        geom_bar(width = 1, stat = "identity") + 
        coord_polar(theta = "y",start = 0,direction = 1) + 
        geom_text(aes(label = round(Percent,2)), size= labSize, position = position_stack(vjust = 0.5)) +
        theme(axis.text.x =element_blank()) + theme_classic()
    return(ggPie)
}

# GGplot alternative to pairs function (additionally it fits linear models to all pair-wise comparisons)
ggPairs <- function(DF, alpha = 1){
    iCol <- colnames(DF)
    matD <- combinations(n = length(iCol), r = 2, v = 1:length(iCol))  
    ggSC <- lapply(1:nrow(matD), function(x){
        tmpL <- lm(DF[,matD[x,2]] ~ DF[,matD[x,1]])
        if(tmpL$coefficients[1]>=0){
            linModEq = paste("y = x *",tmpL$coefficients[2] %>% signif(2),
             "+", tmpL$coefficients[1]  %>% signif(2))
        }else if(tmpL$coefficients[1]<0){linModEq = paste("y = x *",
         signif(tmpL$coefficients[2],2), "-", tmpL$coefficients[1] %>% 
         signif(2) %>% abs)}
        tmpC <- cor(DF[,matD[x,1]], DF[,matD[x,2]], use = "p") %>% round(4)
        tmpP <- cor.test(DF[,matD[x,1]], DF[,matD[x,2]], use = "p")$p.value %>% signif(4)
        tmpC2 <- cor(DF[,matD[x,1]], DF[,matD[x,2]], use = "p", method = "spearman") %>% round(4)
        tmpP2 <- cor.test(DF[,matD[x,1]], DF[,matD[x,2]], use = "p", method = "spearman")$p.value %>% signif(4)
        ggplot(DF, aes(x= DF[,matD[x,1]], y= DF[,matD[x,2]])) + 
            geom_point(alpha = alpha, shape = 16) + 
            geom_smooth(method = lm, se=FALSE, fullrange= TRUE, aes(group=1)) +
            geom_abline(intercept = 0, slope = 1, colour = "gray") + 
            theme_classic() + xlab(iCol[matD[x,1]]) + ylab(iCol[matD[x,2]]) + 
            ggtitle(paste("R =", tmpC, "p = ", tmpP, "\nrho =", tmpC2, "p =", tmpP2), 
                    subtitle = linModEq) +
            coord_cartesian(ylim = range(DF, na.rm = T), xlim = range(DF, na.rm = T))
    })
    ggLabs <- lapply(iCol, function(x){
        df <- data.frame(x = 1, y = 1, text = x)
        ggO <- ggplot(df, aes(x, y)) +
            geom_text(aes(label = text), size = 5) + theme_classic() +
            theme(panel.border = element_rect(colour = 1, fill = NA), axis.line = element_line())+
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) + 
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())
        return(ggO)
    })  
    ggCf <- lapply(1:nrow(matD), function(x){
        return()
    })  
    lOut <- matrix(NA, ncol = ncol(DF), nrow = ncol(DF))
    for(i in 1:nrow(matD)){lOut[matD[i,2], matD[i,1]] <- i}
    for(i in 1:length(iCol)){lOut[i, i] <- length(ggSC) + i}
    for(i in 1:nrow(matD)){lOut[matD[i,1], matD[i,2]] <- length(ggSC) + length(iCol) + i}
    multiplot(plotlist = c(ggSC, ggLabs), layout = lOut)
}

# Simple Barplot function
ggBarplot <- function(x, ci = NA, title = NULL, subt = NULL, xLab = "Names", yLab = "Values"){
    if(is.null(names(x))){names(x) <- 1:length(x)}
    df <- data.frame(names = names(x), value=x, CI = ci)
    outGG <- ggplot(data=df, aes(x=names, y=value)) + 
        geom_bar(stat="identity") + theme_classic() + 
        geom_errorbar(aes(ymin=value-CI, ymax=value+CI), width=.2, position=position_dodge(.9)) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(title, subt) +
        ylab(yLab) + xlab(xLab)
    return(outGG)
}

# Scatterplot with linear model fitted line and calculates correlation
ggScattLinePlot <- function(x, y, title = "", xLab = "", yLab = "", alpha = 1){
    tmpC <- cor(x, y, use = "p") %>% round(4)
    tmpP <- cor.test(x, y, use = "p")$p.value %>% signif(3)
    tmpC2 <- cor(x, y, use = "p", method = "spearman") %>% round(4)
    tmpP2 <- cor.test(x, y, use = "p", method = "spearman")$p.value %>% signif(3)
    tmpDF <- data.frame(var1 = x, var2 = y)
    ggSCLINE <-  ggplot(tmpDF, aes(x = var1, y = var2)) + geom_point(alpha = alpha) + 
        geom_smooth(method = lm, se=FALSE, fullrange= TRUE, aes(group=1)) +
        ggtitle(title, paste("R =", tmpC, "p = ", tmpP, "\nrho =", tmpC2, "p =", tmpP2)) +
        ylab(yLab) + xlab(xLab) + theme_classic()
    return(ggSCLINE)
}

# ggplot heatmap
ggHeatmap <- function(x, y, logTransform = T, nBins = 100){
    tmpC <- cor(x, y, use = "p") %>% round(4)
    tmpP <- cor.test(x, y, use = "p")$p.value %>% signif(3)
    tmpC2 <- cor(x, y, use = "p", method = "spearman") %>% round(4)
    tmpP2 <- cor.test(x, y, use = "p", method = "spearman")$p.value %>% signif(3)
    tmpDF <- data.frame(var1 = x, var2 = y)
    if(logTransform == T){
        ggplot(tmpDF, aes(x = var1, y = var2)) + geom_bin2d(bins = nBins) +
            scale_fill_gradientn(trans = "log", colours = rev(brewer.pal(9, "Spectral"))) +
            theme_classic() + ggtitle("", paste("R =", tmpC, "p = ", tmpP, "\nrho =", tmpC2, "p =", tmpP2)) +
            geom_smooth(method = lm, se=FALSE, fullrange= TRUE, aes(group=1)) + theme(legend.position = "bottom")
    }else{
        ggplot(tmpDF, aes(x = var1, y = var2)) + geom_bin2d(bins = nBins) +
            scale_fill_gradientn(colours = rev(brewer.pal(9, "Spectral"))) +
            theme_classic() + ggtitle("", paste("R =", tmpC, "p = ", tmpP, "\nrho =", tmpC2, "p =", tmpP2)) +
            geom_smooth(method = lm, se=FALSE, fullrange= TRUE, aes(group=1)) + theme(legend.position = "bottom")
    }
    
}

# ggboxplot with variables cut into ordered categories by Interval
ggBoxPlot_cutInterval <- function(x, y, nCut_x){
    tmpDF <- data.frame(var1 = x, var2 = y)
    tmpDF$cut_x <- cut_interval(tmpDF$var1, nCut_x)
    ggplot(tmpDF, aes(y= var2, x = cut_x)) + geom_boxplot(outlier.colour = NA) + theme_classic() +
        stat_n_text(size = 3, angle = 90) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# ggboxplot with variables cut into ordered categories by number of data points
ggBoxPlot_cutNumber <- function(x, y, nCut_x){
    tmpDF <- data.frame(var1 = x, var2 = y)
    tmpDF$cut_x <- cut_number(tmpDF$var1, nCut_x)
    ggplot(tmpDF, aes(y= var2, x = cut_x)) + geom_boxplot(outlier.colour = NA) + theme_classic() +
        stat_n_text(size = 3, angle = 90) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Load/install dependencies ####################################################
# Load/install CRAN packages
installLoad_CRAN <- function(package){
    if (!require(package, character.only = T)) {
        install.packages(package, dependencies = TRUE)
        library(package, character.only = T, quietly = T)
    }
}

CRAN_packs <- c("magrittr", "plyr", "ggplot2", "grid", "gtools", "reshape2", "EnvStats")
sapply(CRAN_packs, function(x) installLoad_CRAN(x))
