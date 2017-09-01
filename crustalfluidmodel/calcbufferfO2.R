# Calculate fO2 buffered by mineral assemblages.  Calculations at T, P used in
# CRUSTALFLUIDMODEL.R.  

TP = expand.grid(T = unique(conds$T), P = unique(conds$P))
buffernames = c("FMQ", "PPM", "MH", "NNO", "IW", "WM", "IM", "QIF", "FFM")

buffO2 = array(numeric(), dim = c(dim(TP)[1],length(buffernames)))
colnames(buffO2) = buffernames
buffO2 <- data.frame(buffO2)

buffO2$FMQ <- -subcrt(c("O2", "fayalite","magnetite","quartz"),c(-1,-3,2,3),c("g","cr","cr","cr"), T=TP$T, P=TP$P)$out$logK
buffO2$PPM <- -subcrt(c("O2", "pyrrhotite", "pyrite", "magnetite"), c(-1, -3, 1.5, 0.5), c("g","cr","cr","cr"), T=TP$T, P=TP$P)$out$logK
buffO2$MH  <- -subcrt(c("O2", "magnetite","hematite"), c(-1, -4, 6), c("g", "cr", "cr"), T=TP$T, P=TP$P)$out$logK
buffO2$NNO <- -subcrt(c("O2", "Ni","NiO"),c("g", "cr","cr"),c(-1,-2,2), T=TP$T, P=TP$P)$out$logK
buffO2$IW  <- -subcrt(c("O2", "Fe", "FeO"), c("g", "cr", "cr"), c(-1, -2, 2), T=TP$T, P=TP$P)$out$logK
buffO2$WM  <- -subcrt(c("O2", "FeO", "magnetite"), c("g", "cr", "cr"), c(-1, -6, 2), T=TP$T, P=TP$P)$out$logK
buffO2$IM  <- -subcrt(c("O2", "Fe", "magnetite"), c("g", "cr", "cr"), c(-1, -3/2, 1/2), T=TP$T, P=TP$P)$out$logK
buffO2$QIF <- -subcrt(c("O2", "quartz", "Fe", "fayalite"), c("g", "cr", "cr", "cr"), c(-1, -1, -2, 1), T=TP$T, P=TP$P)$out$logK
buffO2$FFM <- -subcrt(c("O2", "fayalite","ferrosilite","magnetite"), c(-1,-6,6,2), c("g","cr","cr","cr"), T=TP$T, P=TP$P)$out$logK

buffO2$T = TP$T
buffO2$P = TP$P

buffO2 <- cbind(buffO2[-(1:length(buffernames))], buffO2[1:length(buffernames)])

write.csv(x = buffO2, file = "buffers.csv")

