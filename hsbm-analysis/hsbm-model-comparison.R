library(tidyverse)
library(reticulate)

fit <- py_load_object("fit-dict2.pickle", pickle="pickle")
fit_non_dc <- py_load_object("fit-dict-non-dc2.pickle", pickle="pickle")

dls <- bind_rows(
    tibble(model="corr", dls=fit$dls, iter=seq_along(fit$dls)),
    tibble(model="non-corr", dls=fit_non_dc$dls, iter=seq_along(fit_non_dc$dls))
)

gg <- ggplot(dls, aes(iter, dls, col=factor(model, labels=c("Degree-corrected", "Non-corrected")))) +
    geom_line(size=1.4, alpha=0.9) +
    labs(x="Iterations", y="Minimum Description Length $\\Sigma$", col=NULL) +
    scale_color_manual(values=c(`Degree-corrected`="#eb9e33", `Non-corrected`="#6f9dcf")) +
    theme_bw() +
    theme(legend.position=c(0.78, 0.85), legend.background=element_blank())

tikz_plot(gg, "model-comparison", w=4, h=4)
