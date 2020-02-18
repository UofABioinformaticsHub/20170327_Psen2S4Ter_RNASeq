library(tidyverse)
library(mvtnorm)
library(harmonicmeanp)

S <- hmFry %>%
    dplyr::select(gs_name, fry = PValue) %>%
    left_join(
        dplyr::select(hmCamera, gs_name, camera = PValue)
    ) %>%
    left_join(
        dplyr::select(hmGsea, gs_name, gsea = PValue)
    ) %>%
    left_join(
        dplyr::select(hmRiboGoseq, gs_name, goseq = PValue)
    ) %>%
    as.data.frame() %>% 
    column_to_rownames("gs_name") %>%
    as.matrix() %>% 
    qnorm() %>% 
    cov(use = "pairwise.complete.obs")

hmFry %>%
    dplyr::select(gs_name, fry = PValue) %>%
    left_join(
        dplyr::select(hmCamera, gs_name, camera = PValue)
    ) %>%
    left_join(
        dplyr::select(hmGsea, gs_name, gsea = PValue)
    ) %>%
    left_join(
        dplyr::select(hmRiboGoseq, gs_name, goseq = PValue)
    ) %>%
    as.data.frame() %>% 
    column_to_rownames("gs_name") %>% 
    as.matrix() %>% 
    qnorm() %>% 
    apply(MARGIN = 1, FUN = function(x){
        pmvnorm(upper = x[!is.na(x)], sigma = S[!is.na(x), !is.na(x)])
        }) %>% 
    as.data.frame() %>% 
    set_names("p") %>% 
    rownames_to_column("gs_name") %>% 
    as_tibble() %>% 
    arrange(p) %>% 
    mutate(
        FDR = p.adjust(p, "fdr"),
        adjP = p.adjust(p, "bonf")
        ) %>%
    dplyr::filter(FDR < 0.05)

## Alternatively, use Kost & McDermott
## This is not strictly correct, as the correlations really need to be based
## on the underlying variables, which are assumed to be N(mu, sigma)
## Here we are simply transforming the p-values using qnorm.
## As it stands, this method is more conservative than the mvnorm above
## This is all based on 
## https://www.sciencedirect.com/science/article/abs/pii/S0167715202003103
kostP <- function(p, cors){
    i <- !is.na(p)
    k <- sum(i)
    cors <- cors[i,i]
    x <- cors[upper.tri(cors)]
    ePsi <- 2*k
    varPsi <- sum(3.263*x + 0.710*x^2 + 0.027*x^3)*2 + 4*k
    f <- ePsi*ePsi/varPsi
    c <- 0.5*varPsi/ePsi
    newP <- sum(-2*log(p[i])) / c
    pchisq(newP, df = 2*f, lower.tail = FALSE)
}
hmCors <- hmFry %>%
    dplyr::select(gs_name, fry = PValue) %>%
    left_join(
        dplyr::select(hmCamera, gs_name, camera = PValue)
    ) %>%
    left_join(
        dplyr::select(hmGsea, gs_name, gsea = PValue)
    ) %>%
    left_join(
        dplyr::select(hmRiboGoseq, gs_name, goseq = PValue)
    ) %>%
    as.data.frame() %>%
    column_to_rownames("gs_name") %>%
    as.matrix() %>%
    qnorm() %>%
    cor(use = "pairwise.complete.obs")
hmFry %>%
    dplyr::select(gs_name, fry = PValue) %>%
    left_join(
        dplyr::select(hmCamera, gs_name, camera = PValue)
    ) %>%
    left_join(
        dplyr::select(hmGsea, gs_name, gsea = PValue)
    ) %>%
    left_join(
        dplyr::select(hmRiboGoseq, gs_name, goseq = PValue)
    ) %>%
    as.data.frame() %>% 
    column_to_rownames("gs_name") %>%
    apply(1, kostP, cors = hmCors) %>%
    as.data.frame() %>% 
    set_names("p") %>% 
    rownames_to_column("gs_name") %>% 
    as_tibble() %>% 
    arrange(p) %>% 
    mutate(
        FDR = p.adjust(p, "fdr"),
        adjP = p.adjust(p, "bonf")
    ) %>%
    dplyr::filter(FDR < 0.05)

kgCors <- kgFry %>%
    dplyr::select(gs_name, fry = PValue) %>%
    left_join(
        dplyr::select(kgCamera, gs_name, camera = PValue)
    ) %>%
    left_join(
        dplyr::select(kgGsea, gs_name, gsea = PValue)
    ) %>%
    left_join(
        dplyr::select(kgRiboGoseq, gs_name, goseq = PValue)
    ) %>%
    as.data.frame() %>%
    column_to_rownames("gs_name") %>%
    as.matrix() %>%
    qnorm() %>%
    cor(use = "pairwise.complete.obs")
kgFry %>%
    dplyr::select(gs_name, fry = PValue) %>%
    left_join(
        dplyr::select(kgCamera, gs_name, camera = PValue)
    ) %>%
    left_join(
        dplyr::select(kgGsea, gs_name, gsea = PValue)
    ) %>%
    left_join(
        dplyr::select(kgRiboGoseq, gs_name, goseq = PValue)
    ) %>%
    as.data.frame() %>% 
    column_to_rownames("gs_name") %>%
    apply(1, kostP, cors = kgCors) %>%
    as.data.frame() %>% 
    set_names("p") %>% 
    rownames_to_column("gs_name") %>% 
    as_tibble() %>% 
    arrange(p) %>% 
    mutate(
        FDR = p.adjust(p, "fdr"),
        adjP = p.adjust(p, "bonf")
    ) %>%
    dplyr::filter(FDR < 0.05)

## Try using the harmonic mean
## This seems to give way too low values for some though
hmFry %>%
    dplyr::select(gs_name, fry = PValue) %>%
    left_join(
        dplyr::select(hmCamera, gs_name, camera = PValue)
    ) %>%
    left_join(
        dplyr::select(hmGsea, gs_name, gsea = PValue)
    ) %>%
    left_join(
        dplyr::select(hmRiboGoseq, gs_name, goseq = PValue)
    ) %>%
    as.data.frame() %>% 
    column_to_rownames("gs_name") %>% 
    as.matrix() %>% 
    apply(MARGIN = 1, FUN = hmp.stat) %>% 
    as.data.frame() %>% 
    set_names("p") %>% 
    rownames_to_column("gs_name") %>% 
    as_tibble() %>% 
    arrange(p) %>% 
    mutate(
        FDR = p.adjust(p, "fdr"),
        adjP = p.adjust(p, "bonf")
    ) %>%
    dplyr::filter(FDR < 0.05)
