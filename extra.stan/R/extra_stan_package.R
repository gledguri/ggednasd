#' Plot ridges for a latent variable with no index
#'
#' This function plots a density plot for a latent variable in stan that has no index.
#' For example a variable in stan written "tau ~ normal(0,10)"
#' @param latent_var The number of colours to be produced
#' @return density plot (using ggplot2)
#' @export
#' @examples
#' stan.ridge.plots.f0("tau")
#' 
stan.ridge.plots.f0 <- function(latent_var){
  param <- as.data.frame(extract(stanMod, par = latent_var))
  tparam <- trans(param)
  g <- ggplot(tparam,aes(x = `Z`, y = `Y`))+
    geom_density_ridges(fill="#4CA6FF",scale = 2,size = 0.2,rel_min_height = 0.01,alpha = .8)+
    labs(title = paste("Posterior probability distribution of",latent_var))+
    ylab("Density")+
    theme_minimal()+
    xlab(latent_var)
  return(g)
}

#' Collapse the latent variable to 1 index
#'
#' This function merges the latent variable from stan output onto its given index and prepares
#' the data for ridged plots. Example of a latent variable written in stan is "mu[sample_idx]".
#' @param latent_var The name of the variable (i.e. "mu")
#' @param f1_colname The column name where the factor is as a vector (i.e. "sample_names")
#' @param idx_f1_colname The column name where the index is as a vector (i.e. "sample_idx")
#' @param input_data The dataframe where all the colums above exists
#' @param i Logical vector for stating the indexing of the parameters needs adjusting
#' @return arranged latent parameter output from stan
#' @export
#' @examples
#' stan_param_collapse_1f("sigma","species_name","species_idx",qpcr_df,i=F)
#' stan.plot.data <- stan_param_collapse_1f("phi_1","species_name","species_idx",st_qpcr)

stan_param_collapse_1f <- function (latent_var, f1_colname, idx_f1_colname, input_data,i=F)
{
  fa1 <- input_data[,f1_colname]
  faidx <- input_data[,idx_f1_colname]
  if(i==F)fa1 <- collapse_similar_rows(fa1, faidx)$factor
  param <- as.data.frame(extract(stanMod, par = latent_var))
  colnames(param) <- paste0(fa1)
  if(i==T)param <- collapse_similar_columns(param)
  tparam <- trans(param)
  tparam <- as.data.frame(tparam[,2:3])
  # tparam <- (tparam[!is.na(tparam$Z), ])
  tparam$Y <- factor((tparam$Y), level = rev(unique(tparam$Y)))
  tt1 <<- f1_colname
  parameter <<- latent_var
  return(tparam)
}
