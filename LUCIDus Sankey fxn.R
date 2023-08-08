# This function plot top 25% of the feature in each omics layer


reorder_lucid <- function(model,
                          order) {
  # record parameters 
  pars <- model$pars
  beta <- pars$beta
  mu <- pars$mu
  sigma <- pars$sigma
  gamma <- pars$gamma
  K <- model$K
  
  # 1 - reorder exposure effect
  # use the estimate from the 
  ref_cluster <- order[1]
  beta <- t(t(beta) - beta[ref_cluster, ])[order, ]
  
  # 2 - reorder omic effect
  mu <- mu[order, ]
  var <- vector(mode = "list", length = K)
  for (i in 1:K) {
    var[[i]] <- sigma[[order[i]]]
  }
  
  # 3 - reorder outcome effect
  gamma$beta[1:K] <- gamma$beta[order]
  gamma$sigma <- gamma$sigma[order]
  
  model$pars <- list(beta = beta,
                     mu = mu,
                     sigma = sigma,
                     gamma = gamma)
  
  return(model)
}

# Get sankey dataframe ----
get_sankey_df <- function(x,
                          G_color = "dimgray", 
                          X_color = "#eb8c30",
                          Z_color = "#2fa4da", 
                          Y_color = "#afa58e", 
                          pos_link_color = "#67928b", 
                          neg_link_color = "#d1e5eb", 
                          fontsize = 10) {
  K <- x$K
  var.names <- x$var.names
  pars <- x$pars
  dimG <- length(var.names$Gnames)
  dimZ <- length(var.names$Znames)
  valueGtoX <- as.vector(t(x$pars$beta[, -1]))
  valueXtoZ <- as.vector(t(x$pars$mu))
  valueXtoY <- as.vector(x$pars$gamma$beta)[1:K]
  
  # GtoX
  GtoX <- data.frame(
    source = rep(x$var.names$Gnames, K), 
    target = paste0("Latent Cluster", 
                    as.vector(sapply(1:K, function(x) rep(x, dimG)))), 
    value = abs(valueGtoX), 
    group = as.factor(valueGtoX > 0))
  
  # XtoZ
  XtoZ <- data.frame(
    source = paste0("Latent Cluster", 
                    as.vector(sapply(1:K, 
                                     function(x) rep(x, dimZ)))), 
    target = rep(var.names$Znames, 
                 K), value = abs(valueXtoZ),
    group = as.factor(valueXtoZ > 
                        0))
  
  # subset top 25% of each omics layer
  # top25<- XtoZ %>% 
  #   filter(source == "Latent Cluster2") %>%
  #   arrange(desc(value)) %>%
  #   head(n = 21) %>%
  #   select(target)
  # 
  # XtoZ_sub<- XtoZ %>% 
  #   filter(target %in% top25$target)
    
  
  # XtoY
  XtoY <- data.frame(source = paste0("Latent Cluster", 1:K), 
                     target = rep(var.names$Ynames, K), value = abs(valueXtoY), 
                     group = as.factor(valueXtoY > 0))
  # links <- rbind(GtoX, XtoZ_sub, XtoY)
  links <- rbind(GtoX, XtoZ, XtoY)
  
  nodes <- data.frame(
    name = unique(c(as.character(links$source), 
                    as.character(links$target))), 
    # group = as.factor(c(rep("exposure", 
    #                         dimG), rep("lc", K), rep("biomarker", dimZ * 0.25), "outcome")))
    group = as.factor(c(rep("exposure", 
                            dimG), rep("lc", K), rep("biomarker", dimZ), "outcome")))
  ## the following two lines were used to exclude covars from the plot
  links <- links %>% filter(!grepl("cohort", source))
  nodes <- nodes %>% filter(!grepl("cohort", name)) 
  
  links$IDsource <- match(links$source, nodes$name) - 1
  links$IDtarget <- match(links$target, nodes$name) - 1
  
  color_scale <- data.frame(
    domain = c("exposure", "lc", "biomarker", 
               "outcome", "TRUE", "FALSE"), 
    range = c(G_color, X_color, 
              Z_color, Y_color, pos_link_color, neg_link_color))
  
  sankey_df = list(links = links, 
                   nodes = nodes)
  return(sankey_df)
}




# Sankey Function ----

# lucid_fit1 <- fit1

sankey_diagram_plotly <- function(lucid_fit1, text_size = 15, exposure) {
  # 1. Get sankey dataframes ----
  sankey_dat <- get_sankey_df(lucid_fit1)
  n_omics <- length(lucid_fit1$var.names$Znames)
  # link data
  links <- sankey_dat[["links"]] 
  # node data
  nodes <- sankey_dat[["nodes"]] 
  
  nodes1 <- nodes %>% 
    mutate(
          # group = case_when(str_detect(name,"Cluster") ~ "lc",
          #                    str_detect(group, "cg") ~ "CpG",
          #                    str_detect(name, "TC") ~ "TC",
          #                    str_detect(name, "pro") ~ "Protein",
          #                    str_detect(name, "G1") ~ "exposure",
          #                    str_detect(name, "outcome") ~ "outcome"),
           name = ifelse(name == "G1", exposure,name))
  links1 <- links %>%
    mutate(source = ifelse(source == "G1", exposure,source))

  ## 6.1 Set Node Color Scheme: ----
  color_pal_sankey <- matrix(
    c("exposure", sankey_colors$range[sankey_colors$domain == "exposure"],
      "lc",       "#b3d8ff",
      "biomarker",  sankey_colors$range[sankey_colors$domain == "biomarker"],
      # "TC",      sankey_colors$range[sankey_colors$domain == "layer2"],
      # "Protein", sankey_colors$range[sankey_colors$domain == "layer3"],
      "outcome",  sankey_colors$range[sankey_colors$domain == "Outcome"]), 
    ncol = 2, byrow = TRUE) |>
    as_tibble(.name_repair = "unique") |> 
    janitor::clean_names() |>
    dplyr::rename(group = x1, color = x2)
  
  # Add color scheme to nodes
  nodes_new_plotly <- nodes1 |> 
    tidylog::left_join(color_pal_sankey) |>
    mutate(
      x = case_when(
        group == "exposure" ~ 0,
        str_detect(name, "Cluster") ~ 1/3,
        str_detect(group, "biomarker")|
          str_detect(name, "outcome")~ 2/3
      ))
  
  
  ## 6.2 Get links for Plotly, set color ----
  links_new <- links1  |>
    mutate(
      link_color = case_when(
        # Ref link color
        value == 0 ~     "#f3f6f4",
        # # Cluster 
        # str_detect(source, "Cluster1") &  group == TRUE  ~  "#706C6C",
        # str_detect(source, "Cluster1") &  group == FALSE ~  "#D3D3D3",
        # str_detect(source, "Cluster2") &  group == TRUE  ~  "#706C6C",
        # str_detect(source, "Cluster2") &  group == FALSE ~  "#D3D3D3",
        ##############
        # Exposure
        str_detect(source, exposure) &  group == TRUE  ~  "red",
            # Outcome
        str_detect(target, "outcome") &  group == TRUE  ~  "red",
        # Methylation 
        # str_detect(target, "_") &  group == TRUE  ~  "#bf9000",
        # str_detect(target, "_") &  group == FALSE ~  "#ffd966",
        # Transcriptome
        # str_detect(target, "TC") &  group == TRUE  ~  "#38761d",
        # str_detect(target, "TC") &  group == FALSE ~  "#b6d7a8",
        # # proteome
        str_detect(target, "_") &  group == TRUE  ~  "#a64d79",
        str_detect(target, "_") &  group == FALSE ~  "#ead1dc",
        ##
        group == FALSE ~ "#D3D3D3", # Negative association
        group == TRUE ~  "#706C6C")) # Positive association
  
  plotly_link <- list(
    source = links_new$IDsource,
    target = links_new$IDtarget,
    value = links_new$value+.00000000000000000000001, 
    color = links_new$link_color)  
  
  # Get list of nodes for Plotly
  plotly_node <- list(
    label = nodes_new_plotly$name, 
    color = nodes_new_plotly$color,
    pad = 15,
    thickness = 20,
    line = list(color = "black",width = 0.5),
    x = nodes_new_plotly$x, 
    # y = c(0.01, 
    #       0.3, 0.7, # clusters
    #       seq(from = .01, to = 1, by = 0.04)[1:(dimZ * 0.25)], # biomaker
    #       .95
    y = c(0.01,
          0.3, 0.75, # clusters
          seq(from = 0.1, to = 1, by = 0.08)[1:4],
          seq(from = 0.1+0.08*5, to = 1, by = 0.16)[1:2],
          0.98
  ))
  
  
  ## 6.3 Plot Figure ----
  (fig <- plot_ly(
    type = "sankey",
    domain = list(
      x =  c(0,1),
      y =  c(0,1)),
    orientation = "h",
    node = plotly_node,
    link = plotly_link))
  
  (fig <- fig %>% layout(
    # title = "Basic Sankey Diagram",
    font = list(
      size = text_size
    ))
  )
  
  return(fig)
}
