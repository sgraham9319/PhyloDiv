

#===================================
# Evaluate habitat preference traits
#===================================

hab_pref_summ <- function(ag_dat, fen_dat, pas_dat){
  output <- ag_dat %>%
    select(species) %>%
    mutate(agriculture = ag_dat$Agriculture / ag_dat$total,
           ag_adj = (ag_dat$Agriculture + 0.5) / (ag_dat$total + 1),
           fenced = fen_dat$Fenced / fen_dat$total,
           fen_adj = (fen_dat$Fenced + 0.5) / (fen_dat$total + 1),
           pastoral = pas_dat$Pastoral / pas_dat$total,
           pas_adj = (pas_dat$Pastoral + 0.5) / (pas_dat$total + 1),
           con_ag = ag_dat$Conserved / ag_dat$total,
           con_fen = fen_dat$Conserved / fen_dat$total,
           con_pas = pas_dat$Conserved / pas_dat$total,
           con_ag_adj = (ag_dat$Conserved + 0.5) / (ag_dat$total + 1),
           con_fen_adj = (fen_dat$Conserved + 0.5) / (fen_dat$total + 1),
           con_pas_adj = (pas_dat$Conserved + 0.5) / (pas_dat$total + 1))
  output$conserved <- apply(output[, 8:10], 1, mean, na.rm = T)
  output$con_adj <- apply(output[, 11:13], 1, mean, na.rm = T)
  output <- output %>%
    mutate(agriculture_logit = log(ag_adj / (1 - ag_adj)),
           fenced_logit = log(fen_adj / (1 - fen_adj)),
           pastoral_logit = log(pas_adj / (1 - pas_adj)),
           conserved_logit = log(con_adj / (1 - con_adj))) %>%
    select(-ag_adj, -fen_adj, -pas_adj, -con_ag, -con_fen, -con_pas,
           -con_ag_adj, -con_fen_adj, -con_pas_adj, -con_adj)
  output[which(is.na(output$agriculture)), "agriculture_logit"] <- NaN
  output[which(is.na(output$fenced)), "fenced_logit"] <- NaN
  output[which(is.na(output$pastoral)), "pastoral_logit"] <- NaN
  output[which(is.na(output$conserved)), "conserved_logit"] <- NaN
  output
}

