library(WDI)
library(tidyverse)
library(fixest)
library(modelsummary)
install.packages("wbstats")
library(wbstats)

# 2. Liste des pays africains (codes ISO-2)
afrique <- c("DZ","AO","BJ","BW","BF","BI","CM","CV","CF","TD","KM","CG","CD",
             "CI","DJ","EG","GQ","ER","SZ","ET","GA","GM","GH","GN","GW","KE",
             "LS","LR","LY","MG","MW","ML","MR","MU","MA","MZ","NA","NE","NG",
             "RW","ST","SN","SC","SL","SO","ZA","SS","SD","TZ","TG","TN","UG",
             "ZM","ZW")

# 3. Téléchargement WDI (version wbstats)
df <- wb_data(
  indicator = c(
    acces_rural   = "EG.ELC.ACCS.RU.ZS",
    acces_urbain  = "EG.ELC.ACCS.UR.ZS",
    pib_pc        = "NY.GDP.PCAP.KD",
    urbanisation  = "SP.URB.TOTL.IN.ZS"
  ),
  country = afrique,
  start_date = 2000,
  end_date = 2022
)

# 4. Nettoyage et création de la variable d'écart
df_clean <- df %>%
  mutate(ecart = acces_urbain - acces_rural) %>%
  filter(!is.na(ecart)) %>%
  filter(!is.na(pib_pc), !is.na(urbanisation))

head(df_clean)

# 5. Estimation avec effets fixes pays et année
modele <- feols(
  ecart ~ acces_rural | iso2c + date, 
  data = df_clean,
  cluster = ~ iso2c
)

summary(modele)

# 6. Création du tableau de résultats propre
modelsummary(
  list("Modèle 1 (Effets Fixes)" = modele),
  stars = TRUE,
  gof_omit = "IC|Log|Adj|RMSE",
  title = "Impact de l'accès rural sur l'écart d'accès à l'électricité"
)

# 7. Modèle avec contrôles (Robustesse)
modele2 <- feols(
  ecart ~ acces_rural + pib_pc + urbanisation | iso2c + date, 
  data = df_clean,
  cluster = ~ iso2c
)

summary(modele2)

# 8. Création du tableau de résultats propre modele 2
modelsummary(
  list("Modèle 2 (Effets Fixes)" = modele2),
  stars = TRUE,
  gof_omit = "IC|Log|Adj|RMSE",
  title = "Impact de l'accès rural sur l'écart d'accès à l'électricité"
)

# 9. Visualisation 1 : Evolution de l'écart pour quelques pays
df_clean %>%
  filter(iso2c %in% c("RW", "KE", "NG", "ET", "ZA")) %>%
  ggplot(aes(x = date, y = ecart, color = country)) +
  geom_line(size = 1) +
  labs(
    title = "Évolution de l'écart d'accès à l'électricité",
    subtitle = "Pays sélectionnés (2000-2022)",
    x = "Année", 
    y = "Écart urbain - rural (points de %)",
    color = "Pays"
  ) +
  theme_minimal()

# 10. Visualisation 2 : Relation entre accès rural et écart
df_clean %>%
  ggplot(aes(x = acces_rural, y = ecart)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", color = "darkred") +
  labs(
    title = "Accès rural et écart d'accès",
    x = "Accès à l'électricité en zone rurale (%)",
    y = "Écart urbain - rural (points de %)"
  ) +
  theme_minimal()

ggsave("graphique_1_evolution.png", width = 8, height = 5, dpi = 300)
ggsave("graphique_2_nuage.png", width = 8, height = 5, dpi = 300)

df_clean %>%
  filter(iso2c %in% c("RW", "KE", "NG", "ET", "ZA")) %>%
  ggplot(aes(x = date, y = ecart, color = country)) +
  geom_line(size = 1) +
  labs(
    title = "Évolution de l'écart d'accès à l'électricité",
    subtitle = "Pays sélectionnés (2000-2022)",
    x = "Année", 
    y = "Écart urbain - rural (points de %)",
    color = "Pays"
  ) +
  theme_minimal()