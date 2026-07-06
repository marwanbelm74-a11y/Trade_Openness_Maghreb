# Quelle ouverture compte pour la croissance au Maghreb ?
# Algérie, Maroc, Tunisie, 1990-2022
# Version 2 : période allongée, ouverture décomposée (commerce, IDE,
# exportations manufacturières), dates de rupture propres à chaque pays
# Toutes les données viennent de la Banque mondiale (package WDI)

library(WDI)
library(dplyr)
library(ggplot2)
library(tseries)
library(ARDL)
library(lmtest)
library(vars)
library(plm)
library(modelsummary)

set.seed(2024)
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

debut <- 1990
fin <- 2022
pays_noms <- c(DZ = "Algérie", MA = "Maroc", TN = "Tunisie")


# ---- données ----
# ouv   : commerce (exports + imports) en % du PIB
# ide   : entrées nettes d'IDE en % du PIB
# manuf : exportations manufacturières en % des exportations de marchandises
# fuel  : exportations de combustibles en % des exportations de marchandises

indic <- c(pib_hab = "NY.GDP.PCAP.KD",
           ouv     = "NE.TRD.GNFS.ZS",
           ide     = "BX.KLT.DINV.WD.GD.ZS",
           manuf   = "TX.VAL.MANF.ZS.UN",
           fuel    = "TX.VAL.FUEL.ZS.UN",
           inv     = "NE.GDI.FTOT.ZS",
           infl    = "FP.CPI.TOTL.ZG",
           pop     = "SP.POP.GROW",
           rente   = "NY.GDP.PETR.RT.ZS")

# petits trous au milieu d'une série : interpolation linéaire
interp <- function(x) {
  if (sum(!is.na(x)) < 2) return(x)
  zoo::na.approx(x, na.rm = FALSE)
}

dat <- WDI(country = names(pays_noms), indicator = indic,
           start = debut, end = fin) |>
  mutate(across(all_of(names(indic)), as.numeric),
         pays = unname(pays_noms[iso2c]),
         ly = log(pib_hab)) |>
  arrange(pays, year) |>
  group_by(pays) |>
  mutate(across(c(manuf, fuel), interp)) |>
  ungroup() |>
  as.data.frame()

# nombre d'années disponibles par variable : à regarder avant de continuer,
# manuf et fuel peuvent commencer après 1990
dat |>
  group_by(pays) |>
  summarise(across(c(ouv, ide, manuf, fuel, inv, rente), ~ sum(!is.na(.x))))

# dates de rupture
# entrée en vigueur de l'accord d'association avec l'UE
ue <- c("Tunisie" = 1998, "Maroc" = 2000, "Algérie" = 2005)
dat$post_ue <- as.numeric(dat$year >= unname(ue[dat$pays]))
# Printemps arabe
dat$post11 <- as.numeric(dat$year >= 2011)

datasummary(pib_hab + ouv + ide + manuf + fuel + inv + infl + rente ~ pays * (Mean + SD),
            data = dat, output = "tables/stats_desc.html")


# ---- graphiques ----

courbe <- function(var, ylab, fichier) {
  g <- ggplot(dat, aes(year, .data[[var]], colour = pays)) +
    geom_line(linewidth = 1) +
    labs(x = NULL, y = ylab, colour = NULL) +
    theme_minimal()
  ggsave(paste0("figures/", fichier, ".png"), g, width = 7, height = 4, dpi = 300)
  g
}

g1 <- courbe("ouv", "Commerce extérieur (% du PIB)", "ouverture")
g2 <- courbe("pib_hab", "PIB par habitant (USD constants 2015)", "pib_hab")
g3 <- courbe("manuf", "Exportations manufacturières (% des exportations)", "manuf")
g4 <- courbe("fuel", "Exportations de combustibles (% des exportations)", "fuel")
g5 <- courbe("ide", "Entrées nettes d'IDE (% du PIB)", "ide")


# ---- racines unitaires ----

racine <- function(x) {
  x <- na.omit(x)
  tryCatch(
    suppressWarnings(c(adf_niveau = adf.test(x)$p.value,
                       adf_diff = adf.test(diff(x))$p.value,
                       kpss_niveau = kpss.test(x)$p.value)),
    error = function(e) c(adf_niveau = NA, adf_diff = NA, kpss_niveau = NA))
}

vars_test <- c("ly", "ouv", "ide", "manuf", "inv", "rente")

ur <- lapply(split(dat, dat$pays), function(d) round(t(sapply(d[vars_test], racine)), 3))
ur   # ADF rejette la racine unitaire si p < 0.05, KPSS rejette la stationnarité si p < 0.05

sink("tables/racines_unitaires.txt")
print(ur)
sink()


# ---- ARDL pays par pays ----
# spec A : commerce + IDE (+ rentes pour l'Algérie)
# spec B : composition des exportations (manufacturier) + IDE

specs <- list(
  A = list("Algérie" = ly ~ ouv + ide + inv + rente,
           "Maroc"   = ly ~ ouv + ide + inv,
           "Tunisie" = ly ~ ouv + ide + inv),
  B = list("Algérie" = ly ~ manuf + ide + inv,
           "Maroc"   = ly ~ manuf + ide + inv,
           "Tunisie" = ly ~ manuf + ide + inv)
)

modeles <- list()
lt <- data.frame()
bornes <- data.frame()

for (s in names(specs)) {
  for (p in names(specs[[s]])) {
    f <- specs[[s]][[p]]
    v <- all.vars(f)
    
    # on garde les années où toutes les variables du modèle existent
    d <- dat[dat$pays == p, ]
    d <- d[complete.cases(d[, v]), ]
    if (any(diff(d$year) != 1)) stop("trou dans la série : ", s, " ", p)
    x <- ts(d[, v], start = min(d$year))
    
    sel <- auto_ardl(f, data = x, max_order = rep(2, length(v)), selection = "AIC")
    m <- sel$best_model
    modeles[[paste(s, p)]] <- m
    
    cat("\n-----", s, "-", p, "(", min(d$year), "-", max(d$year), ") -----\n")
    cat("ordre retenu :", sel$best_order, "\n")
    
    bf <- bounds_f_test(m, case = 3)
    print(bf)
    
    mult <- multipliers(m)
    print(mult)
    
    print(summary(recm(m, case = 3)))
    
    bg <- bgtest(m, order = 2)
    print(bg)
    print(bptest(m))
    print(shapiro.test(residuals(m)))
    
    bornes <- rbind(bornes, data.frame(
      spec = s, pays = p,
      periode = paste(range(d$year), collapse = "-"),
      ordre = paste(sel$best_order, collapse = ","),
      F_stat = round(as.numeric(bf$statistic), 3),
      p_bornes = round(bf$p.value, 3),
      p_autocorr = round(bg$p.value, 3)))
    
    mult$spec <- s
    mult$pays <- p
    lt <- rbind(lt, mult)
  }
}

bornes
write.csv(bornes, "tables/bornes.csv", row.names = FALSE)
write.csv(lt, "tables/long_terme.csv", row.names = FALSE)

# les p-values du package sont asymptotiques ; avec T autour de 30,
# comparer aussi la stat F aux valeurs de Narayan (2005)


# ---- causalité de Toda-Yamamoto ----
# d = ordre d'intégration max supposé (1)

ty <- function(y, x, k, d = 1) {
  y <- as.numeric(y)
  x <- as.numeric(x)
  p <- k + d
  Y <- embed(y, p + 1)
  X <- embed(x, p + 1)
  dat_ty <- data.frame(dep = Y[, 1], Y[, -1], X[, -1])
  names(dat_ty) <- c("dep", paste0("y", 1:p), paste0("x", 1:p))
  complet <- lm(dep ~ ., data = dat_ty)
  restreint <- lm(dep ~ ., data = dat_ty[, !(names(dat_ty) %in% paste0("x", 1:k))])
  waldtest(restreint, complet, test = "Chisq")$`Pr(>Chisq)`[2]
}

res_ty <- data.frame()
for (p in names(pays_noms) |> unname()) {
  for (var in c("ouv", "ide", "manuf")) {
    d <- dat[dat$pays == p, ]
    d <- d[complete.cases(d[, c("ly", var)]), ]
    k <- unname(VARselect(as.matrix(d[, c("ly", var)]), lag.max = 3,
                          type = "const")$selection["AIC(n)"])
    res_ty <- rbind(res_ty, data.frame(
      pays = p, variable = var, retards = k,
      p_var_vers_pib = round(ty(d$ly, d[[var]], k), 3),
      p_pib_vers_var = round(ty(d[[var]], d$ly, k), 3)))
  }
}
res_ty
write.csv(res_ty, "tables/toda_yamamoto.csv", row.names = FALSE)


# ---- panel ----
# variable dépendante : croissance du PIB par habitant
# attention : plm masque lag() de dplyr, d'où le dplyr:: explicite

dat <- dat |>
  group_by(pays) |>
  mutate(dly = ly - dplyr::lag(ly)) |>
  ungroup() |>
  as.data.frame()

summary(dat$dly)   # doit être non nul, avec 3 NA (une année de départ par pays)

pd <- pdata.frame(dat, index = c("pays", "year"))

p1 <- plm(dly ~ ouv + ide + inv + infl + pop, data = pd, model = "within")
p2 <- plm(dly ~ manuf + ide + inv + infl + pop, data = pd, model = "within")
p3 <- plm(dly ~ ouv + ide + inv + infl + pop + post_ue + ouv:post_ue + ide:post_ue,
          data = pd, model = "within")
p4 <- plm(dly ~ ouv + ide + inv + infl + pop + post11 + ouv:post11 + ide:post11,
          data = pd, model = "within")
p5 <- plm(dly ~ ouv + rente + ouv:rente + ide + inv + infl + pop,
          data = pd, model = "within")

panels <- list("1 Baseline" = p1,
               "2 Manufacturier" = p2,
               "3 Accord UE" = p3,
               "4 Printemps arabe" = p4,
               "5 Rentes" = p5)

tab <- function(t, nom) {
  data.frame(modele = nom, variable = rownames(t), round(t[, 1:4], 4),
             row.names = NULL)
}

# erreurs-types Driscoll-Kraay
res_panel <- do.call(rbind, lapply(names(panels), function(n) {
  m <- panels[[n]]
  t <- coeftest(m, vcov = vcovSCC(m, type = "HC1"))
  cat("\n", n, "\n")
  print(t)
  tab(t, n)
}))

write.csv(res_panel, "tables/panel.csv", row.names = FALSE)