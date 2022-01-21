### Aula 13 - Estatística Espacial I
### Rafael Erbisti


## Exemplo 1

check_sigma = function(rho) {
  Sigma_inv = matrix(c(1,-rho,0,-rho,2,-rho,0,-rho,1), 
                     ncol=3, byrow=TRUE)
  solve(Sigma_inv)
}

check_sigma(rho=0)
check_sigma(rho=0.5)
check_sigma(rho=-0.6)
check_sigma(rho=1)
check_sigma(rho=-1)
check_sigma(rho=1.5)
check_sigma(rho=-1.5)



## Exemplo 2: dados da síndrome da morte súbita infantil (SIDS) 
## na Carolina do Norte (EUA)

library(dplyr)
library(sf)
library(ggplot2)
library(spdep)
library(spatialreg)


# Carregando o dado e o shape (esta tudo no mesmo objeto)
sids <- st_read(dsn = system.file("shapes/sids.shp", package = "spData"), crs = 4326 )

# Visualizando a base
head(sids)

# Mapa dos nascidos vivos
ggplot(sids)+
  geom_sf(aes(fill = BIR79))+
  scale_fill_viridis_c()

# Histogram dos nascimentos
ggplot(sids)+
  geom_histogram(aes(x = BIR79))


# Mapa dos óbitos por síndrome da morte súbita infantil
ggplot(sids)+
  geom_sf(aes(fill = SID79))+
  scale_fill_viridis_c()

# Histogram dos óbitos por síndrome da morte súbita infantil
ggplot(sids)+
  geom_histogram(aes(x = SID79))


# Criando a variável de taxa de morte a cada mil nascimentos de 1979 a 1984
sids['sids_rate79'] <- (sids['SID79'] * 1000)/ sids['BIR79']

# Mapa da taxa de SIDS por 1000 nascimentos
ggplot(sids)+
  geom_sf(aes(fill = sids_rate79))+
  scale_fill_viridis_c()

# Histogram da taxa de SIDS por 1000 nascimentos
ggplot(sids)+
  geom_histogram(aes(x = sids_rate79))


## MODELAGEM

# Calculando a matriz de vizinhança
W = st_touches(sids, sparse=FALSE)
listW = mat2listw(W)


# Calculando os centroides e as conexões (ja exitem na base lon e lat)

sids_coords = sids %>% 
              st_centroid() %>% 
              st_coordinates()

plot(st_geometry(sids))
plot(listW, sids_coords, add=TRUE, col=4, pch=16)


## MODELO CAR


# Criando a variável de taxa de morte a cada mil nascimentos de 1974 a 1978
# para usar como covariavel
sids['sids_rate74'] <- (sids['SID74'] * 1000)/ sids['BIR74']

x11()
ggplot(sids)+
  geom_sf(aes(fill = sids_rate79))+
  scale_fill_viridis_c()

nc_car <- spautolm(formula = sids_rate79 ~ sids_rate74,
                   data = sids,
                   listw = listW, family = "CAR")

summary(nc_car)


## MODELO SAR

nc_sar <- spautolm(formula = sids_rate79 ~ sids_rate74, data = sids,
                   listw = listW, family = "SAR")

summary(nc_sar)


nc_ind <- lm(sids_rate79 ~ sids_rate74, data = sids)


summary(nc_ind)
AIC(nc_ind)


## Avaliando os modelos

# residuos

sids['r_car'] <- nc_car$fit$residuals
sids['r_sar'] <- nc_sar$fit$residuals
sids['r_ind'] <- nc_ind$residuals


ggplot(sids)+
  geom_sf(aes(fill = r_car))+
  scale_fill_viridis_c()

ggplot(sids)+
  geom_sf(aes(fill = r_sar))+
  scale_fill_viridis_c()


ggplot(sids)+
  geom_sf(aes(fill = r_ind))+
  scale_fill_viridis_c()


# valores preditos

sids['fit_car'] <- nc_car$fit$fitted.values
sids['fit_sar'] <- nc_sar$fit$fitted.values


ggplot(sids)+
  geom_sf(aes(fill = fit_car))+
  scale_fill_viridis_c()

ggplot(sids)+
  geom_sf(aes(fill = fit_sar))+
  scale_fill_viridis_c()


# qqplot

ggplot(sids,aes(sample=r_car)) + 
  stat_qq() +
  stat_qq_line(col=4, lwd=1.1) +
  labs(x='quantis teóricos',y='quantis amostrais',title="CAR") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))


ggplot(sids,aes(sample=r_sar)) + 
  stat_qq() +
  stat_qq_line(col=4, lwd=1.1) +
  labs(x='quantis teóricos',y='quantis amostrais',title="SAR") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))

# residuos CAR vs SAR
ggplot(sids,aes(x = r_car , y = r_sar)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col=2, lwd=1.1) +
  labs(x='resíduos CAR',y='resíduos SAR', title= "CAR vs SAR - resíduos") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))


# residuos CAR vs IND
ggplot(sids,aes(x = r_car , y = r_ind)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col=2, lwd=1.1) +
  labs(x='resíduos CAR',y='resíduos IND', title= "CAR vs IND - resíduos") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))

