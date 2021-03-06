### Aula 13 - Estat�stica Espacial I
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
check_sigma(rho=1)    # correla��o muito forte
check_sigma(rho=-1)   # correla��o muito forte
check_sigma(rho=1.5)  # correla��o n�o vi�vel
check_sigma(rho=-1.5) # correla��o n�o vi�vel



## Exemplo 2: dados da s�ndrome da morte s�bita infantil (SIDS) 
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
# ponta esquerda e bordas podem ter alguma rela��o

# Histograma dos nascimentos
ggplot(sids)+
  geom_histogram(aes(x = BIR79))
# assimetria na distribui��o

# Mapa dos �bitos por s�ndrome da morte s�bita infantil
ggplot(sids)+
  geom_sf(aes(fill = SID79))+
  scale_fill_viridis_c()
# mais �bitos nos lugares com mais nascimentos
# bordas com comportamento parecido com os nascimentos, nasce e morre menos

# Histograma dos �bitos por s�ndrome da morte s�bita infantil
ggplot(sids)+
  geom_histogram(aes(x = SID79))
# assimetria na distribui��o

# Modelo Poisson pode ser usado para modelar �bitos: contagem

# Modelar taxa de �bitos pelo modelo Normal
# Criando a vari�vel de taxa de morte a cada mil nascimentos de 1979 a 1984
sids['sids_rate79'] <- (sids['SID79'] * 1000)/ sids['BIR79']

# Mapa da taxa de SIDS por 1000 nascimentos
ggplot(sids)+
  geom_sf(aes(fill = sids_rate79))+
  scale_fill_viridis_c()

# Histograma da taxa de SIDS por 1000 nascimentos
ggplot(sids)+
  geom_histogram(aes(x = sids_rate79))
# Se aproxima mais de uma distribui��o cont�nua normal

## MODELAGEM

# Calculando a matriz de vizinhan�a
W = st_touches(sids, sparse=FALSE)
listW = mat2listw(W)
listW
# nonzero links: 490 - quantidade de vizinhos

# Calculando os centroides e as conex�es (ja existem na base lon e lat)

sids_coords = sids %>% 
              st_centroid() %>% 
              st_coordinates()

plot(st_geometry(sids))
plot(listW, sids_coords, add=TRUE, col=4, pch=16)


## MODELO CAR
# Regress�o n�o espacial: incorporar latitude e longitude como covari�veis 
# (informa��es do espa�o) e taxa do per�odo anterior ao estudado

# Criando a vari�vel de taxa de morte a cada mil nascimentos de 1974 a 1978
# para usar como covariavel
sids['sids_rate74'] <- (sids['SID74'] * 1000)/ sids['BIR74']

# x11()
ggplot(sids)+
  geom_sf(aes(fill = sids_rate74))+
  scale_fill_viridis_c()

## MODELO CAR com lon e lat
ggplot(sids)+
  geom_sf(aes(fill = sids_rate79))+
  scale_fill_viridis_c()

nc_car <- spautolm(formula = sids_rate79 ~ sids_rate74+ lon + lat,
                   data = sids,
                   listw = listW, family = "CAR")
summary(nc_car)
# lambda (rho) estimado pr�ximo de 0, dependencia espacial proxima de 0, 
# quase independente


## MODELO CAR sem lon e lat
nc_car <- spautolm(formula = sids_rate79 ~ sids_rate74,
                   data = sids,
                   listw = listW, family = "CAR")
summary(nc_car)
# lambda (rho) estimado pr�ximo de 0, dependencia espacial proxima de 0, 
# quase independente


## MODELO SAR
nc_sar <- spautolm(formula = sids_rate79 ~ sids_rate74, data = sids,
                   listw = listW, family = "SAR")

summary(nc_sar)
# Resultados muito parecidos com as do CAR, lambda menor, indicando modelo mais
# independente que o CAR


## MODELO LM
nc_ind <- lm(sids_rate79 ~ sids_rate74, data = sids)

summary(nc_ind)
AIC(nc_ind)
# dado n�o possui depend�ncia espacial
# rho= par�metro que estima o grau de depend�ncia espacial � pr�ximo de 0



## Avaliando os modelos
# y chapeu e residuos muito parecidos entre CAR, SAR e modelo independente 


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
# sem diferen�a nos qqplots, assimetria nos pontos extremos
ggplot(sids,aes(sample=r_car)) + 
  stat_qq() +
  stat_qq_line(col=4, lwd=1.1) +
  labs(x='quantis te�ricos',y='quantis amostrais',title="CAR") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))


ggplot(sids,aes(sample=r_sar)) + 
  stat_qq() +
  stat_qq_line(col=4, lwd=1.1) +
  labs(x='quantis te�ricos',y='quantis amostrais',title="SAR") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))

# residuos CAR vs SAR
ggplot(sids,aes(x = r_car , y = r_sar)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col=2, lwd=1.1) +
  labs(x='res�duos CAR',y='res�duos SAR', title= "CAR vs SAR - res�duos") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# residuos iguais, em cima da reta

# residuos CAR vs IND
ggplot(sids,aes(x = r_car , y = r_ind)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col=2, lwd=1.1) +
  labs(x='res�duos CAR',y='res�duos IND', title= "CAR vs IND - res�duos") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# residuos iguais, em cima da reta

