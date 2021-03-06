## Aula 12
## Estat�stica Espacial I
## Rafael Erbisti


## EXEMPLO: DADOS CENSIT�RIOS DO AMAZONAS

library(rgdal)
library(maptools)
library(sp)
library(spdep)
library(spatialreg)

# Lendo o shapefile no rgdal
Amazonas <- readOGR(dsn="C:\\Users\\RAFAEL\\Dropbox\\UFF\\CURSOS\\2021\\Semestre 2021.2\\GET00153 - Estat�stica Espacial I\\Aulas\\Amazonas",layer="AM_Mun97_region")

# ou usando o maptools (alterando setwd para onde o arquivo shape est�)
Amazonas <- readShapeSpatial("AM_Mun97_region")

plot(Amazonas)

# Lendo os dados (cuidar ordem e eliminar espa�os entre palavras)
dados <- read.table("AM_Mun97_region.txt",header=T)

# Adicionando os dados ao objeto que cont�m o shapefile (n�o precisa, mas ajuda)
Amazonas@data = cbind(Amazonas@data,dados)

### encontrando quem � vizinho de quem
nb <- poly2nb(Amazonas,queen=TRUE)

### ponderacao (linhas somando 1)
nbw <- nb2listw(nb,style="W")
names(nbw)



### Regress�o Espacial
### Vari�veis da regress�o: IDH versus Estabelecimentos de Sa�de per capita

y <- Amazonas@data$�ndice_de_Desenvolvimento_Humano_Municipal_._2010_.IDHM_2010.
x <- Amazonas@data$Estabelecimentos_de_Sa�de_SUS/Amazonas@data$PIB_per_capita_a_pre�os_correntes

SEM <- errorsarlm(y~x, data=data.frame(cbind(x,y)), nbw, method="eigen")
SAR <- lagsarlm(y~x, data=data.frame(cbind(x,y)), nbw, method="eigen")
summary(SEM)
summary(SAR)

# An�lise dos res�duos

Amazonas@data$residuos_SEM <- SEM$residuals
Amazonas@data$residuos_SAR <- SAR$residuals
Amazonas@data$residuos_DIF <- SEM$residuals-SAR$residuals

spplot(Amazonas,"residuos_SEM",col.regions=terrain.colors(256))
spplot(Amazonas,"residuos_SAR",col.regions=terrain.colors(256))
spplot(Amazonas,"residuos_DIF",col.regions=terrain.colors(256))

# cuidado com hip�tese alternativa
imoran.SEM <- moran.mc(SEM$residuals, nbw, nsim=999, alternative="less")
imoran.SAR <- moran.mc(SAR$residuals, nbw, nsim=999, alternative="less")
imoran.SEM
imoran.SAR

----------------------------------------------
  
# OBS:
#   
# Spatial simultaneous autoregressive error model estimation (SEM):
# y = X beta + u
# u = lambda W u + e
# 
# Spatial simultaneous autoregressive lag model estimation (SAR)
# y = rho W y + X beta + e
