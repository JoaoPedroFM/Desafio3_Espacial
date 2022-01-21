## Aula 12
## Estatística Espacial I
## Rafael Erbisti


## EXEMPLO: DADOS CENSITÁRIOS DO AMAZONAS

library(rgdal)
library(maptools)
library(sp)
library(spdep)
library(spatialreg)

# Lendo o shapefile no rgdal
Amazonas <- readOGR(dsn="C:\\Users\\RAFAEL\\Dropbox\\UFF\\CURSOS\\2021\\Semestre 2021.2\\GET00153 - Estatística Espacial I\\Aulas\\Amazonas",layer="AM_Mun97_region")

# ou usando o maptools (alterando setwd para onde o arquivo shape está)
Amazonas <- readShapeSpatial("AM_Mun97_region")

plot(Amazonas)

# Lendo os dados (cuidar ordem e eliminar espaços entre palavras)
dados <- read.table("AM_Mun97_region.txt",header=T)

# Adicionando os dados ao objeto que contém o shapefile (não precisa, mas ajuda)
Amazonas@data = cbind(Amazonas@data,dados)

### encontrando quem é vizinho de quem
nb <- poly2nb(Amazonas,queen=TRUE)

### ponderacao (linhas somando 1)
nbw <- nb2listw(nb,style="W")
names(nbw)



### Regressão Espacial
### Variáveis da regressão: IDH versus Estabelecimentos de Saúde per capita

y <- Amazonas@data$Índice_de_Desenvolvimento_Humano_Municipal_._2010_.IDHM_2010.
x <- Amazonas@data$Estabelecimentos_de_Saúde_SUS/Amazonas@data$PIB_per_capita_a_preços_correntes

SEM <- errorsarlm(y~x, data=data.frame(cbind(x,y)), nbw, method="eigen")
SAR <- lagsarlm(y~x, data=data.frame(cbind(x,y)), nbw, method="eigen")
summary(SEM)
summary(SAR)

# Análise dos resíduos

Amazonas@data$residuos_SEM <- SEM$residuals
Amazonas@data$residuos_SAR <- SAR$residuals
Amazonas@data$residuos_DIF <- SEM$residuals-SAR$residuals

spplot(Amazonas,"residuos_SEM",col.regions=terrain.colors(256))
spplot(Amazonas,"residuos_SAR",col.regions=terrain.colors(256))
spplot(Amazonas,"residuos_DIF",col.regions=terrain.colors(256))

# cuidado com hipótese alternativa
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
