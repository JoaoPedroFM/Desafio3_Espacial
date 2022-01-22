### Desafio 3 - Estatística Espacial I

## Exemplo: Dados obtidos a partir do banco de dados Scottish Statistics
## (http://statistics.gov.scot), mas também estão incluídos no pacote 
## CARBayesdata no R. Dados da região de Greater Glasgow & Clyde. 

{library(CARBayesdata)
  library(sp)
  library(dplyr)
  library(sf)
  library(ggplot2)
  library(spdep)
  library(spatialreg)
  library(rgdal)
  library(leaflet)
  library(tidyverse)
  library(readxl)}

# Carregando os dados 
RJ= st_read("RJ_Municipios_2020.shp")
analfabetismo = read_excel("analfabetismo_rj.xlsx")
desemprego = read_excel("desemprego_rj.xlsx")
trab_infantil = read_excel("trab_infantil_rj.xlsx")
nascimentos = read_excel("nascimentos_rj.xlsx") 
  
base = nascimentos |> 
  left_join(analfabetismo) |> 
  left_join(desemprego) |> 
  left_join(trab_infantil) 

## Analise dos dados

# Histogram do nascimentos dos imoveis
ggplot(base)+
  geom_histogram(aes(x = nascimentos))

# Histogram do log do nascimentos dos imoveis
base <- base |> mutate(lognascimentos = log(base$nascimentos))
ggplot(base)+
  geom_histogram(aes(x = lognascimentos))
# variavel real continua, não possui assimetria

# log preco vs taxa analfabetismo
ggplot(base,aes(x = tx_analfabetismo , y = lognascimentos)) + 
  geom_point() +
  labs(x='taxa analfabetismo',y='log nascimentos') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# parece haver relação negativa entre as variáveis

# log preco vs taxa desemprego
ggplot(base,aes(x = tx_desemprego , y = lognascimentos)) + 
  geom_point() +
  labs(x='taxa desemprego',y='log nascimentos') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# grafico adequado

# log preco vs taxa trabalho infantil
ggplot(base,aes(x = tx_trab_infant , y = lognascimentos)) + 
  geom_point() +
  labs(x='taxa trabalho infantil',y='log nascimentos') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# não é possível concluir que há algum tipo de relação



## Análise espacial

# juntando as informacoes numa base unica
base.sp <- merge(base, RJ, by=c("CD_MUN","NM_MUN"), all.x=F, all.y=F)

# O pacote leaflet precisa das coordenadas em lat/long
# Precisamos transformar o sistema de coordenadas
base.sp <- spTransform(base.sp,CRS("+proj=longlat +datum=WGS84 +no_defs"))

# fazendo o mapa do nascimentos
colours <- colorNumeric(palette = "YlOrRd", domain = base.sp@data$nascimentos)
map1 <- leaflet(data=base.sp) |>
  addTiles() |>
  addPolygons(fillColor = ~colours(nascimentos), color="", weight=1,
              fillOpacity = 0.8) |>
  addLegend(pal = colours, values = base.sp@data$nascimentos, opacity = 1,
            title="nascimentos") |>
  addScaleBar(position="bottomleft")
map1
# mais escuros, maiores os nascimentoss
# localmente há áreas que formam clusters em relação aos nascimentoss


## Regressao Nao Espacial

form <- lognascimentos ~ tx_analfabetismo + tx_desemprego + tx_trab_infant 
model <- lm(formula=form, data=base.sp@data)
summary(model) # todas as variáveis significativas
AIC(model) 

# Matriz de vizinhança (W)
W.nb <- poly2nb(base.sp, row.names = rownames(base.sp@data))
W.list <- nb2listw(W.nb, style="B")
W.list

# Indice de Moran Global
globalMoran <- moran.test(residuals(model), W.list)
globalMoran
# H0: I=0 versus H1: I!=0
# Estatistica de teste 0.3 e p-valor pequeno -> rejeita H0, indicando que
# há estrutura de dependencia espacial nos residuos do modelo

# Indice de Moral Local (LISA)
localMoran <- localmoran(residuals(model), W.list)
localMoran

base.sp@data <- cbind(base.sp@data,LISA = localMoran[,4])


ggplot(base.sp@data)+
  geom_histogram(aes(x = LISA))
# as áreas em torno do 0 indicam não dependência
# as áreas com LISA maior ou menor indicam dependência


colours <- colorNumeric(palette = "Greens", domain = base.sp@data$LISA)
map2 <- leaflet(data=base.sp) |>
  addTiles() |>
  addPolygons(fillColor = ~colours(LISA), color="", weight=1,
              fillOpacity = 0.8) |>
  addLegend(pal = colours, values = base.sp@data$LISA, opacity = 1,
            title="nascimentos") |>
  addScaleBar(position="bottomleft")
map2
# valores nos grupos indicam indice alto ou baixo -> clusters, dependencia local
# dependencia global não tão alta

## Regressao Espacial

## MODELO CAR

form <- lognascimentos ~ tx_analfabetismo + tx_desemprego + tx_trab_infant 
nc_car <- spautolm(formula = form,
                   data = base.sp@data,
                   listw =  W.list,
                   family = "CAR")

summary(nc_car)
# todas as variaveis são significativas, exceto driveshop
# rho= 0.14788 (não 0 nem 1) dependencia mais local que global
# p-valor proximo de 0, o modelo não é independente
# AIC penaliza quanto ao número de parâmetros, mais completo que o modelo independente
# AIC menor que o do modelo independente-> melhor ajuste

## MODELO SAR

nc_sar <- spautolm(formula = form,
                   data = base.sp@data,
                   listw =  W.list,
                   family = "SAR")

summary(nc_sar)
# se aproxima do modelo CAR: variaveis são significativas, exceto driveshop
# rho= 0.10615 (não 0 nem 1) dependencia mais local que global
# p-valor prÓximo de 0, o modelo não é independente
# AIC menor que o do modelo independente, mas maior que o do modelo CAR
# Ajuste melhor para o modelo CAR

## MODELO INDEPENDENTE

summary(model)
AIC(model)



## Avaliando os modelos

# residuos

base.sp@data['r_car'] <- nc_car$fit$residuals
base.sp@data['r_sar'] <- nc_sar$fit$residuals
base.sp@data['r_ind'] <- model$residuals

colours <- colorNumeric(palette = "Blues", domain = c(-1,1))
map3 <- leaflet(data=base.sp) |>
  addTiles() |>
  addPolygons(fillColor = ~colours(r_car), color="", weight=1,
              fillOpacity = 0.8) |>
  addLegend(pal = colours, values = base.sp@data$r_car, opacity = 1,
            title="CAR") |>
  addScaleBar(position="bottomleft")
map3

map4 <- leaflet(data=base.sp) |>
  addTiles() |>
  addPolygons(fillColor = ~colours(r_sar), color="", weight=1,
              fillOpacity = 0.8) |>
  addLegend(pal = colours, values = base.sp@data$r_sar, opacity = 1,
            title="SAR") |>
  addScaleBar(position="bottomleft")
map4

map5 <- leaflet(data=base.sp) |>
  addTiles() |>
  addPolygons(fillColor = ~colours(r_ind), color="", weight=1,
              fillOpacity = 0.8) |>
  addLegend(pal = colours, values = base.sp@data$r_ind, opacity = 1,
            title="IND") |>
  addScaleBar(position="bottomleft")
map5

# Residuos do modelo CAR com magnitude menor do que dos outros modelos porque se
# ajustou melhor: Y ajustado próximo do Y verdadeiro
# Mapas parecidos, estruturas de dependencia presentes nas regioes mais escuras,
# efeitos maiores, blocos locais visto nas outras análises

# residuos CAR vs SAR
ggplot(base.sp@data,aes(x = r_car , y = r_sar)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col=2, lwd=1.1) +
  labs(x='resíduos CAR',y='resíduos SAR', title= "CAR vs SAR - resíduos") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# mais pontos acima da reta, residuos do SAR maiores que o do CAR
# Ajuste do Y pelo SAR foi mais distante do verdadeiro que o do CAR

# residuos CAR vs IND
ggplot(base.sp@data,aes(x = r_car , y = r_ind)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col=2, lwd=1.1) +
  labs(x='resíduos CAR',y='resíduos IND', title= "CAR vs IND - resíduos") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# Pontos abaixo da reta- indicando que o modelo independente se ajusta melhor 
# do que o CAR em regiões onde não tem estrutura de dependência espacial, modelos
# equivalentes nessas zonas
# Há mais pontos acima da reta- residuos do independente maiores que o do CAR
# Ajuste do Y pelo independente foi mais distante do verdadeiro que o do CAR, 
# pelas medidas vistas (AIC)


  
  