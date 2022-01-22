### Aula sincrona - Estatística Espacial I
### 20/01/2022
### Rafael Erbisti


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
library(leaflet)}

# Carregando os dados 
data(GGHB.IG) # area de estudo dividida em 271 zonas geograficas
head(GGHB.IG)
data(pricedata) # dados de preços de imóveis de 270 das 271 zonas
head(pricedata)

# price: preço mediano dos imóveis vendidos em 2008 (em milhares de libras)
# crime: taxa de criminalidade por 10 mil habitantes
# rooms: número médio de quartos em um imóvel
# sales: número de imóveis vendidos em um ano
# driveshop: tempo médio gasto dirigindo até o centro comercial mais próximo
# type: tipo de propriedade de maior prevalencia em cada área - apartamento (flat),
#       terraço (terraced), geminada (semi) e não geminada (detached).


## Analise dos dados

# Histogram do preço dos imoveis
ggplot(pricedata)+
  geom_histogram(aes(x = price))

# Histogram do log do preço dos imoveis
pricedata <- pricedata |> mutate(logprice = log(pricedata$price))
ggplot(pricedata)+
  geom_histogram(aes(x = logprice))
# variavel real continua, não possui assimetria

# log preco vs taxa criminalidade
ggplot(pricedata,aes(x = crime , y = logprice)) + 
  geom_point() +
  labs(x='taxa criminalidade',y='log preço') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# parece haver relação negativa entre as variáveis

# log preco vs numero de quartos
ggplot(pricedata,aes(x = rooms , y = logprice, fill=rooms)) + 
  geom_point() +
  labs(x='número de quartos',y='log preço') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# grafico adequado

ggplot(pricedata,aes(x = as.factor(rooms) , y = logprice, fill=rooms)) + 
  geom_boxplot() +
  labs(x='número de quartos',y='log preço') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# grafico para ver o comportamento das categorias: 
# parece haver relação positiva entre as variáveis

# log preco vs numero de imoveis vendidos
ggplot(pricedata,aes(x = sales , y = logprice)) + 
  geom_point() +
  labs(x='imoveis vendidos',y='log preço') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# não é possível concluir que há algum tipo de relação

# log preco vs tempo ate o centro comercial
ggplot(pricedata,aes(x = driveshop , y = logprice)) + 
  geom_point() +
  labs(x='tempo até comércio',y='log preço') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# não é possível concluir que há algum tipo de relação

# log preco vs tipo de imovel
ggplot(pricedata,aes(x = type , y = logprice, fill=type)) + 
  geom_boxplot(col="gray40") +
  labs(x='tipo de imóvel',y='log preço') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# de forma geral, há diferença nos preços dos imóveis: mais caras- detached e 
# mais baratas- flat e terrace


## Análise espacial


# juntando as informacoes numa base unica
pricedata.sp <- merge(x=GGHB.IG, y=pricedata, by="IG", all.x=FALSE)


# O pacote leaflet precisa das coordenadas em lat/long
# Precisamos transformar o sistema de coordenadas
pricedata.sp <- spTransform(pricedata.sp,CRS("+proj=longlat +datum=WGS84 +no_defs"))

# fazendo o mapa do preço
colours <- colorNumeric(palette = "YlOrRd", domain = pricedata.sp@data$price)
map1 <- leaflet(data=pricedata.sp) |>
        addTiles() |>
        addPolygons(fillColor = ~colours(price), color="", weight=1,
                    fillOpacity = 0.8) |>
        addLegend(pal = colours, values = pricedata.sp@data$price, opacity = 1,
                  title="Price") |>
        addScaleBar(position="bottomleft")
map1
# mais escuros, maiores os preços
# localmente há áreas que formam clusters em relação aos preços


## Regressao Nao Espacial

form <- logprice ~ crime + rooms + sales + factor(type) + driveshop
model <- lm(formula=form, data=pricedata.sp@data)
summary(model) # todas as variáveis significativas
AIC(model) 

# Matriz de vizinhança (W)
W.nb <- poly2nb(pricedata.sp, row.names = rownames(pricedata.sp@data))
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

pricedata.sp@data <- cbind(pricedata.sp@data,LISA = localMoran[,4])


ggplot(pricedata.sp@data)+
  geom_histogram(aes(x = LISA))
# as áreas em torno do 0 indicam não dependência
# as áreas com LISA maior ou menor indicam dependência


colours <- colorNumeric(palette = "Greens", domain = pricedata.sp@data$LISA)
map2 <- leaflet(data=pricedata.sp) |>
  addTiles() |>
  addPolygons(fillColor = ~colours(LISA), color="", weight=1,
              fillOpacity = 0.8) |>
  addLegend(pal = colours, values = pricedata.sp@data$LISA, opacity = 1,
            title="Price") |>
  addScaleBar(position="bottomleft")
map2
# valores nos grupos indicam indice alto ou baixo -> clusters, dependencia local
# dependencia global não tão alta

## Regressao Espacial

## MODELO CAR

form <- logprice ~ crime + rooms + sales + factor(type) + driveshop
nc_car <- spautolm(formula = form,
                   data = pricedata.sp@data,
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
                   data = pricedata.sp@data,
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

pricedata.sp@data['r_car'] <- nc_car$fit$residuals
pricedata.sp@data['r_sar'] <- nc_sar$fit$residuals
pricedata.sp@data['r_ind'] <- model$residuals

colours <- colorNumeric(palette = "Blues", domain = c(-1,1))
map3 <- leaflet(data=pricedata.sp) |>
  addTiles() |>
  addPolygons(fillColor = ~colours(r_car), color="", weight=1,
              fillOpacity = 0.8) |>
  addLegend(pal = colours, values = pricedata.sp@data$r_car, opacity = 1,
            title="CAR") |>
  addScaleBar(position="bottomleft")
map3

map4 <- leaflet(data=pricedata.sp) |>
  addTiles() |>
  addPolygons(fillColor = ~colours(r_sar), color="", weight=1,
              fillOpacity = 0.8) |>
  addLegend(pal = colours, values = pricedata.sp@data$r_sar, opacity = 1,
            title="SAR") |>
  addScaleBar(position="bottomleft")
map4

map5 <- leaflet(data=pricedata.sp) |>
  addTiles() |>
  addPolygons(fillColor = ~colours(r_ind), color="", weight=1,
              fillOpacity = 0.8) |>
  addLegend(pal = colours, values = pricedata.sp@data$r_ind, opacity = 1,
            title="IND") |>
  addScaleBar(position="bottomleft")
map5

# Residuos do modelo CAR com magnitude menor do que dos outros modelos porque se
# ajustou melhor: Y ajustado próximo do Y verdadeiro
# Mapas parecidos, estruturas de dependencia presentes nas regioes mais escuras,
# efeitos maiores, blocos locais visto nas outras análises

# residuos CAR vs SAR
ggplot(pricedata.sp@data,aes(x = r_car , y = r_sar)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col=2, lwd=1.1) +
  labs(x='resíduos CAR',y='resíduos SAR', title= "CAR vs SAR - resíduos") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# mais pontos acima da reta, residuos do SAR maiores que o do CAR
# Ajuste do Y pelo SAR foi mais distante do verdadeiro que o do CAR

# residuos CAR vs IND
ggplot(pricedata.sp@data,aes(x = r_car , y = r_ind)) + 
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

