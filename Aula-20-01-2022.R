### Aula sincrona - Estat�stica Espacial I
### 20/01/2022
### Rafael Erbisti


## Exemplo: Dados obtidos a partir do banco de dados Scottish Statistics
## (http://statistics.gov.scot), mas tamb�m est�o inclu�dos no pacote 
## CARBayesdata no R. Dados da regi�o de Greater Glasgow & Clyde. 

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
data(pricedata) # dados de pre�os de im�veis de 270 das 271 zonas
head(pricedata)

# price: pre�o mediano dos im�veis vendidos em 2008 (em milhares de libras)
# crime: taxa de criminalidade por 10 mil habitantes
# rooms: n�mero m�dio de quartos em um im�vel
# sales: n�mero de im�veis vendidos em um ano
# driveshop: tempo m�dio gasto dirigindo at� o centro comercial mais pr�ximo
# type: tipo de propriedade de maior prevalencia em cada �rea - apartamento (flat),
#       terra�o (terraced), geminada (semi) e n�o geminada (detached).


## Analise dos dados

# Histogram do pre�o dos imoveis
ggplot(pricedata)+
  geom_histogram(aes(x = price))

# Histogram do log do pre�o dos imoveis
pricedata <- pricedata |> mutate(logprice = log(pricedata$price))
ggplot(pricedata)+
  geom_histogram(aes(x = logprice))
# variavel real continua, n�o possui assimetria

# log preco vs taxa criminalidade
ggplot(pricedata,aes(x = crime , y = logprice)) + 
  geom_point() +
  labs(x='taxa criminalidade',y='log pre�o') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# parece haver rela��o negativa entre as vari�veis

# log preco vs numero de quartos
ggplot(pricedata,aes(x = rooms , y = logprice, fill=rooms)) + 
  geom_point() +
  labs(x='n�mero de quartos',y='log pre�o') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# grafico adequado

ggplot(pricedata,aes(x = as.factor(rooms) , y = logprice, fill=rooms)) + 
  geom_boxplot() +
  labs(x='n�mero de quartos',y='log pre�o') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# grafico para ver o comportamento das categorias: 
# parece haver rela��o positiva entre as vari�veis

# log preco vs numero de imoveis vendidos
ggplot(pricedata,aes(x = sales , y = logprice)) + 
  geom_point() +
  labs(x='imoveis vendidos',y='log pre�o') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# n�o � poss�vel concluir que h� algum tipo de rela��o

# log preco vs tempo ate o centro comercial
ggplot(pricedata,aes(x = driveshop , y = logprice)) + 
  geom_point() +
  labs(x='tempo at� com�rcio',y='log pre�o') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# n�o � poss�vel concluir que h� algum tipo de rela��o

# log preco vs tipo de imovel
ggplot(pricedata,aes(x = type , y = logprice, fill=type)) + 
  geom_boxplot(col="gray40") +
  labs(x='tipo de im�vel',y='log pre�o') +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# de forma geral, h� diferen�a nos pre�os dos im�veis: mais caras- detached e 
# mais baratas- flat e terrace


## An�lise espacial


# juntando as informacoes numa base unica
pricedata.sp <- merge(x=GGHB.IG, y=pricedata, by="IG", all.x=FALSE)


# O pacote leaflet precisa das coordenadas em lat/long
# Precisamos transformar o sistema de coordenadas
pricedata.sp <- spTransform(pricedata.sp,CRS("+proj=longlat +datum=WGS84 +no_defs"))

# fazendo o mapa do pre�o
colours <- colorNumeric(palette = "YlOrRd", domain = pricedata.sp@data$price)
map1 <- leaflet(data=pricedata.sp) |>
        addTiles() |>
        addPolygons(fillColor = ~colours(price), color="", weight=1,
                    fillOpacity = 0.8) |>
        addLegend(pal = colours, values = pricedata.sp@data$price, opacity = 1,
                  title="Price") |>
        addScaleBar(position="bottomleft")
map1
# mais escuros, maiores os pre�os
# localmente h� �reas que formam clusters em rela��o aos pre�os


## Regressao Nao Espacial

form <- logprice ~ crime + rooms + sales + factor(type) + driveshop
model <- lm(formula=form, data=pricedata.sp@data)
summary(model) # todas as vari�veis significativas
AIC(model) 

# Matriz de vizinhan�a (W)
W.nb <- poly2nb(pricedata.sp, row.names = rownames(pricedata.sp@data))
W.list <- nb2listw(W.nb, style="B")
W.list

# Indice de Moran Global
globalMoran <- moran.test(residuals(model), W.list)
globalMoran
# H0: I=0 versus H1: I!=0
# Estatistica de teste 0.3 e p-valor pequeno -> rejeita H0, indicando que
# h� estrutura de dependencia espacial nos residuos do modelo

# Indice de Moral Local (LISA)
localMoran <- localmoran(residuals(model), W.list)
localMoran

pricedata.sp@data <- cbind(pricedata.sp@data,LISA = localMoran[,4])


ggplot(pricedata.sp@data)+
  geom_histogram(aes(x = LISA))
# as �reas em torno do 0 indicam n�o depend�ncia
# as �reas com LISA maior ou menor indicam depend�ncia


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
# dependencia global n�o t�o alta

## Regressao Espacial

## MODELO CAR

form <- logprice ~ crime + rooms + sales + factor(type) + driveshop
nc_car <- spautolm(formula = form,
                   data = pricedata.sp@data,
                   listw =  W.list,
                   family = "CAR")

summary(nc_car)
# todas as variaveis s�o significativas, exceto driveshop
# rho= 0.14788 (n�o 0 nem 1) dependencia mais local que global
# p-valor proximo de 0, o modelo n�o � independente
# AIC penaliza quanto ao n�mero de par�metros, mais completo que o modelo independente
# AIC menor que o do modelo independente-> melhor ajuste

## MODELO SAR

nc_sar <- spautolm(formula = form,
                   data = pricedata.sp@data,
                   listw =  W.list,
                   family = "SAR")

summary(nc_sar)
# se aproxima do modelo CAR: variaveis s�o significativas, exceto driveshop
# rho= 0.10615 (n�o 0 nem 1) dependencia mais local que global
# p-valor pr�ximo de 0, o modelo n�o � independente
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
# ajustou melhor: Y ajustado pr�ximo do Y verdadeiro
# Mapas parecidos, estruturas de dependencia presentes nas regioes mais escuras,
# efeitos maiores, blocos locais visto nas outras an�lises

# residuos CAR vs SAR
ggplot(pricedata.sp@data,aes(x = r_car , y = r_sar)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col=2, lwd=1.1) +
  labs(x='res�duos CAR',y='res�duos SAR', title= "CAR vs SAR - res�duos") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# mais pontos acima da reta, residuos do SAR maiores que o do CAR
# Ajuste do Y pelo SAR foi mais distante do verdadeiro que o do CAR

# residuos CAR vs IND
ggplot(pricedata.sp@data,aes(x = r_car , y = r_ind)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col=2, lwd=1.1) +
  labs(x='res�duos CAR',y='res�duos IND', title= "CAR vs IND - res�duos") +
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5))
# Pontos abaixo da reta- indicando que o modelo independente se ajusta melhor 
# do que o CAR em regi�es onde n�o tem estrutura de depend�ncia espacial, modelos
# equivalentes nessas zonas
# H� mais pontos acima da reta- residuos do independente maiores que o do CAR
# Ajuste do Y pelo independente foi mais distante do verdadeiro que o do CAR, 
# pelas medidas vistas (AIC)

