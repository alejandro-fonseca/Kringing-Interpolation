## Kriging Interpolation and semivariogram Analysis based on IDSTA data taken by MODIS Satellite Image
#Croatia - 2008



library(Rcmdr)
library(RcmdrMisc)
library(sp)


setwd("C:/Users")
load("IDSTA.ov.rda") 

attributes(IDSTA.ov)
HDEM<-IDSTA.ov$HRdem
Dsea<-IDSTA.ov$HRdsea
TWI<-IDSTA.ov$HRtwi
latitud<-IDSTA.ov$Lat
longitud<-IDSTA.ov$Lon
LST<-IDSTA.ov$LST2008_02_02
BASE_DATOS <- data.frame(longitud,latitud,LST,HDEM, Dsea, TWI)
 
#remover NA
BASE_DATOS2 <- na.omit(BASE_DATOS)

#convertir dataframe a geodata y graficar
library(geoR)
base_datos.geo<-as.geodata(BASE_DATOS2)
plot.geodata(base_datos.geo) 


#NORMALIDAD Y TENDENCIA
#ANALISIS DE NORMALIDAD  DATOS
library(mvnormtest)
library(stats)
library(MVN) 
#graficos
hist(BASE_DATOS2$LST)
boxplot(BASE_DATOS2$LST)
ks.test(BASE_DATOS2$LST,pnorm)
#prueba de normalidad multivariante
mshapiro.test(t(BASE_DATOS2$LST))
hzTest(BASE_DATOS2, cov = TRUE, qqplot = FALSE )

#Corregir normalidad
library(RcmdrPlugin.epack)
lambda<-bc2(BASE_DATOS2$LST+1)
LST.T<-(((BASE_DATOS2$LST)^lambda )-1 )/lambda
#2da prueba 
mshapiro.test(t(LST.T))
#NO FUNCIONO
          
#Correcion de normalidad anamorfosis gaussiana
library(RGeostats)
library(gstat)
mdb = BASE_DATOS2[,c("longitud","latitud","LST")] 
mdb.rgdb = db.create(mdb,ndim=2,autoname=F)
mdb.herm = anam.fit(mdb.rgdb,name="LST",type="gaus")
mdb.hermtrans = anam.z2y(mdb.rgdb,names="LST",anam=mdb.herm)
LST.TG = mdb.hermtrans@items$Gaussian.LST 
#Sin transformar
mshapiro.test(t(BASE_DATOS2$LST))
#Transformación Box-Cox
mshapiro.test(t(LST.T))
#Con transformación: Anamorfosis Gaussiana
mshapiro.test(t(LST.TG))

#AGREGAR DATOS NORMALES A BASE DE DATOS
BASE_DATOS.T<-data.frame(BASE_DATOS2$longitud,BASE_DATOS2$latitud,LST.TG)

#ANALISIS DE TENDENCIA
x11()
par(mfrow=c(2,3))
plot(longitud,LST,pch=20,xlab="longitud",ylab="Temperatura", xlim=c(12,20))
plot(latitud,LST,pch=20,xlab="Latitud",ylab="Temperatura", ylim=c(0,12))
plot(HDEM,LST,pch=20,xlab="Modelo Digital de elevación",ylab="Temperatura", xlim=c(1,1453))
plot(Dsea,LST,pch=20,xlab="Distancia costera",ylab="Temperatura", xlim=c(0,253)) 
plot(TWI,LST,pch=20,xlab="Indice de humedad topografico",ylab="Temperatura", xlim=c(12,22))

     
# Calculo de superficies de tendencia
library(spatial)             
basedatos.ls <- surf.ls(4,BASE_DATOS2$longitud,BASE_DATOS2$latitud,BASE_DATOS2$LST)
bd.trsurf <- trmat(basedatos.ls, 10, 20, 42, 47, 1)

# Representación de superficies de tendencia
par(pty="s",mar=c(2,2,2,2))
contour(bd.trsurf)
points(longitud,latitud,pch=20)
par(mar=c(0,0,0,0))
image(bd.trsurf)
points(longitud,latitud,pch=20)
par(mfrow=c(1,1))
persp(bd.trsurf)
persp(bd.trsurf,theta=210,phi=20,col=2,ltheta=-50,shade=0.25,xlab="longitud",ylab="latitud")
         
#Analisis estadistico de la tendencia
attach(BASE_DATOS2)
summary(lm(LST~longitud+latitud+HDEM+ Dsea+ TWI))
MODELO<-lm(LST~longitud+latitud+HDEM)
summary(MODELO)
       
#ANISOTROPIA 
install.packages("intamap")
library(intamap)
library(sp)

xy<- base_datos.geo$coords
z<-base_datos.geo$data
pts = data.frame(xy,LST.TG,BASE_DATOS2$HDEM, z=z)
coordinates(pts) <- ~longitud+latitud
SpatialPoints(coordinates(pts))
class(pts)
estimateAnisotropy(pts)
# r = 1 y un theta = 0. Por lo tanto no hay Anisotropia

#CALCULAR VARIOGRAMA
library(gstat)
attach(pts)
#Experimental robusto
t.svgm <- variogram(LST.TG~longitud+latitud+pts$BASE_DATOS2.HDEM,data= pts, cressie=TRUE)
#Experimental clasico
t.svgm2 <- variogram(LST.TG~longitud+latitud+pts$BASE_DATOS2.HDEM,data= pts, cressie=FALSE)
plot(t.svgm$dist, t.svgm$gamma,xlim=c(0,2.5), xlab="h", ylab="Gamma", main="Robusto" )
plot(t.svgm2$dist, t.svgm2$gamma,xlim=c(0,2.5), xlab="h", ylab="Gamma", main="Clasico" )
#esferico
ESF<-vgm(0.5,"Sph",0.7, 0)
#exponencial
EXP<-vgm(0.5,"Exp",0.7, 0)
#estable
EST<-vgm(0.5,"Exc",0.7, 0)
plot(t.svgm,ESF, xlab="h", ylab="Gamma", main="Esferico")
plot(t.svgm,EXP, xlab="h", ylab="Gamma", main="Exponencial")
plot(t.svgm,EST, xlab="h", ylab="Gamma", main="Estable")

 #METODO DE AJUSTE FUNCION GSTAT           
Esf.gstat <- gstat(id="LST.TG", formula=LST.TG~longitud+latitud+HDEM, model=ESF, data=pts)
Exp.gstat <- gstat(id="LST.TG", formula=LST.TG~longitud+latitud+HDEM, model=EXP, data=pts)
Est.gstat <- gstat(id="LST.TG", formula=LST.TG~longitud+latitud+HDEM, model=EST, data=pts)

#Metodo de ajuste Minimos cuadrados ordianarios

ESF.ols <- fit.variogram(t.svgm, vgm(0.5, "Sph", 0.7, 0),fit.method = 6)
EXP.ols <- fit.variogram(t.svgm, vgm(0.5, "Exp", 0.7, 0),fit.method = 6)
EST.ols <- fit.variogram(t.svgm, vgm(0.5, "Exc", 0.7, 0),fit.method = 6)

ESF.A<-vgm(0.367,"Sph",0.4814, 0.11)
EXP.A<-vgm(0.4811,"Exp",0.1577, 0)
EST.A<-vgm(0.4963,"Exc",0.1144, 0)

plot(t.svgm,ESF.A, xlab="h", ylab="Gamma", main="Esferico Ajustado")
plot(t.svgm,EXP.A, xlab="h", ylab="Gamma", main="Exponencial Ajustado")
plot(t.svgm,EST.A, xlab="h", ylab="Gamma", main="Estable Ajustado")

Dist<-t.svgm$dist
ESF.OLS <- variogramLine(ESF.A, min=0, dist_vector=t.svgm$dist)
EXP.OLS <- variogramLine(EXP.A, min=0, dist_vector=t.svgm$dist)
EST.OLS <- variogramLine(EST.A, min=0, dist_vector=t.svgm$dist)

#RESIDUOS
 resi.ESF <- sum((t.svgm$gamma-ESF.OLS$gamma)^2)
 resi.EXP <- sum((t.svgm$gamma-EXP.OLS$gamma)^2)
 resi.EST <- sum((t.svgm$gamma-EST.OLS$gamma)^2)

