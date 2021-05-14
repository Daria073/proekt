library("psych")
library("dplyr")
library("ggplot2")
library("GGally")
library("pracma")
library("miscTools")
library("klaR")
library("ggpubr")
library(shiny)

yadr1 <- function(x,n,m){
  FF <- matrix(data=NA, nrow = n, ncol = m)
  for (j in 1:m) {
    for (i in 1:n){
      FF[i,j] <- as.vector(t(x[i,]))%*%as.numeric(x[j,])
    }
  }
  return(FF)
}
yadr2 <- function(x,n,m){
  p <- 1
  # p <- 2
  FF <- matrix(data=NA, nrow = n, ncol = m)
  E <- matrix(data=0, nrow = n, ncol = m)
  E <- diag(rep(1,n))
  for (j in 1:m) {
    for (i in 1:n){
      FF[i,j] <- (E[i,j] + as.vector(t(x[i,]))%*%as.numeric(x[j,]))^p
    }
  }
  return(FF)
}
yadr3 <- function(x,n,m){
  gamma <- 1
  FF <- matrix(data=NA, nrow = n, ncol = m)
  for (j in 1:m) {
    for (i in 1:n){
      norma <- sum((x[i,]-x[j,])*(x[i,]-x[j,]))
      FF[i,j] <- exp((-gamma/2)*(norma)) 
    }
  }
  return(FF)
}

train <- function(x,FF,t){
  n <- nrow(x)
  m <- n
  wmp <- t
  gamma <- 1
  alfa <- rep(1,m)
  s <- c()
  alfabound <- 10^12
  weifgtbound <- 10^(-6)
  numofiterations <- 100
  l =TRUE
  for (k in 1:numofiterations) { 
    A <- diag(alfa)
    while (l == TRUE) {
      for (i in 1:n) {
        s[i] <- 1/(1+exp(t[i]*(wmp%*%t(FF[i,])))) 
      }
      R <- diag(c(s*(1-s))) 
      z <- as.vector(FF%*%wmp + solve(R)%*%(s-t)) 
      sigma <- solve((t(FF)%*%R%*%FF)+A)
      wmp.new <- as.vector((sigma%*%t(FF)%*%R)%*%z)  
      #wmp <- wmp.new
      if (Norm(wmp.new - wmp)<0.1) {
        l = FALSE
      } else {
        wmp <- wmp.new
      }
    }
    for (j in 1:m){
      if ((abs(wmp[j])<weifgtbound) || (alfa[j]>alfabound)){
        wmp[j] <- 0
        alfa[j] <- 10^13
        gamma <- 0 
      } else {
        alfa[j] <- (1-alfa[j]*sigma[j,j])/(wmp[j]^2)
      }
    }
  }
  return(wmp)
}

wplot <- function(wmp,m){
  wmp.0 <- wmp 
  for (j in 1:m){
    if (wmp.0[j] != 0) {
      wmp.0[j] = 1
    } else {
      wmp.0[j] = NA
    }
  }
  return(wmp.0)
}

vttrain <- function(wmp,FF,m){
  t.tr <- c()
  for (j in 1:m){
    t.tr[j] <- wmp%*%FF[,j]
    if (t.tr[j]<0) {
      t.tr[j] <- 1 
    } else {
      t.tr[j] <- -1
    }
  }
  return(t.tr)
}

vtvertrain <- function(wmp,FF,m){
  t.itog <- c()
  for (j in 1:m){
    t.itog[j] <- wmp%*%FF[,j]
    if (t.itog[j]<0) {
      t.itog[j] <- 1 
    } else {
      t.itog[j] <- -1
    }
  }
  return(t.itog)
}

rak <- read.table("data/rak.data", head = FALSE, as.is = TRUE, sep =',')
rak1 <- rak
rak <- rak[3:12]
ill<-read.table("data/index.csv", head=TRUE, sep = ";")
ill <- na.omit(ill)
ill1 <- ill
ill <- ill[2:10]

shinyServer(function(input, output) {
  
  datasetInput <- reactive({
    switch(input$dataset,
           "rak" = rak,
           "ill" = ill)
  })
  NormaInput <- reactive({
    switch(input$normalise, "normalise" = {data.frame(scale(datasetInput()))})
  })
  nInput <- reactive({ nrow(datasetInput()[input$obv[1]:input$obv[2],]) })
  mInput <- reactive({ nrow(datasetInput()[input$obv[1]:input$obv[2],]) })
  yaderInput <- reactive({
    switch(input$yader,
           "Linear_core" = yadr1(NormaInput()[input$obv[1]:input$obv[2],],nInput(),mInput()),
           "Polynomial_core" = yadr2(NormaInput()[input$obv[1]:input$obv[2],],nInput(),mInput()),
           "Gaussian_core" = yadr3(NormaInput()[input$obv[1]:input$obv[2],],nInput(),mInput()))
  })
  tInput <- reactive({ 
    if (input$dataset == "rak") {
    rak1$V2[input$obv[1]:input$obv[2]]
  } else {
    ill1$X11[input$obv[1]:input$obv[2]]
  }
    })
  verInput <- reactive({
    NormaInput()[-(input$obv[1]:input$obv[2]),]
  })
  tverInput <- reactive ({
    if (input$dataset == "rak") {
      rak1$V2[-(input$obv[1]:input$obv[2])]
    } else {
      ill1$X11[-(input$obv[1]:input$obv[2])]
    }
  })
  yaderverInput <- reactive({
    switch(input$yader,
           "Linear_core" = yadr1(verInput(),nInput(),length(tverInput())),
           "Polynomial_core" = yadr2(verInput(),nInput(),length(tverInput())),
           "Gaussian_core" = yadr3(verInput(),nInput(),length(tverInput())))
  })

  output$view <- renderTable({
    if(input$normalise == TRUE){
      head(NormaInput()[input$obv[1]:input$obv[2],])
    } else{
      head(datasetInput()[input$obv[1]:input$obv[2],])
    }
  })
  output$tr <- renderText ({
    paste("Данная выборка имеет ", nrow(datasetInput()), " наблюдений.")
  })
  
  output$plot1 <- renderPlot({
    #if (input$result == TRUE){
    if(input$normalise == TRUE){
    wmp <- train(NormaInput()[input$obv[1]:input$obv[2],],yaderInput(),tInput())
    pp <- data.frame(t = tInput(),x = NormaInput()[input$obv[1]:input$obv[2],],wmp = wplot(wmp,mInput()))
    if (input$dataset == "rak") {
      g1 <- ggplot(data = pp, aes(x = rak$V3[input$obv[1]:input$obv[2]],
                                  y = rak$V4[input$obv[1]:input$obv[2]],
                                  shape = factor(t), color = factor(t), size = 0.02)) + geom_point()+ xlab ("Радиус") + ylab ("Текстура")
      g3 <- ggplot(data = pp, aes(x = rak$V3[input$obv[1]:input$obv[2]],
                                  y = rak$V4[input$obv[1]:input$obv[2]],
                                  shape = factor(wmp), color = factor(wmp), size = 0.08)) + geom_point(color = "green")+ xlab ("Радиус") + ylab ("Текстура")
      ggarrange(g1,g3)
    } else {
      g1 <- ggplot(data = pp, aes(x = ill$X5[input$obv[1]:input$obv[2]],
                                  y = ill$X8[input$obv[1]:input$obv[2]],
                                  shape = factor(t), color = factor(t), size = 0.02)) + geom_point()+ xlab ("ALT") + ylab ("Общий белок")
      g3 <- ggplot(data = pp, aes(x = ill$X5[input$obv[1]:input$obv[2]],
                                  y = ill$X8[input$obv[1]:input$obv[2]],
                                  shape = factor(wmp), color = factor(wmp), size = 0.08)) + geom_point(color = "green")+ xlab ("ALT") + ylab ("Общий белок")
      ggarrange(g1,g3)
   # }
    } }
  })  
  output$plot2 <- renderPlot({
   # if (input$result == TRUE){
    if(input$normalise == TRUE){
    wmp <- train(NormaInput()[input$obv[1]:input$obv[2],],yaderInput(),tInput())
    pp <- data.frame(t = tInput(),x = NormaInput()[input$obv[1]:input$obv[2],],t.tr = vttrain(wmp,yaderInput(),mInput()))
    if (input$dataset == "rak") {
        g1 <- ggplot(data = pp, aes(x = rak$V3[input$obv[1]:input$obv[2]],
                                    y = rak$V4[input$obv[1]:input$obv[2]],
                                    shape = factor(t), color = factor(t), size = 0.02)) + geom_point()+ xlab ("Радиус") + ylab ("Текстура")
        g3 <- ggplot(data = pp, aes(x = rak$V3[input$obv[1]:input$obv[2]],
                                    y = rak$V4[input$obv[1]:input$obv[2]],
                                    shape = factor(t.tr), color = factor(t.tr), size = 0.04)) + geom_point()+ xlab ("Радиус") + ylab ("Текстура")
        ggarrange(g1,g3)
      } else {
        g1 <- ggplot(data = pp, aes(x = ill$X5[input$obv[1]:input$obv[2]],
                                    y = ill$X8[input$obv[1]:input$obv[2]],
                                    shape = factor(t), color = factor(t), size = 0.02)) + geom_point()+ xlab ("ALT") + ylab ("Общий белок")
        g3 <- ggplot(data = pp, aes(x = ill$X5[input$obv[1]:input$obv[2]],
                                    y = ill$X8[input$obv[1]:input$obv[2]],
                                    shape = factor(t.tr), color = factor(t.tr), size = 0.04)) + geom_point()+ xlab ("ALT") + ylab ("Общий белок")
        ggarrange(g1,g3)
    #  }
    } }
  })
  output$plot3 <- renderPlot({
    #if (input$result == TRUE){
    if(input$normalise == TRUE){
    wmp <- train(NormaInput()[input$obv[1]:input$obv[2],],yaderInput(),tInput())
    pp <- data.frame(t = tverInput(),x = verInput(),t.itog = vtvertrain(wmp,yaderverInput(),length(tverInput())))
     if (input$dataset == "rak") {
        g1 <- ggplot(data = pp, aes(x = rak$V3[-(input$obv[1]:input$obv[2])],
                                    y = rak$V4[-(input$obv[1]:input$obv[2])],
                                    shape = factor(t), color = factor(t))) + geom_point()+ xlab ("Радиус") + ylab ("Текстура")
        g3 <- ggplot(data = pp, aes(x = rak$V3[-(input$obv[1]:input$obv[2])],
                                    y = rak$V4[-(input$obv[1]:input$obv[2])],
                                    shape = factor(t.itog), color = factor(t.itog))) + geom_point()+ xlab ("Радиус") + ylab ("Текстура")
        ggarrange(g1,g3)
      } else {
        g1 <- ggplot(data = pp, aes(x = ill$X5[-(input$obv[1]:input$obv[2])],
                                    y = ill$X8[-(input$obv[1]:input$obv[2])],
                                    shape = factor(t), color = factor(t))) + geom_point()+ xlab ("ALT") + ylab ("Общий белок")
        g3 <- ggplot(data = pp, aes(x = ill$X5[-(input$obv[1]:input$obv[2])],
                                    y = ill$X8[-(input$obv[1]:input$obv[2])],
                                    shape = factor(t.itog), color = factor(t.itog))) + geom_point()+ xlab ("ALT") + ylab ("Общий белок")
        ggarrange(g1,g3)
    #  }
    } 
    }
  })
  
  output$summary1 <- renderTable({
   # if (input$result == TRUE){
    wmp <- train(NormaInput()[input$obv[1]:input$obv[2],],yaderInput(),tInput())
    Acc0 <- mean(vttrain(wmp,yaderInput(),mInput()) == tInput())
    paste("Точность=",round(100*Acc0,2),"%")
   # }
  })
  output$summary2 <- renderTable({
   # if (input$result == TRUE){
    wmp <- train(NormaInput()[input$obv[1]:input$obv[2],],yaderInput(),tInput())
    (table(Факт = tverInput(), Прогноз = vtvertrain(wmp,yaderverInput(),length(tverInput()))))
    #}
  })
  
  output$summary3 <- renderTable({
   # if (input$result == TRUE){
    wmp <- train(NormaInput()[input$obv[1]:input$obv[2],],yaderInput(),tInput())
    Acc1 <- mean(vtvertrain(wmp,yaderverInput(),length(tverInput())) == tverInput())
    paste("Точность=",round(100*Acc1,2),"%")
   # }
  })
  
})
