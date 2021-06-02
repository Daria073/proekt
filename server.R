library("psych")
library("dplyr")
library("ggplot2")
library("GGally")
library("pracma")
library("miscTools")
library("klaR")
library("ggpubr")
library("e1071")
library("caret")
library(shiny)

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

methodsvm <- function(rak1,t){
  svm1 <- svm(formula = t ~ V3+V4+V5+V6+V7+V8+V9+V10+V11+V12, data = rak1, 
      cross = 10, kernel = "radial", gamma = 0.5, cost = 4, type = "C-classification")
  return(svm1)
}

rak <- read.table("data/rak.data", head = FALSE, as.is = TRUE, sep =',')
rak1 <- rak
rak <- rak[3:12]

shinyServer(function(input, output) {
  
  datasetInput <- reactive({
    rak
  })
  NormaInput <- reactive({
    switch(input$normalise, "normalise" = {data.frame(scale(datasetInput()))})
  })
  nInput <- reactive({ nrow(datasetInput()[input$obv[1]:input$obv[2],]) })
  mInput <- reactive({ nrow(datasetInput()[input$obv[1]:input$obv[2],]) })
  
  yaderInput <- reactive({
    yadr3(NormaInput()[input$obv[1]:input$obv[2],],nInput(),mInput())
  })
  tInput <- reactive({ 
    rak1$V2[input$obv[1]:input$obv[2]]
    })
  verInput <- reactive({
    NormaInput()[-(input$obv[1]:input$obv[2]),]
  })
  tverInput <- reactive ({
      rak1$V2[-(input$obv[1]:input$obv[2])]
  })
  yaderverInput <- reactive({
      yadr3(verInput(),nInput(),length(tverInput()))
  })

  methodInput <- reactive({
    methodsvm(NormaInput()[input$obv[1]:input$obv[2],],tInput())
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
  
  output$plot2 <- renderPlot({
   # if (input$result == TRUE){
    if(input$normalise == TRUE){
    wmp <- train(NormaInput()[input$obv[1]:input$obv[2],],yaderInput(),tInput())
    pp <- data.frame(t = tInput(),x = NormaInput()[input$obv[1]:input$obv[2],],t.tr = vttrain(wmp,yaderInput(),mInput()))
    pp2 <- data.frame(x = NormaInput()[input$obv[1]:input$obv[2],],t = tInput(),s = methodInput()$fitted)
    if (input$method == "RVM") {
      g1 <- ggplot(data = pp, aes(x = rak$V3[input$obv[1]:input$obv[2]],
                                  y = rak$V4[input$obv[1]:input$obv[2]],
                                  shape = factor(t), color = factor(t))) + geom_point()+ xlab ("Радиус") + ylab ("Текстура")
      g3 <- ggplot(data = pp, aes(x = rak$V3[input$obv[1]:input$obv[2]],
                                  y = rak$V4[input$obv[1]:input$obv[2]],
                                  shape = factor(t.tr), color = factor(t.tr))) + geom_point()+ xlab ("Радиус") + ylab ("Текстура")
      ggarrange(g1,g3)
    } else {
      
      g1 <- ggplot(data = pp2, aes(x = rak$V3[input$obv[1]:input$obv[2]],
                                   y = rak$V4[input$obv[1]:input$obv[2]],
                                   shape = factor(t), color = factor(t))) + geom_point()+ xlab ("Радиус") + ylab ("Текстура")
      g3 <- ggplot(data = pp2, aes(x = rak$V3[input$obv[1]:input$obv[2]],
                                   y = rak$V4[input$obv[1]:input$obv[2]],
                                   shape = factor(s), color = factor(s))) + geom_point()+ xlab ("Радиус") + ylab ("Текстура")
      ggarrange(g1,g3)
      # }
    } }
  })
  output$plot3 <- renderPlot({
    #if (input$result == TRUE){
    if(input$normalise == TRUE){
    wmp <- train(NormaInput()[input$obv[1]:input$obv[2],],yaderInput(),tInput())
    pp <- data.frame(t = tverInput(),x = verInput(),t.itog = vtvertrain(wmp,yaderverInput(),length(tverInput())))
    pp2 <- data.frame(x = verInput(),t = tverInput(),s.itog = predict(methodInput(), verInput()))
    if (input$method == "RVM") {
      g1 <- ggplot(data = pp, aes(x = rak$V3[input$obv[1]:input$obv[2]],
                                  y = rak$V4[input$obv[1]:input$obv[2]],
                                  shape = factor(t), color = factor(t))) + geom_point()+ xlab ("Радиус") + ylab ("Текстура")
      g3 <- ggplot(data = pp, aes(x = rak$V3[input$obv[1]:input$obv[2]],
                                  y = rak$V4[input$obv[1]:input$obv[2]],
                                  shape = factor(t.itog), color = factor(t.itog))) + geom_point()+ xlab ("Радиус") + ylab ("Текстура")
      ggarrange(g1,g3)
    } else {
      
      g1 <- ggplot(data = pp2, aes(x = rak$V3[input$obv[1]:input$obv[2]],
                                   y = rak$V4[input$obv[1]:input$obv[2]],
                                   shape = factor(t), color = factor(t))) + geom_point()+ xlab ("Радиус") + ylab ("Текстура")
      g3 <- ggplot(data = pp2, aes(x = rak$V3[input$obv[1]:input$obv[2]],
                                   y = rak$V4[input$obv[1]:input$obv[2]],
                                   shape = factor(s.itog), color = factor(s.itog))) + geom_point() + xlab ("Радиус") + ylab ("Текстура")
      ggarrange(g1,g3)
      # }
    } } 
  })
  
  output$summary1 <- renderTable({
   # if (input$result == TRUE){
    if (input$method == "RVM") {
      wmp <- train(NormaInput()[input$obv[1]:input$obv[2],],yaderInput(),tInput())
      Acc0 <- mean(vttrain(wmp,yaderInput(),mInput()) == tInput())
    } else {
      Acc0 <- mean(methodInput()$fitted == tInput())
      # }
    } 
    paste("Точность=",round(100*Acc0,2),"%")
   # }
  })
  output$summary2 <- renderTable({
   # if (input$result == TRUE){
    if (input$method == "RVM") {
      wmp <- train(NormaInput()[input$obv[1]:input$obv[2],],yaderInput(),tInput())
      (table(Факт = tverInput(), Прогноз = vtvertrain(wmp,yaderverInput(),length(tverInput()))))
    } else {
      (table(Факт = tverInput(), Прогноз = predict(methodInput(), verInput())))
      # }
    } 
  })
  
  output$summary3 <- renderTable({
   # if (input$result == TRUE){
    if (input$method == "RVM") {
      wmp <- train(NormaInput()[input$obv[1]:input$obv[2],],yaderInput(),tInput())
      Acc1 <- mean(vtvertrain(wmp,yaderverInput(),length(tverInput())) == tverInput())
    } else {
      Acc1 <- mean(predict(methodInput(), verInput()) == tverInput())
      # }
    } 
    paste("Точность=",round(100*Acc1,2),"%")
   # }
  })
  
})
