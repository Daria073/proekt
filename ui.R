library(shiny)
library(shinythemes)
setwd("C:/mydipl")
rak <- read.table("data/rak.data", head = FALSE, as.is = TRUE, sep =',')
ill<-read.table("data/index.csv", head=TRUE, sep = ";")

shinyUI(fluidPage( theme = shinytheme("united"),
  titlePanel("Relevance Vector Machine"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Выбери набор данных:",
                  choices = c("rak", "ill")),
      tags$hr(),
      checkboxInput("normalise", label = "Нормализовать данные", value = FALSE),
      sliderInput("obv", label = "Обучающая выборка:",
                  min = 1, max = 600, value = c(1, 150)),
      
      helpText("Замечание: обучающая выборка не должна превышать всю выборку."),
      selectInput("yader", label = "Выбери ядровую функцию:", 
                  choices = c("Linear_core", "Polynomial_core",
                                 "Gaussian_core"),
                  selected = "Gaussian_core")
     # actionButton("result", "Результат")
      ),

    mainPanel(
      tabsetPanel(type = "tabs",
                   tabPanel("Table", h4("Набор данных:"),tableOutput("view"),
                          textOutput("tr")),
                  tabPanel("Plot", h4("Диаграммы рассеяния классов и релевантных векторов в проекциях на факторы по обучающей выборке:"),plotOutput("plot1"),
                           h4("Диаграммы рассеяния классов исходной выборки и после обучения выборки:"),plotOutput("plot2"),
                           h4("Сравнение результатов классификации RVM с исходными классами для тестовой выборки:"),plotOutput("plot3")),
                  tabPanel("Summary", h4("Оценка точности для обучающей выборки:"),tableOutput("summary1"),
                           h4("Сравнительная таблица количества фактической и прогнозируемой переменной:"),tableOutput("summary2"),
                           h4("Оценка точности классификации:"),tableOutput("summary3"))
      )
    )
  )
))