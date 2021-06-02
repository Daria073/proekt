library(shiny)
library(shinythemes)
setwd("C:/pole")
rak <- read.table("data/rak.data", head = FALSE, as.is = TRUE, sep =',')

shinyUI(fluidPage( theme = shinytheme("united"),
  titlePanel("Программная система по диагностированию рака молочной железы"),
  
  sidebarLayout(
    sidebarPanel(
      #selectInput("dataset", "Набор данных:",
      #            choices = c("rak")),
      #tags$hr(),
      checkboxInput("normalise", label = "Нормализовать данные", value = FALSE),
      sliderInput("obv", label = "Обучающая выборка:",
                  min = 1, max = 600, value = c(1, 400)),
      
      helpText("Замечание: обучающая выборка не должна превышать всю выборку."),
      
      selectInput("method", label = "Выбери метод:", 
                  choices = c("SVM", "RVM"),
                  selected = "RVM"),
      actionButton("result", "Результат")
      ),

    mainPanel(
      tabsetPanel(type = "tabs",
                   tabPanel("Table", h4("Набор данных:"),tableOutput("view"),
                          textOutput("tr")),
                  tabPanel("Plot", h4("Диаграммы рассеяния классов исходной выборки и после обучения выборки:"),plotOutput("plot2"),
                          h4("Сравнение результатов классификации метода с исходными классами для тестовой выборки:"),plotOutput("plot3")),
                  tabPanel("Summary", h4("Оценка точности для обучающей выборки:"),tableOutput("summary1"),
                           h4("Сравнительная таблица количества фактической и прогнозируемой переменной:"),tableOutput("summary2"),
                           h4("Оценка точности классификации:"),tableOutput("summary3"))
      )
    )
  )
))