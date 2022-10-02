packages <- c("ggplot2", "dplyr", "plotly", "tidyverse", "shiny")

install.packages(setdiff(packages, rownames(installed.packages()))) 

library(shiny)
library(ggplot2)
library(plotly)
library(tidyverse)
library(dplyr)

##########phase response curve######################################
pprc <- function(k,phic) {
  
  phi.values = seq(0, 1, by = 0.001)
  dat = data.frame(x=phi.values)
  dat$y = with(dat, ifelse((0<=x & x < phic),k*x+1, ifelse((phic< x & x<1),k*x+1-k, NA)))
  p <- ggplot(dat, aes(x,y)) + geom_line() + 
    labs(y="f(phi) = T/te",x="phi")
  
  return(ggplotly(p))
}

##########cobweb plot###############################################
p_cobweb<-function(ts,te,theta,k,phic,phi_o,N)
{
  par(mfrow=c(2,1))
  
  ######relationship between phase and time of sinus beat
  x <- rep(NA,N)
  x[1] <- phi_o
  for(i in 2:N)
  {
    if (0 <= x[i-1] & x[i-1] <(ts-theta)/te)
    {
      x[i]<-(x[i-1]+ts/te) %% 1
    }
    else if (x[i-1] >= (ts-theta)/te & x[i-1] < phic )
    {
      x[i] <- (x[i-1] - ((k*x[i-1]) +1) + ts/te) %% 1
    }
    else if (x[i-1] >= phic & x[i-1] <1)
    {
      x[i] <- (x[i-1] - (k*(x[i-1]-1) + 1) + ts/te) %% 1
    }
  }
  
  plot(x[1:N],type='l',xlab='t',ylab='phi')
  
  ######cobweb plot 
  #initial phi_{i}
  x <- seq(0,1,0.001)
  #create an empty vector to store value of phi_{i+1}
  x_np1<-rep(NA,length(x))
  #phase of the (i+1) th sinus beat function:
  for (i in 1:length(x))
  {
    if (0 <= x[i] & x[i] <(ts-theta)/te)
    {
      x_np1[i]<-(x[i]+ts/te) %% 1
    } else if (x[i] >= (ts-theta)/te & x[i] < phic )
    {
      x_np1[i] <- (x[i] - ((k*x[i]) +1) + ts/te) %% 1
    } else if (x[i] >= phic & x[i] <1)
    {
      x_np1[i] <- (x[i] - (k*(x[i]-1) + 1) + ts/te) %% 1
    }
  }
  
  plot(x,x_np1,xlab="phi_i",ylab="phi_i+1",ylim=c(0,1),xlim=c(0,1))
  lines(c(0,1),c(0,1),type='l',lty=2)
  
  start=phi_o
  bound=FALSE
  
  #sinus phase function
  fsinus <- function(start){
    if (0 <= start & start <(ts-theta)/te)
    {
      end <-(start+ts/te) %% 1
    } else if (start >= (ts-theta)/te & start < phic )
    {
      end <- (start - ((k*start) +1) + ts/te) %% 1
    } else if (start >= phic & start <1)
    {
      end <- (start - (k*(start-1) + 1) + ts/te) %% 1
    }
    return(end) 
  }
  
  #the first y point
  lines(x=c(start,start),y=c(0,fsinus(start)),col='blue')
  #iteration to get cobweb plot
  for(i in 1:(2*N))
  {
    if(bound)
    {
      lines(x=c(start,start),y=c(start,fsinus(start)),col='blue')
      bound=FALSE
    }
    else
    {
      #this point goes to y=x line
      lines(x=c(start,fsinus(start)),
            y=c(fsinus(start),fsinus(start)),
            col='blue')
      bound=TRUE
      start=fsinus(start)
    }
  }
}

##############shiny app ########################################################
ui <- fluidPage(
  tags$head(
    # to display latex value in shiny app
    tags$link(rel="stylesheet", 
              href="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.css", 
              integrity="sha384-dbVIfZGuN1Yq7/1Ocstc1lUEm+AT+/rCkibIcC/OmWo5f0EA48Vf8CytHzGrSwbQ",
              crossorigin="anonymous"),
    HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/katex.min.js" integrity="sha384-2BKqo+exmr9su6dir+qCw08N2ZKRucY4PrGQPPWU1A7FtlCGjmEGFqXCv5nyM5Ij" crossorigin="anonymous"></script>'),
    HTML('<script defer src="https://cdn.jsdelivr.net/npm/katex@0.10.1/dist/contrib/auto-render.min.js" integrity="sha384-kWPLUVMOks5AQFrykwIup5lo0m3iMkkHrD0uJ4H5cjeGihAutqP0yW0J6dpFiVkI" crossorigin="anonymous"></script>'),
    HTML('
    <script>
      document.addEventListener("DOMContentLoaded", function(){
        renderMathInElement(document.body, {
          delimiters: [{left: "$", right: "$", display: false}]
        });
      })
    </script>')
  ),
  titlePanel("QLSC 600 Entrainment zones for modulated parasystole"),
  helpText('phase of the (i+1)th sinus beat is given by'),
  helpText('$\\phi_{i+1}$ = $\\phi_{i}$ + $\\frac{t_s}{t_e}$ (mod1)    for 0$\\le$ $\\phi_i$ < $\\frac{t_s-\\theta}{te}$'),
  helpText('$\\phi_{i+1}$ = $\\phi_{i}$ - f($\\phi_i$)+ $\\frac{t_s}{t_e}$ (mod1)    for  $\\frac{t_s-\\theta}{te}$ $\\le$ $\\phi_i$ < 1'),
  helpText('The Phase Response Curve is given by:'),
  helpText('f($\\phi$) = k$\\phi$+1 for 0$\\le$ $\\phi$ < $\\phi_{c}$'),
  helpText('f($\\phi$) = k($\\phi$-1)+1 for $\\phi_{c}$ $\\le$ $\\phi$ < 1'),
  helpText('$t_{e}$ = 1.5'),
  helpText('$\\phi_{c}$ = 0.5'),
  helpText('$\\theta$ = 0.4'),
  sidebarLayout(
    sidebarPanel(
      sliderInput("ts", "ts",value=0.64, min = 0.1, max = 1.0, step=0.01),
      sliderInput("k", "k",min = 0, max = 1.0, value = 0.4,step=0.01),
      sliderInput("phio", "phi_o",min = 0.1, max = 1, value = 0.2, step=0.01),
      sliderInput("N", "Iteration",min = 2, max = 1000, value = 100, step=1)
    ),
    mainPanel(
      plotlyOutput("prcplot", height = "250", width = "250"),
      plotOutput("plot")
    )
  )
)

server <- function(input, output, session) {
  output$prcplot <- renderPlotly(
    pprc(input$k,0.5)
  )
  output$plot <- renderPlot({
    print(input$phio)
    validate(
      need(input$phio != 1, "Please select initial condition less than 1")
    )
    input$newplot
    p_cobweb(input$ts,1.5,0.4,input$k,0.5,input$phio,input$N)
  })
}
shinyApp(ui, server)
