library(shiny)
source('model.R', local=TRUE)
cols <- c(S='#0072B2',E='#B57900',I='#D55E00',R='#009E73',D='#595959')
ui <- fluidPage(
  tags$head(tags$style(HTML('body{background:#f5f7fa;color:#193047;font-size:16px} .container-fluid{max-width:1400px} .well,.tab-content{background:white;border-radius:10px} .tab-content{padding:20px} h1,h2,h3{color:#123954} .btn-primary{background:#156078} .note{padding:14px;background:#eaf3f7;border-left:4px solid #156078;margin:12px 0} pre{white-space:pre-wrap} .shiny-plot-output{max-width:100%}'))),
  titlePanel('When does an infection spread?'),
  p('ENSCI 3120 | Predict, run, compare, explain.'),
  sidebarLayout(sidebarPanel(width=3,
    selectInput('preset','Classroom scenario',c('Baseline','Lower transmission','Faster recovery')),
    actionButton('load','Load scenario settings'),
    helpText('Loading a scenario changes the controls. Predict first, then press Run.'),
    selectInput('model','Model',c('SIR: lasting immunity'='SIR','SIS: susceptible again'='SIS',
      'SEIR: exposed before infectious'='SEIR','SIRD: disease mortality'='SIRD')),
    numericInput('S','Initially susceptible (S)',99,min=0,max=2000),
    numericInput('I','Initially infectious (I)',1,min=0,max=2000),
    conditionalPanel("input.model == 'SEIR'",numericInput('E','Initially exposed (E)',0,min=0,max=2000)),
    conditionalPanel("input.model != 'SIS'",numericInput('R','Initially recovered / immune (R)',0,min=0,max=2000)),
    numericInput('beta','Transmission coefficient (beta)',0.004,min=0,max=.02,step=.0005),
    numericInput('gamma','Recovery rate (gamma, per day)',.2,min=0,max=2,step=.05),
    conditionalPanel("input.model == 'SEIR'",numericInput('a','Progression E to I (a, per day)',.25,min=0,max=2,step=.05)),
    conditionalPanel("input.model == 'SIRD'",numericInput('mu','Disease mortality (mu, per day)',.05,min=0,max=2,step=.01)),
    tags$details(tags$summary('Births, deaths, and vaccination'),
      numericInput('m','Background birth/death rate (m, per day)',0,min=0,max=2,step=.01),
      conditionalPanel("input.model != 'SIS'",numericInput('v','Fraction of newborns vaccinated (v)',0,min=0,max=1,step=.1)),
      helpText('Births balance background deaths. v acts only on newborns; when m = 0 it has no effect.')),
    sliderInput('days','Duration (days)',min=1,max=180,value=60,step=1),
    actionButton('run','Run simulation',class='btn-primary'),
    tags$hr(),actionButton('save','Keep current run for comparison'),
    actionButton('clear','Clear comparison'),
    downloadButton('csv','Download current data'),
    helpText('The downloaded data include the settings used for the displayed run.')
  ), mainPanel(width=9,
    div(class='note',textOutput('status')),
    tabsetPanel(
      tabPanel('Investigate',h3('Make a prediction'),
        p('Will the infectious group initially grow or shrink? Which flow should be larger?'),
        textAreaInput('prediction','Your prediction (optional; stays in this session)',rows=2,width='100%'),
        plotOutput('timeplot',height='380px'),
        verbatimTextOutput('summary'),
        h3('What changes the infectious group?'),
        plotOutput('flowplot',height='280px'),
        p('Compare entries to and exits from I. A falling infectious curve does not necessarily mean that transmission has stopped.'),
        tags$details(tags$summary('Show the threshold and initial growth explanation'),verbatimTextOutput('threshold'))),
      tabPanel('Phase plot',p('Follow the trajectory through susceptible and infectious host counts. The filled circle marks the start.'),
        plotOutput('phase',height='430px'),p('The dashed vertical line is the susceptible threshold when it falls within the graph. For SEIR this is a transmission-chain threshold, not the boundary for instantaneous growth of I.')),
      tabPanel('Model and assumptions',uiOutput('equations'),
        h3('What the symbols mean'),
        tags$ul(tags$li('S: susceptible; I: infectious; R: recovered or otherwise immune.'),
          tags$li('E: exposed but not yet infectious (SEIR). D: cumulative disease deaths (SIRD).'),
          tags$li('N is the living population. Transmission is beta*S*I, not beta*S*I/N.'),
          tags$li('beta has units hosts^-1 day^-1 in a fixed study area. All other rates have units day^-1.'),
          tags$li('Homogeneous mixing, constant rates, and continuous host counts. This is a teaching model, not a forecast of a particular disease.'),
          tags$li('Background births enter S or, if vaccinated, R. D does not count background deaths.'),
          tags$li('No disease-induced death except in SIRD. Immunity persists except in SIS.')),
        div(class='note','Notation bridge: the reading uses m for total loss from the infected group. Here gamma is recovery, m is background mortality, and mu is additional disease mortality. For the classroom SIR example, m = 0 and the reading\'s loss rate equals gamma.'),
        p('Original R implementation based on the equations and model sequence in ',
          a('EcoEvoApps',href='https://ecoevoapps.gitlab.io/docs/articles/sir-dynamics.html',target='_blank'),
          '. Runs locally in your browser when published with Shinylive.')),
      tabPanel('Classroom activity',h3('1. Baseline'),p('Load Baseline. Sketch a prediction, then run. When I starts to decline, have new infections stopped? Use the two flow curves.'),
        h3('2. Lower transmission'),p('Keep the baseline run for comparison. Load Lower transmission. Change only beta from 0.004 to 0.0015. Predict, run, and explain the difference.'),
        h3('3. Faster recovery'),p('Load Faster recovery. Compare with the baseline: gamma changes from 0.20 to 0.40 per day. What happens to the susceptible threshold and mean infectious duration?'),
        h3('Extensions'),p('Return to baseline and change the model. How do loss of immunity, a latent period, or disease mortality change your interpretation? Keep track of which assumptions changed.'),
        p('This app does not submit answers to Top Hat. Submit your explanation there when prompted.'))
    )
  ))
)
server <- function(input,output,session) {
  controls <- reactive(canonical_config(list(model=input$model,S=input$S,E=input$E,I=input$I,R=input$R,
    beta=input$beta,gamma=input$gamma,m=input$m,v=input$v,a=input$a,mu=input$mu,days=input$days)))
  current <- reactiveVal(NULL); previous <- reactiveVal(NULL); error <- reactiveVal(NULL)
  observeEvent(input$load, {
    p <- default_config()
    if(input$preset=='Lower transmission') p$beta <- .0015
    if(input$preset=='Faster recovery') p$gamma <- .4
    updateSelectInput(session,'model',selected=p$model)
    for(n in setdiff(names(p),c('model','days'))) updateNumericInput(session,n,value=p[[n]])
    updateSliderInput(session,'days',value=p$days)
  })
  observeEvent(input$run, {
    result <- tryCatch(simulate_model(controls()),error=function(e)e)
    if(inherits(result,'error')) { error(conditionMessage(result)); return() }
    current(result); error(NULL)
  })
  observeEvent(input$save,{ req(current()); previous(current()) })
  observeEvent(input$clear,previous(NULL))
  output$status <- renderText({
    if(!is.null(error())) return(paste('Unable to run:',error(),'Any graphs shown are from the last successful run.'))
    if(is.null(current())) return('Choose settings and make a prediction. Press Run simulation to see the results.')
    if(!identical(controls(),current()$config)) return('Settings have changed. Graphs still show the last completed run. Press Run to update them.')
    'Graphs show the completed run. Keep this run to compare it with your next experiment.'
  })
  result <- reactive({ validate(need(!is.null(current()),'Make your prediction, then press Run simulation.')); current() })
  label <- function(r) {
    p <- r$config
    sprintf('%s: beta=%g, gamma=%g, m=%g, v=%g; S=%g, I=%g, E=%g, R=%g; a=%g, mu=%g; %g days',
      p$model,p$beta,p$gamma,p$m,p$v,p$S,p$I,p$E,p$R,p$a,p$mu,p$days)
  }
  output$summary <- renderText({
    r <- result(); d <- r$data; peak <- which.max(d$I)
    paste(label(r),sprintf('\nPeak infectious on sampled curve: %.2f hosts at day %.2f.\nFinal infectious: %.3f hosts.',d$I[peak],d$time[peak],tail(d$I,1)),
      if(!is.null(previous())) paste('\nComparison (dashed):',label(previous())) else '')
  })
  output$timeplot <- renderPlot({
    r <- result(); old <- previous(); d <- r$data; states <- r$states
    ymax <- max(d[states],if(!is.null(old)) old$data[old$states] else 0)
    xmax <- max(d$time,if(!is.null(old)) old$data$time else 0)
    par(mar=c(4,4,2,1)); matplot(d$time,d[states],type='l',lty=1,lwd=3,col=cols[states],
      xlim=c(0,xmax),ylim=c(0,max(1,ymax)),xlab='Time (days)',ylab='Hosts',main='Host states over time')
    if(!is.null(old)) for(s in old$states) lines(old$data$time,old$data[[s]],col=cols[s],lty=2,lwd=2)
    shown <- union(states,if(!is.null(old)) old$states else character())
    legend('topright',legend=shown,col=cols[shown],lty=1,lwd=3,bg='white')
  },res=110)
  output$flowplot <- renderPlot({
    r <- result(); d <- r$data
    par(mar=c(4,4,2,1)); matplot(d$time,d[c('infectious_entries','infectious_exits')],type='l',
      col=c('#0072B2','#D55E00'),lty=c(1,2),lwd=3,xlab='Time (days)',ylab='Hosts per day',main='Entries to and exits from I')
    legend('topright',legend=c(if(r$config$model=='SEIR') 'E to I: a*E' else 'New infections: beta*S*I','All exits from I'),
      col=c('#0072B2','#D55E00'),lty=c(1,2),lwd=3,bg='white')
  },res=110)
  output$phase <- renderPlot({
    r <- result(); d <- r$data; old <- previous()
    xs <- c(d$S,if(!is.null(old)) old$data$S); ys <- c(d$I,if(!is.null(old)) old$data$I)
    par(mar=c(4,4,2,1)); plot(d$S,d$I,type='l',lwd=3,col=cols['I'],
      xlim=c(0,max(1,xs)),ylim=c(0,max(1,ys)),xlab='Susceptible hosts',ylab='Infectious hosts')
    points(d$S[1],d$I[1],pch=16,col=cols['I'],cex=1.5)
    if(!is.null(old)) lines(old$data$S,old$data$I,lty=2,lwd=2,col='grey40')
    st <- threshold(r$config); if(is.finite(st)) abline(v=st,lty=2,col=cols['S'])
  },res=110)
  output$threshold <- renderText({
    r <- result(); p <- r$config; st <- threshold(p)
    if(p$model=='SEIR') return(sprintf('SEIR: instantaneous dI/dt = a*E - (gamma+m)*I = %.4f hosts/day initially.\nTransmission-chain threshold S* = (gamma+m)*(a+m)/(beta*a) = %s.\nAn initial fall in I can precede an outbreak as exposed hosts become infectious.',r$initial_slope,format(st,digits=5)))
    sprintf('For I > 0, infectious hosts increase when beta*S > gamma+m+mu.\nS* = (gamma+m+mu)/beta = %s susceptible hosts.\nInitial dI/dt = %.4f hosts/day.\nIf I = 0, no infection is introduced spontaneously. With beta = 0 there is no transmission.',format(st,digits=5),r$initial_slope)
  })
  output$equations <- renderUI({
    r <- result(); p <- r$config
    eq <- switch(p$model,
      SIR='dS/dt = m*N*(1-v) - m*S - beta*S*I\ndI/dt = beta*S*I - (gamma+m)*I\ndR/dt = gamma*I - m*R + m*v*N',
      SIS='dS/dt = m*N - m*S - beta*S*I + gamma*I\ndI/dt = beta*S*I - (gamma+m)*I',
      SEIR='dS/dt = m*N*(1-v) - m*S - beta*S*I\ndE/dt = beta*S*I - (a+m)*E\ndI/dt = a*E - (gamma+m)*I\ndR/dt = gamma*I - m*R + m*v*N',
      SIRD='dS/dt = m*N*(1-v) - m*S - beta*S*I\ndI/dt = beta*S*I - (gamma+m+mu)*I\ndR/dt = gamma*I - m*R + m*v*N\ndD/dt = mu*I')
    tagList(h3(paste('Displayed run:',p$model)),tags$pre(eq))
  })
  output$csv <- downloadHandler(filename=function() 'infection-simulation.csv',content=function(file) {
    r <- result(); d <- r$data
    for(n in names(r$config)) d[[paste0('setting_',n)]] <- r$config[[n]]
    write.csv(d,file,row.names=FALSE)
  })
}
shinyApp(ui,server)
