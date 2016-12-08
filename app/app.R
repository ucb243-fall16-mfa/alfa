#shiny app
#default loading text
dftext=c("Cat Pee", "Passion Fruit", "Green Pepper", "Mineral", "Smoky", "Cirtrus", "Cat Pee", "Passion Fruit", "Green Pepper", "Mineral","Tropical", "Leafy", "Cat Pee", "Passion Fruit", "Green Pepper", "Mineral", "Grassy", "Flinty", "Cat Pee", "Passion Fruit", "Green Pepper", "Mineral", "Leafy", "Cat Pee", "Passion Fruit", "Green Pepper", "Mineral", "Vegetal", "Hay", "Cat Pee", "Passion Fruit", "Green Pepper", "Mineral", "Melon", "Cat Pee", "Passion Fruit", "Green Pepper", "Mineral", "Cat Pee", "Passion Fruit", "Green Pepper", "Mineral", "Grass", "Smoky", "Cat Pee", "Passion Fruit", "Green Pepper", "Mineral", "Peach", "Cat Pee", "Passion Fruit", "Green Pepper", "Mineral")
dftext=paste(dftext,collapse = ',')
#load our mfa package
install.packages("mfa_0.1.tar.gz", repos = NULL, type = "source")
library(shiny)
library(mfa)
url <- "https://raw.githubusercontent.com/ucb-stat243/stat243-fall-2016/master/problem-sets/final-project/data/wines.csv"
raw_data <- read.csv(url, header = TRUE, sep = ",")
data <- raw_data[,-1]

# Define UI for application that draws a histogram
ui <- fluidPage(
   # Application title
   titlePanel("Multi-factor Analysis"),
   sidebarLayout(
      sidebarPanel(
        checkboxInput("deffile","Use default dataset (wine dataset)?",TRUE),
        conditionalPanel(condition="!input.deffile",
                         fileInput('file1', 'Choose file to upload',
                                   accept = c(
                                     'text/csv',
                                     'text/comma-separated-values',
                                     'text/tab-separated-values',
                                     'text/plain',
                                     '.csv',
                                     '.tsv'
                                   )
                         )
        ),
        textInput('vec1','Enter a vector (comma delimmited) indicates the column width of each block (subtable)',
                  "6,6,6,5,6,5,4,6,5,4"),
        numericInput('ncomp','The number of component to keep',1,min=1),
        textInput('group','Enter a vector (comma delimmited) indicates the number of observations in  each group',"4,4,4"),
        checkboxInput('Check5','Customize group names?',FALSE),
        conditionalPanel(condition = 'input.Check5',
                        textInput('groupname','Enter a vector(comma delimmited) indicates the names of each group',
                                  "New Zealand,France,Canada")
                        ),
        checkboxInput('Check1','Bar-chart of Eigenvalues',FALSE),
        checkboxInput('Check2','Scatterplot of Factor Scores',FALSE),
        conditionalPanel(condition = 'input.Check2',
                         textInput('dimen','Choose two dimensions to plot(Comma delimmited)')),
        checkboxInput('Check3','Scatterplot of Partial Factors Scores',FALSE),
        conditionalPanel(condition='input.Check3',
                         textInput('dimen1','Choose two dimensions to plot(Comma delimmited)'),
                         textInput('var','Choose which variable(s) (comma delimmited) to display (show all if left empty)')),
        checkboxInput('Check4','Scatterplot of the Loadings',FALSE),
        conditionalPanel(condition='input.Check4',
                         sliderInput('scalex',label=h5('Choose your Scale for loadings for x-axis'),
                                     min=0.1, max =5,value =1, step=0.1),
                         sliderInput('scaley',label=h5('Choose your Scale for loadings for y-axis'),
                                     min=0.1, max =5,value =1, step=0.1),
                         textInput('dimen2','Choose two dimensions to plot(Comma delimmited)'),
                         textInput('var1','Choose which variable(s) (comma delimmited) to display (show all if left empty)'),
                         textInput('loadnames','Enter a vector (comma delimmited) indicates the names of each loading',
                                   dftext)
                         ),
        checkboxInput("Check6","Scatterplot of the Partial Factor Score together with Loadings",FALSE),
        conditionalPanel(condition = 'input.Check6',
                         sliderInput('scale_x',label=h5('Choose your Scale for loadings for x-axis'),
                                     min=0.1, max =5,value =1, step=0.1),
                         sliderInput('scale_y',label=h5('Choose your Scale for loadings for y-axis'),
                                     min=0.1, max =5,value =1, step=0.1),
                         textInput('dimen3','Choose two dimensions to plot(Comma delimmited)'),
                         textInput('var2','Choose which variable(s) (comma delimmited) to display (show all if left empty)'),
                         textInput('loadnames2','Enter a vector (comma delimmited) indicates the names of each loading',
                                   dftext)
                         )
        ),
      mainPanel(
         conditionalPanel(condition='input.Check1',
                          plotOutput('bar')),
         conditionalPanel(condition='input.Check2',
                          plotOutput('scatter1')),
         conditionalPanel(condition='input.Check3',
                          plotOutput('scatter2')),
         conditionalPanel(condition='input.Check4',
                          plotOutput('scatter3')),
         conditionalPanel(condition='input.Check6',
                          plotOutput('scatter4')
                          )
      )
   )
)

# Define server logic required to make plots
server <- function(input, output) {
  #mfa object
  result=reactive({
    if(input$deffile){
      prodata=data
    }
    else{
      infile=input$file1
      prodata=read.csv(infile$datapath,header=TRUE)
    }
    len=as.numeric(unlist(strsplit(input$vec1,',')))
    len=cumsum(len)
    left=c(1,len[-length(len)]+1)
    set=list()
    for (i in seq_along(len)){
      set[[i]]=(left[i]:len[i])
    }
    mfa_gen(data=prodata,sets=set,ncomps=input$ncomp)
  })
  #group name
  gname=reactive({
    if (input$group==""){
      return(NULL)
    }
    else{
      len1=as.numeric(unlist(strsplit(input$group,',')))
      len1=cumsum(len1)
      left1=c(1,len1[-length(len1)]+1)
      setgroup=list()
      for (i in seq_along(len1)){
        setgroup[[i]]=(left1[i]:len1[i])
      }
      if (input$Check5){
        names(setgroup)=unlist(strsplit(input$groupname,','))      
      }
      return(setgroup)
    }
  })
  #Eigenvalue plot
  output$bar=renderPlot({
    res=result()
    ncomps=input$ncomp
    values=round(res$eigenvalues[1:ncomps],3)
    ggplot(data=NULL) +
      geom_linerange(aes(x=1:ncomps,ymin=0,ymax=values))+
      xlab('n Components')+
      ylab('EigenValues')+
      scale_x_continuous(breaks=1:ncomps,label=1:ncomps) +
      geom_text(aes(x=1:ncomps,y=values+0.03,label = values))
  })
  #Scatterplot for factor score
  output$scatter1=renderPlot({
    dimension=as.numeric(unlist(strsplit(input$dimen,',')))
    validate(
      need(length(dimension)==2,'Please input only TWO dimensions')
    )
    plot(mfa = result(), dim1 =dimension[1],
             dim2 =dimension[2], type = 1, cat = gname())
  })
  #Scatterplot for partial factors scores
  output$scatter2=renderPlot({
    dimension=as.numeric(unlist(strsplit(input$dimen1,',')))
    validate(
      need(length(dimension)==2,'Please input only TWO dimensions')
    )
    if (input$var==''){
      plot(mfa=result(), dim1=dimension[1], dim2=dimension[2], 
           type=2, text = NULL, cat = gname(), var =NULL)
    }
    else {
      variable=as.numeric(unlist(strsplit(input$var,',')))
      plot(mfa=result(), dim1=dimension[1], dim2=dimension[2], 
         type=2, text = NULL, cat = gname(), var =variable)
      }
    
  })
  #Scatterplot of the loadings
  output$scatter3=renderPlot({
    dimension=as.numeric(unlist(strsplit(input$dimen2,',')))
    validate(
      need(length(dimension)==2,'Please input only TWO dimensions')
    )
    lnames=NULL
    if (input$loadnames!=""){
      lnames=unlist(strsplit(input$loadnames,","))
    }
    if (input$var1==''){
      plot(mfa=result(), dim1=dimension[1], dim2=dimension[2], 
           type=3, text = lnames, cat = gname(), var =NULL,
           scale_x=input$scalex,scale_y=input$scaley)
    }
    else {
      variable=as.numeric(unlist(strsplit(input$var1,',')))
      plot(mfa=result(), dim1=dimension[1], dim2=dimension[2], 
           type=3, text = lnames, cat = gname(), var =variable,
           scale_x=input$scalex,scale_y=input$scaley)
    }
  })
  
  #scatter plot for PFS+loadings
  output$scatter4=renderPlot({
    dimension=as.numeric(unlist(strsplit(input$dimen3,',')))
    validate(
      need(length(dimension)==2,'Please input only TWO dimensions')
    )
    lnames=NULL
    if (input$loadnames2!=""){
      lnames=unlist(strsplit(input$loadnames2,","))
    }
    if (input$var2==''){
      plot(mfa=result(), dim1=dimension[1], dim2=dimension[2], 
           type=4, text = lnames, cat = gname(), var =NULL,
           scale_x=input$scale_x,scale_y=input$scale_y)
    }
    else {
      variable=as.numeric(unlist(strsplit(input$var2,',')))
      plot(mfa=result(), dim1=dimension[1], dim2=dimension[2], 
           type=4, text = lnames, cat = gname(), var =variable,
           scale_x=input$scale_x,scale_y=input$scale_y)
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

