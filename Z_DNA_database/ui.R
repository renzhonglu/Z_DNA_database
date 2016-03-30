library(shiny)
library(shinythemes)
library(markdown)
library(rmarkdown)
library(DT)

shinyUI(fluidPage(#theme = shinytheme("Readable"),
  titlePanel(h2(em("Z-Catcher2", style = "color: red"), 
                                 "for Predicting Potential Z-DNA Forming Regions in the Whole Genome",align = "center"), windowTitle = "Z-Catcher2"),
  
  
  navbarPage("Content",theme = shinytheme("Readable"),
                   #br(),
                   #br(),
                   tabPanel("Introduction", 
                            mainPanel(width = 10,
                                      h4("Z-DNA Introduction"),
                                      #p("Z-DNA is a "),
                                      includeMarkdown("zdna_intro.md"),
                                      #includeHTML("zdna_intro.html"),
                                      
                                      br(),
                                      img(src = "GraphicalAbstracts.jpg", alt = "abstracts", width = 1280, height = 480),
                                      br(),
                                      h4(em("Z-Cather2"), "Introduction"),
                                      includeMarkdown("zcatcher2_intro.md")
                                      )
                            ),
                   tabPanel("Summary/Download",
                            sidebarLayout(
                              sidebarPanel(width = 3,
                                           helpText("All species genomes were download from UCSC:",br(),
                                                    "Human ",tags$i("(H.sapiens: hg38)"),br(),
                                                    "Chimpanzee ",tags$i("(P.troglodytes: panTro4)"),br(),
                                                    "Mouse ",tags$i("(M.musculus: mm10)"), br(),
                                                    "Rattus",tags$i("(R.norvegicus: rn6)"), br(),
                                                    "Zebrafish", tags$i("(D.rerio: danRer10)"), br(),
                                                    "FruitFly", tags$i("(D.melanogaster: dm6)"), br(),
                                                    "Celegans", tags$i("(C.elegans: ce10)"),br(),
                                                    "Saccerevisiae", tags$i("(S.cerevisiae: sacCer3)"),br(),
                                                    "Arabidopsis", tags$i("(A.thaliana: TAIR10)"), br()
                                           ),
                                           selectInput("species","Select One Species:",choices = NULL),
                                           htmlOutput("selectUI", inline = T),
                                           selectInput("chromosome","Select One Chromosome or ALL:", choices = htmlOutput("varselect")),
                                           helpText("Please select Species firstly, then click Dowload button!"),
                                           downloadButton('downloadData1', 'Download')
                              ),
                              mainPanel(width = 9,
                                        fluidRow(
                                          column(4,
                                                 h4(textOutput("text1")),
                                                 tableOutput("summary")
                                          ),
                                          column(8,
                                                 h4(textOutput("text2")), #h4("Length of ZDRs in selected chromosome"),
                                                 plotOutput("histplot1"),
                                                 h4(textOutput("text3")),
                                                 plotOutput("histplot2")
                                          )
                                        )
                                        #h4("Observations"),
                                        #tableOutput("view")
                              )
                            )
                            ),
             
             
                   tabPanel("Search",
                            sidebarLayout(
                              sidebarPanel(width = 3,
                                fileInput('file1', label = h4('Choose a FASTA File:'),
                                          accept='text'),
                                tags$hr(),
                                sliderInput('slid1', label = h4("Choose the Sigma0 (nsd):"),
                                            min = -0.1,max=-0.035,step = 0.005,-0.07),
                                textInput("text", label = h4("Output file name:"), 
                                          value = "output.txt"),
                                actionButton("submit","Submit")
                              ),
                              mainPanel(width =9,
                                
                                verbatimTextOutput('text'),
                                DT::dataTableOutput('contents'),
                                downloadButton('downloadData2', 'Download')
                              )
                            )
                            
                            ),
                   tabPanel("Help")
),
tags$footer(em("Department of Bioinformatics"),
                  ", Shool of Basic Medical Sciences, Southern Medical University",
                  img(src = "xiaohui.jpg", alt = "xiaohui", width = 64, height = 64),
                  br(),
                  "E-mail:", tags$a(href = "mailto:jmli@smu.edu.cn","jmli@smu.edu.cn"), align = "center")

)
)