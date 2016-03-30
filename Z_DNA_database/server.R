library(shiny)
library(ggplot2)
library(RSQLite)
library(DT)
shinyServer(function(input, output,session) {
  #
  species <- c("Human","Chimpanzee","Mouse","Rattus","Zebrafish","FruitFly","Celegans","Saccerevisiae","Arabidopsis")
  hschronames <- c(paste0("CHR",c(1:22,"X","Y")))
  ptchronames <- c(paste0("CHR",c("1","2A","2B",3:22,"X","Y")))
  mschronames <- c(paste0("CHR",c(1:19,"X","Y")))
  rnchronames <- c(paste0("CHR",c(1:20,"X","Y")))
  drchronames <- c(paste0("CHR",c(1:25)))
  dmchronames <- c(paste0("CHR",c("2L","2R","3L","3R","4","X","Y")))
  cechronames <- c(paste0("CHR",c("I","II","III","IV","V","X")))
  scchronames <- c(paste0("CHR",c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI")))
  atchronames <- c(paste0("CHR",1:5))
  
  
  #
  updateSelectInput(session, inputId = "species", label = "Select One Species:", choices = c("All",species), selected = "All")
  
  #
  output$selectUI <- renderUI({
  })
  observe({
    updateSelectInput(session, inputId = "chromosome", label = "Select One Chromosome or ALL:",
                      choices = chromosomes())
  })
  
  #
  chromosomes <- reactive({switch (input$species, 
                                   "Human"= c("All",hschronames),
                                   "Chimpanzee" = c("All",ptchronames),
                                   "Mouse" = c("All",mschronames),
                                   "Rattus" = c("All", rnchronames),
                                   "Zebrafish" = c("All", drchronames),
                                   "FruitFly" = c("All", dmchronames),
                                   "Celegans" = c("All",cechronames),
                                   "Saccerevisiae" = c("All", scchronames),
                                   "Arabidopsis" = c("All", atchronames),
                                   "Select..."
  )})
  
  #
  con <- RSQLite::dbConnect(SQLite(),"data/ZDRsDatabase.sqlite")
  sqltables <- RSQLite::dbListTables(con)
  #"SaccZDRs", "arabidopsisZDRs", "chimpanzeeZDRs", "flyZDRs", "humanZDRs", "mouseZDRs", "Rattus","wormZDRs", "zebrafishZDRs"
  summarycomm <- paste("SELECT Chromosome, count(Chromosome) AS ZDR_counts from", sqltables, "group by Chromosome")
  summarydata1 <- reactive({switch(input$species,
                                   "Human" = RSQLite::dbGetQuery(con,summarycomm[5]),
                                   "Chimpanzee" = RSQLite::dbGetQuery(con, summarycomm[3]),
                                   "Mouse" = RSQLite::dbGetQuery(con, summarycomm[6]),
                                   "Rattus" = RSQLite::dbGetQuery(con, summarycomm[7]),
                                   "Zebrafish" = RSQLite::dbGetQuery(con, summarycomm[9]),
                                   "FruitFly" = RSQLite::dbGetQuery(con, summarycomm[4]),
                                   "Celegans" = RSQLite::dbGetQuery(con, summarycomm[8]),
                                   "Saccerevisiae" = RSQLite::dbGetQuery(con, summarycomm[1]),
                                   "Arabidopsis" = RSQLite::dbGetQuery(con, summarycomm[2])
  )})
  selectcomm.hs <- paste0("SELECT Chromosome, ZDRlength FROM ", sqltables[5], " Where Chromosome = ", "'",hschronames,"'")
  selectcomm.pt <- paste0("SELECT Chromosome, ZDRlength FROM ", sqltables[3], " Where Chromosome = ", "'",ptchronames,"'")
  selectcomm.mm <- paste0("SELECT Chromosome, ZDRlength FROM ", sqltables[6], " Where Chromosome = ", "'",mschronames,"'")
  selectcomm.rn <- paste0("SELECT Chromosome, ZDRlength FROM ", sqltables[7], " Where Chromosome = ", "'",rnchronames,"'")
  selectcomm.dr <- paste0("SELECT Chromosome, ZDRlength FROM ", sqltables[9], " Where Chromosome = ", "'",drchronames,"'")
  selectcomm.dm <- paste0("SELECT Chromosome, ZDRlength FROM ", sqltables[4], " Where Chromosome = ", "'",dmchronames,"'")
  selectcomm.ce <- paste0("SELECT Chromosome, ZDRlength FROM ", sqltables[8], " Where Chromosome = ", "'",cechronames,"'")
  selectcomm.sc <- paste0("SELECT Chromosome, ZDRlength FROM ", sqltables[1], " Where Chromosome = ", "'",scchronames,"'")
  selectcomm.at <- paste0("SELECT Chromosome, ZDRlength FROM ", sqltables[2], " Where Chromosome = ", "'",atchronames,"'")
  
  
  #
  selectcomm.hs2 <- paste0("SELECT * FROM ", sqltables[5], " Where Chromosome = ", "'",hschronames,"'")
  selectcomm.pt2 <- paste0("SELECT * FROM ", sqltables[3], " Where Chromosome = ", "'",ptchronames,"'")
  selectcomm.mm2 <- paste0("SELECT * FROM ", sqltables[6], " Where Chromosome = ", "'",mschronames,"'")
  selectcomm.rn2 <- paste0("SELECT * FROM ", sqltables[7], " Where Chromosome = ", "'",rnchronames,"'")
  selectcomm.dr2 <- paste0("SELECT * FROM ", sqltables[9], " Where Chromosome = ", "'",drchronames,"'")
  selectcomm.dm2 <- paste0("SELECT * FROM ", sqltables[4], " Where Chromosome = ", "'",dmchronames,"'")
  selectcomm.ce2 <- paste0("SELECT * FROM ", sqltables[8], " Where Chromosome = ", "'",cechronames,"'")
  selectcomm.sc2 <- paste0("SELECT * FROM ", sqltables[1], " Where Chromosome = ", "'",scchronames,"'")
  selectcomm.at2 <- paste0("SELECT * FROM ", sqltables[2], " Where Chromosome = ", "'",atchronames,"'")
  
  ##
  
  #
  totalZDRs <- read.table("data/totalZDR.txt", sep = "\t", header = T)
  
  #
  summarydata2 <- reactive({
    #initial value
    #if(input$species == "All" || input$chromosome == "Select..."){
    #  totalZDRs
    #}
    #human
    if(input$species == "Human") {
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select Chromosome,ZDRlength from ",sqltables[5])),
             'CHR1' = RSQLite::dbGetQuery(con, selectcomm.hs[1]),
             'CHR2' = RSQLite::dbGetQuery(con, selectcomm.hs[2]),
             'CHR3' = RSQLite::dbGetQuery(con, selectcomm.hs[3]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.hs[4]),
             'CHR5' = RSQLite::dbGetQuery(con, selectcomm.hs[5]),
             'CHR6' = RSQLite::dbGetQuery(con, selectcomm.hs[6]),
             'CHR7' = RSQLite::dbGetQuery(con, selectcomm.hs[7]),
             'CHR8' = RSQLite::dbGetQuery(con, selectcomm.hs[8]),
             'CHR9' = RSQLite::dbGetQuery(con, selectcomm.hs[9]),
             'CHR10' = RSQLite::dbGetQuery(con, selectcomm.hs[10]),
             'CHR11' = RSQLite::dbGetQuery(con, selectcomm.hs[11]),
             'CHR12' = RSQLite::dbGetQuery(con, selectcomm.hs[12]),
             'CHR13' = RSQLite::dbGetQuery(con, selectcomm.hs[13]),
             'CHR14' = RSQLite::dbGetQuery(con, selectcomm.hs[14]),
             'CHR15' = RSQLite::dbGetQuery(con, selectcomm.hs[15]),
             'CHR16' = RSQLite::dbGetQuery(con, selectcomm.hs[16]),
             'CHR17' = RSQLite::dbGetQuery(con, selectcomm.hs[17]),
             'CHR18' = RSQLite::dbGetQuery(con, selectcomm.hs[18]),
             'CHR19' = RSQLite::dbGetQuery(con, selectcomm.hs[19]),
             'CHR20' = RSQLite::dbGetQuery(con, selectcomm.hs[20]),
             'CHR21' = RSQLite::dbGetQuery(con, selectcomm.hs[21]),
             'CHR22' = RSQLite::dbGetQuery(con, selectcomm.hs[22]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.hs[23]),
             'CHRY' = RSQLite::dbGetQuery(con, selectcomm.hs[24])
      )
    } else if(input$species == "Chimpanzee") {      #Chimpanzee
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select Chromosome,ZDRlength from ",sqltables[3])),
             'CHR1' = RSQLite::dbGetQuery(con, selectcomm.pt[1]),
             'CHR2A' = RSQLite::dbGetQuery(con, selectcomm.pt[2]),
             'CHR2B' = RSQLite::dbGetQuery(con, selectcomm.pt[3]),
             'CHR3' = RSQLite::dbGetQuery(con, selectcomm.pt[4]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.pt[5]),
             'CHR5' = RSQLite::dbGetQuery(con, selectcomm.pt[6]),
             'CHR6' = RSQLite::dbGetQuery(con, selectcomm.pt[7]),
             'CHR7' = RSQLite::dbGetQuery(con, selectcomm.pt[8]),
             'CHR8' = RSQLite::dbGetQuery(con, selectcomm.pt[9]),
             'CHR9' = RSQLite::dbGetQuery(con, selectcomm.pt[10]),
             'CHR10' = RSQLite::dbGetQuery(con, selectcomm.pt[11]),
             'CHR11' = RSQLite::dbGetQuery(con, selectcomm.pt[12]),
             'CHR12' = RSQLite::dbGetQuery(con, selectcomm.pt[13]),
             'CHR13' = RSQLite::dbGetQuery(con, selectcomm.pt[14]),
             'CHR14' = RSQLite::dbGetQuery(con, selectcomm.pt[15]),
             'CHR15' = RSQLite::dbGetQuery(con, selectcomm.pt[16]),
             'CHR16' = RSQLite::dbGetQuery(con, selectcomm.pt[17]),
             'CHR17' = RSQLite::dbGetQuery(con, selectcomm.pt[18]),
             'CHR18' = RSQLite::dbGetQuery(con, selectcomm.pt[19]),
             'CHR19' = RSQLite::dbGetQuery(con, selectcomm.pt[20]),
             'CHR20' = RSQLite::dbGetQuery(con, selectcomm.pt[21]),
             'CHR21' = RSQLite::dbGetQuery(con, selectcomm.pt[22]),
             'CHR22' = RSQLite::dbGetQuery(con, selectcomm.pt[23]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.pt[24]),
             'CHRY' = RSQLite::dbGetQuery(con, selectcomm.pt[25])
      )
    } else if(input$species == "Mouse") {#Mouse
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select Chromosome,ZDRlength from ",sqltables[6])),
             'CHR1' = RSQLite::dbGetQuery(con, selectcomm.mm[1]),
             'CHR2' = RSQLite::dbGetQuery(con, selectcomm.mm[2]),
             'CHR3' = RSQLite::dbGetQuery(con, selectcomm.mm[3]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.mm[4]),
             'CHR5' = RSQLite::dbGetQuery(con, selectcomm.mm[5]),
             'CHR6' = RSQLite::dbGetQuery(con, selectcomm.mm[6]),
             'CHR7' = RSQLite::dbGetQuery(con, selectcomm.mm[7]),
             'CHR8' = RSQLite::dbGetQuery(con, selectcomm.mm[8]),
             'CHR9' = RSQLite::dbGetQuery(con, selectcomm.mm[9]),
             'CHR10' = RSQLite::dbGetQuery(con, selectcomm.mm[10]),
             'CHR11' = RSQLite::dbGetQuery(con, selectcomm.mm[11]),
             'CHR12' = RSQLite::dbGetQuery(con, selectcomm.mm[12]),
             'CHR13' = RSQLite::dbGetQuery(con, selectcomm.mm[13]),
             'CHR14' = RSQLite::dbGetQuery(con, selectcomm.mm[14]),
             'CHR15' = RSQLite::dbGetQuery(con, selectcomm.mm[15]),
             'CHR16' = RSQLite::dbGetQuery(con, selectcomm.mm[16]),
             'CHR17' = RSQLite::dbGetQuery(con, selectcomm.mm[17]),
             'CHR18' = RSQLite::dbGetQuery(con, selectcomm.mm[18]),
             'CHR19' = RSQLite::dbGetQuery(con, selectcomm.mm[19]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.mm[20]),
             'CHRY' = RSQLite::dbGetQuery(con, selectcomm.mm[21])
      )
    } else if(input$species == "Rattus") {#Mouse
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select Chromosome,ZDRlength from ",sqltables[7])),
             'CHR1' = RSQLite::dbGetQuery(con, selectcomm.rn[1]),
             'CHR2' = RSQLite::dbGetQuery(con, selectcomm.rn[2]),
             'CHR3' = RSQLite::dbGetQuery(con, selectcomm.rn[3]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.rn[4]),
             'CHR5' = RSQLite::dbGetQuery(con, selectcomm.rn[5]),
             'CHR6' = RSQLite::dbGetQuery(con, selectcomm.rn[6]),
             'CHR7' = RSQLite::dbGetQuery(con, selectcomm.rn[7]),
             'CHR8' = RSQLite::dbGetQuery(con, selectcomm.rn[8]),
             'CHR9' = RSQLite::dbGetQuery(con, selectcomm.rn[9]),
             'CHR10' = RSQLite::dbGetQuery(con, selectcomm.rn[10]),
             'CHR11' = RSQLite::dbGetQuery(con, selectcomm.rn[11]),
             'CHR12' = RSQLite::dbGetQuery(con, selectcomm.rn[12]),
             'CHR13' = RSQLite::dbGetQuery(con, selectcomm.rn[13]),
             'CHR14' = RSQLite::dbGetQuery(con, selectcomm.rn[14]),
             'CHR15' = RSQLite::dbGetQuery(con, selectcomm.rn[15]),
             'CHR16' = RSQLite::dbGetQuery(con, selectcomm.rn[16]),
             'CHR17' = RSQLite::dbGetQuery(con, selectcomm.rn[17]),
             'CHR18' = RSQLite::dbGetQuery(con, selectcomm.rn[18]),
             'CHR19' = RSQLite::dbGetQuery(con, selectcomm.rn[19]),
             'CHR20' = RSQLite::dbGetQuery(con, selectcomm.rn[20]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.rn[21]),
             'CHRY' = RSQLite::dbGetQuery(con, selectcomm.rn[22])
      )
    } else if(input$species == "Zebrafish") {#Zebrafish
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select Chromosome,ZDRlength from ",sqltables[9])),
             'CHR1' = RSQLite::dbGetQuery(con, selectcomm.dr[1]),
             'CHR2' = RSQLite::dbGetQuery(con, selectcomm.dr[2]),
             'CHR3' = RSQLite::dbGetQuery(con, selectcomm.dr[3]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.dr[4]),
             'CHR5' = RSQLite::dbGetQuery(con, selectcomm.dr[5]),
             'CHR6' = RSQLite::dbGetQuery(con, selectcomm.dr[6]),
             'CHR7' = RSQLite::dbGetQuery(con, selectcomm.dr[7]),
             'CHR8' = RSQLite::dbGetQuery(con, selectcomm.dr[8]),
             'CHR9' = RSQLite::dbGetQuery(con, selectcomm.dr[9]),
             'CHR10' = RSQLite::dbGetQuery(con, selectcomm.dr[10]),
             'CHR11' = RSQLite::dbGetQuery(con, selectcomm.dr[11]),
             'CHR12' = RSQLite::dbGetQuery(con, selectcomm.dr[12]),
             'CHR13' = RSQLite::dbGetQuery(con, selectcomm.dr[13]),
             'CHR14' = RSQLite::dbGetQuery(con, selectcomm.dr[14]),
             'CHR15' = RSQLite::dbGetQuery(con, selectcomm.dr[15]),
             'CHR16' = RSQLite::dbGetQuery(con, selectcomm.dr[16]),
             'CHR17' = RSQLite::dbGetQuery(con, selectcomm.dr[17]),
             'CHR18' = RSQLite::dbGetQuery(con, selectcomm.dr[18]),
             'CHR19' = RSQLite::dbGetQuery(con, selectcomm.dr[19]),
             'CHR20' = RSQLite::dbGetQuery(con, selectcomm.dr[20]),
             'CHR21' = RSQLite::dbGetQuery(con, selectcomm.dr[21]),
             'CHR22' = RSQLite::dbGetQuery(con, selectcomm.dr[22]),
             'CHR23' = RSQLite::dbGetQuery(con, selectcomm.dr[23]),
             'CHR24' = RSQLite::dbGetQuery(con, selectcomm.dr[24]),
             'CHR25' = RSQLite::dbGetQuery(con, selectcomm.dr[25])
      )
    } else if(input$species == "FruitFly") { #FruitFly
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select Chromosome,ZDRlength from ",sqltables[4])),
             'CHR2L' = RSQLite::dbGetQuery(con, selectcomm.dm[1]),
             'CHR2R' = RSQLite::dbGetQuery(con, selectcomm.dm[2]),
             'CHR3L' = RSQLite::dbGetQuery(con, selectcomm.dm[3]),
             'CHR3R' = RSQLite::dbGetQuery(con, selectcomm.dm[4]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.dm[5]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.dm[6]),
             'CHRY' = RSQLite::dbGetQuery(con, selectcomm.dm[7])
      )
    } else if(input$species == "Celegans") {#Celegans
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select Chromosome,ZDRlength from ",sqltables[8])),
             'CHRI' = RSQLite::dbGetQuery(con, selectcomm.ce[1]),
             'CHRII' = RSQLite::dbGetQuery(con, selectcomm.ce[2]),
             'CHRIII' = RSQLite::dbGetQuery(con, selectcomm.ce[3]),
             'CHRIV' = RSQLite::dbGetQuery(con, selectcomm.ce[4]),
             'CHRV' = RSQLite::dbGetQuery(con, selectcomm.ce[5]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.ce[6])
      )
    } else if(input$species == "Saccerevisiae") {#Saccerevisiae
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select Chromosome,ZDRlength from ",sqltables[1])),
             'CHRI' = RSQLite::dbGetQuery(con, selectcomm.sc[1]),
             'CHRII' = RSQLite::dbGetQuery(con, selectcomm.sc[2]),
             'CHRIII' = RSQLite::dbGetQuery(con, selectcomm.sc[3]),
             'CHRIV' = RSQLite::dbGetQuery(con, selectcomm.sc[4]),
             'CHRV' = RSQLite::dbGetQuery(con, selectcomm.sc[5]),
             'CHRVI' = RSQLite::dbGetQuery(con, selectcomm.sc[6]),
             'CHRVII' = RSQLite::dbGetQuery(con, selectcomm.sc[7]),
             'CHRVIII' = RSQLite::dbGetQuery(con, selectcomm.sc[8]),
             'CHRIX' = RSQLite::dbGetQuery(con, selectcomm.sc[9]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.sc[10]),
             'CHRXI' = RSQLite::dbGetQuery(con, selectcomm.sc[11]),
             'CHRXII' = RSQLite::dbGetQuery(con, selectcomm.sc[12]),
             'CHRXIII' = RSQLite::dbGetQuery(con, selectcomm.sc[13]),
             'CHRXIV' = RSQLite::dbGetQuery(con, selectcomm.sc[14]),
             'CHRXV' = RSQLite::dbGetQuery(con, selectcomm.sc[15]),
             'CHRXVI' = RSQLite::dbGetQuery(con, selectcomm.sc[16])
      )
    } else if(input$species == "Arabidopsis") {#Arabidopsis
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select Chromosome,ZDRlength from ",sqltables[2])),
             'CHR1' = RSQLite::dbGetQuery(con, selectcomm.at[1]),
             'CHR2' = RSQLite::dbGetQuery(con, selectcomm.at[2]),
             'CHR3' = RSQLite::dbGetQuery(con, selectcomm.at[3]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.at[4]),
             'CHR5' = RSQLite::dbGetQuery(con, selectcomm.at[5])
             
      )
    }
    
    
  })
  ####
  
  ###_____________________________________________________________________________
  data2 <- reactive({
    #initial value
    #if(input$species == "All" || input$chromosome == "Select..."){
    #  totalZDRs
    #}
    #human
    if(input$species == "Human") {
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select * from ",sqltables[5])),
             'CHR1' = RSQLite::dbGetQuery(con, selectcomm.hs2[1]),
             'CHR2' = RSQLite::dbGetQuery(con, selectcomm.hs2[2]),
             'CHR3' = RSQLite::dbGetQuery(con, selectcomm.hs2[3]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.hs2[4]),
             'CHR5' = RSQLite::dbGetQuery(con, selectcomm.hs2[5]),
             'CHR6' = RSQLite::dbGetQuery(con, selectcomm.hs2[6]),
             'CHR7' = RSQLite::dbGetQuery(con, selectcomm.hs2[7]),
             'CHR8' = RSQLite::dbGetQuery(con, selectcomm.hs2[8]),
             'CHR9' = RSQLite::dbGetQuery(con, selectcomm.hs2[9]),
             'CHR10' = RSQLite::dbGetQuery(con, selectcomm.hs2[10]),
             'CHR11' = RSQLite::dbGetQuery(con, selectcomm.hs2[11]),
             'CHR12' = RSQLite::dbGetQuery(con, selectcomm.hs2[12]),
             'CHR13' = RSQLite::dbGetQuery(con, selectcomm.hs2[13]),
             'CHR14' = RSQLite::dbGetQuery(con, selectcomm.hs2[14]),
             'CHR15' = RSQLite::dbGetQuery(con, selectcomm.hs2[15]),
             'CHR16' = RSQLite::dbGetQuery(con, selectcomm.hs2[16]),
             'CHR17' = RSQLite::dbGetQuery(con, selectcomm.hs2[17]),
             'CHR18' = RSQLite::dbGetQuery(con, selectcomm.hs2[18]),
             'CHR19' = RSQLite::dbGetQuery(con, selectcomm.hs2[19]),
             'CHR20' = RSQLite::dbGetQuery(con, selectcomm.hs2[20]),
             'CHR21' = RSQLite::dbGetQuery(con, selectcomm.hs2[21]),
             'CHR22' = RSQLite::dbGetQuery(con, selectcomm.hs2[22]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.hs2[23]),
             'CHRY' = RSQLite::dbGetQuery(con, selectcomm.hs2[24])
      )
    } else if(input$species == "Chimpanzee") {      #Chimpanzee
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select * from ",sqltables[3])),
             'CHR1' = RSQLite::dbGetQuery(con, selectcomm.pt2[1]),
             'CHR2A' = RSQLite::dbGetQuery(con, selectcomm.pt2[2]),
             'CHR2B' = RSQLite::dbGetQuery(con, selectcomm.pt2[3]),
             'CHR3' = RSQLite::dbGetQuery(con, selectcomm.pt2[4]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.pt2[5]),
             'CHR5' = RSQLite::dbGetQuery(con, selectcomm.pt2[6]),
             'CHR6' = RSQLite::dbGetQuery(con, selectcomm.pt2[7]),
             'CHR7' = RSQLite::dbGetQuery(con, selectcomm.pt2[8]),
             'CHR8' = RSQLite::dbGetQuery(con, selectcomm.pt2[9]),
             'CHR9' = RSQLite::dbGetQuery(con, selectcomm.pt2[10]),
             'CHR10' = RSQLite::dbGetQuery(con, selectcomm.pt2[11]),
             'CHR11' = RSQLite::dbGetQuery(con, selectcomm.pt2[12]),
             'CHR12' = RSQLite::dbGetQuery(con, selectcomm.pt2[13]),
             'CHR13' = RSQLite::dbGetQuery(con, selectcomm.pt2[14]),
             'CHR14' = RSQLite::dbGetQuery(con, selectcomm.pt2[15]),
             'CHR15' = RSQLite::dbGetQuery(con, selectcomm.pt2[16]),
             'CHR16' = RSQLite::dbGetQuery(con, selectcomm.pt2[17]),
             'CHR17' = RSQLite::dbGetQuery(con, selectcomm.pt2[18]),
             'CHR18' = RSQLite::dbGetQuery(con, selectcomm.pt2[19]),
             'CHR19' = RSQLite::dbGetQuery(con, selectcomm.pt2[20]),
             'CHR20' = RSQLite::dbGetQuery(con, selectcomm.pt2[21]),
             'CHR21' = RSQLite::dbGetQuery(con, selectcomm.pt2[22]),
             'CHR22' = RSQLite::dbGetQuery(con, selectcomm.pt2[23]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.pt2[24]),
             'CHRY' = RSQLite::dbGetQuery(con, selectcomm.pt2[25])
      )
    } else if(input$species == "Mouse") {#Mouse
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select * from ",sqltables[6])),
             'CHR1' = RSQLite::dbGetQuery(con, selectcomm.mm2[1]),
             'CHR2' = RSQLite::dbGetQuery(con, selectcomm.mm2[2]),
             'CHR3' = RSQLite::dbGetQuery(con, selectcomm.mm2[3]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.mm2[4]),
             'CHR5' = RSQLite::dbGetQuery(con, selectcomm.mm2[5]),
             'CHR6' = RSQLite::dbGetQuery(con, selectcomm.mm2[6]),
             'CHR7' = RSQLite::dbGetQuery(con, selectcomm.mm2[7]),
             'CHR8' = RSQLite::dbGetQuery(con, selectcomm.mm2[8]),
             'CHR9' = RSQLite::dbGetQuery(con, selectcomm.mm2[9]),
             'CHR10' = RSQLite::dbGetQuery(con, selectcomm.mm2[10]),
             'CHR11' = RSQLite::dbGetQuery(con, selectcomm.mm2[11]),
             'CHR12' = RSQLite::dbGetQuery(con, selectcomm.mm2[12]),
             'CHR13' = RSQLite::dbGetQuery(con, selectcomm.mm2[13]),
             'CHR14' = RSQLite::dbGetQuery(con, selectcomm.mm2[14]),
             'CHR15' = RSQLite::dbGetQuery(con, selectcomm.mm2[15]),
             'CHR16' = RSQLite::dbGetQuery(con, selectcomm.mm2[16]),
             'CHR17' = RSQLite::dbGetQuery(con, selectcomm.mm2[17]),
             'CHR18' = RSQLite::dbGetQuery(con, selectcomm.mm2[18]),
             'CHR19' = RSQLite::dbGetQuery(con, selectcomm.mm2[19]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.mm2[20]),
             'CHRY' = RSQLite::dbGetQuery(con, selectcomm.mm2[21])
      )
    } else if(input$species == "Rattus") {#Mouse
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select * from ",sqltables[7])),
             'CHR1' = RSQLite::dbGetQuery(con, selectcomm.rn2[1]),
             'CHR2' = RSQLite::dbGetQuery(con, selectcomm.rn2[2]),
             'CHR3' = RSQLite::dbGetQuery(con, selectcomm.rn2[3]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.rn2[4]),
             'CHR5' = RSQLite::dbGetQuery(con, selectcomm.rn2[5]),
             'CHR6' = RSQLite::dbGetQuery(con, selectcomm.rn2[6]),
             'CHR7' = RSQLite::dbGetQuery(con, selectcomm.rn2[7]),
             'CHR8' = RSQLite::dbGetQuery(con, selectcomm.rn2[8]),
             'CHR9' = RSQLite::dbGetQuery(con, selectcomm.rn2[9]),
             'CHR10' = RSQLite::dbGetQuery(con, selectcomm.rn2[10]),
             'CHR11' = RSQLite::dbGetQuery(con, selectcomm.rn2[11]),
             'CHR12' = RSQLite::dbGetQuery(con, selectcomm.rn2[12]),
             'CHR13' = RSQLite::dbGetQuery(con, selectcomm.rn2[13]),
             'CHR14' = RSQLite::dbGetQuery(con, selectcomm.rn2[14]),
             'CHR15' = RSQLite::dbGetQuery(con, selectcomm.rn2[15]),
             'CHR16' = RSQLite::dbGetQuery(con, selectcomm.rn2[16]),
             'CHR17' = RSQLite::dbGetQuery(con, selectcomm.rn2[17]),
             'CHR18' = RSQLite::dbGetQuery(con, selectcomm.rn2[18]),
             'CHR19' = RSQLite::dbGetQuery(con, selectcomm.rn2[19]),
             'CHR20' = RSQLite::dbGetQuery(con, selectcomm.rn2[20]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.rn2[21]),
             'CHRY' = RSQLite::dbGetQuery(con, selectcomm.rn2[22])
      )
    } else if(input$species == "Zebrafish") {#Zebrafish
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select * from ",sqltables[9])),
             'CHR1' = RSQLite::dbGetQuery(con, selectcomm.dr2[1]),
             'CHR2' = RSQLite::dbGetQuery(con, selectcomm.dr2[2]),
             'CHR3' = RSQLite::dbGetQuery(con, selectcomm.dr2[3]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.dr2[4]),
             'CHR5' = RSQLite::dbGetQuery(con, selectcomm.dr2[5]),
             'CHR6' = RSQLite::dbGetQuery(con, selectcomm.dr2[6]),
             'CHR7' = RSQLite::dbGetQuery(con, selectcomm.dr2[7]),
             'CHR8' = RSQLite::dbGetQuery(con, selectcomm.dr2[8]),
             'CHR9' = RSQLite::dbGetQuery(con, selectcomm.dr2[9]),
             'CHR10' = RSQLite::dbGetQuery(con, selectcomm.dr2[10]),
             'CHR11' = RSQLite::dbGetQuery(con, selectcomm.dr2[11]),
             'CHR12' = RSQLite::dbGetQuery(con, selectcomm.dr2[12]),
             'CHR13' = RSQLite::dbGetQuery(con, selectcomm.dr2[13]),
             'CHR14' = RSQLite::dbGetQuery(con, selectcomm.dr2[14]),
             'CHR15' = RSQLite::dbGetQuery(con, selectcomm.dr2[15]),
             'CHR16' = RSQLite::dbGetQuery(con, selectcomm.dr2[16]),
             'CHR17' = RSQLite::dbGetQuery(con, selectcomm.dr2[17]),
             'CHR18' = RSQLite::dbGetQuery(con, selectcomm.dr2[18]),
             'CHR19' = RSQLite::dbGetQuery(con, selectcomm.dr2[19]),
             'CHR20' = RSQLite::dbGetQuery(con, selectcomm.dr2[20]),
             'CHR21' = RSQLite::dbGetQuery(con, selectcomm.dr2[21]),
             'CHR22' = RSQLite::dbGetQuery(con, selectcomm.dr2[22]),
             'CHR23' = RSQLite::dbGetQuery(con, selectcomm.dr2[23]),
             'CHR24' = RSQLite::dbGetQuery(con, selectcomm.dr2[24]),
             'CHR25' = RSQLite::dbGetQuery(con, selectcomm.dr2[25])
      )
    } else if(input$species == "FruitFly") { #FruitFly
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select * from ",sqltables[4])),
             'CHR2L' = RSQLite::dbGetQuery(con, selectcomm.dm2[1]),
             'CHR2R' = RSQLite::dbGetQuery(con, selectcomm.dm2[2]),
             'CHR3L' = RSQLite::dbGetQuery(con, selectcomm.dm2[3]),
             'CHR3R' = RSQLite::dbGetQuery(con, selectcomm.dm2[4]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.dm2[5]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.dm2[6]),
             'CHRY' = RSQLite::dbGetQuery(con, selectcomm.dm2[7])
      )
    } else if(input$species == "Celegans") {#Celegans
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select * from ",sqltables[8])),
             'CHRI' = RSQLite::dbGetQuery(con, selectcomm.ce2[1]),
             'CHRII' = RSQLite::dbGetQuery(con, selectcomm.ce2[2]),
             'CHRIII' = RSQLite::dbGetQuery(con, selectcomm.ce2[3]),
             'CHRIV' = RSQLite::dbGetQuery(con, selectcomm.ce2[4]),
             'CHRV' = RSQLite::dbGetQuery(con, selectcomm.ce2[5]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.ce2[6])
      )
    } else if(input$species == "Saccerevisiae") {#Saccerevisiae
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select * from ",sqltables[1])),
             'CHRI' = RSQLite::dbGetQuery(con, selectcomm.sc2[1]),
             'CHRII' = RSQLite::dbGetQuery(con, selectcomm.sc2[2]),
             'CHRIII' = RSQLite::dbGetQuery(con, selectcomm.sc2[3]),
             'CHRIV' = RSQLite::dbGetQuery(con, selectcomm.sc2[4]),
             'CHRV' = RSQLite::dbGetQuery(con, selectcomm.sc2[5]),
             'CHRVI' = RSQLite::dbGetQuery(con, selectcomm.sc2[6]),
             'CHRVII' = RSQLite::dbGetQuery(con, selectcomm.sc2[7]),
             'CHRVIII' = RSQLite::dbGetQuery(con, selectcomm.sc2[8]),
             'CHRIX' = RSQLite::dbGetQuery(con, selectcomm.sc2[9]),
             'CHRX' = RSQLite::dbGetQuery(con, selectcomm.sc2[10]),
             'CHRXI' = RSQLite::dbGetQuery(con, selectcomm.sc2[11]),
             'CHRXII' = RSQLite::dbGetQuery(con, selectcomm.sc2[12]),
             'CHRXIII' = RSQLite::dbGetQuery(con, selectcomm.sc2[13]),
             'CHRXIV' = RSQLite::dbGetQuery(con, selectcomm.sc2[14]),
             'CHRXV' = RSQLite::dbGetQuery(con, selectcomm.sc2[15]),
             'CHRXVI' = RSQLite::dbGetQuery(con, selectcomm.sc2[16])
      )
    } else if(input$species == "Arabidopsis") {#Arabidopsis
      switch(input$chromosome, 
             'All' = RSQLite::dbGetQuery(con, paste("select * from ",sqltables[2])),
             'CHR1' = RSQLite::dbGetQuery(con, selectcomm.at2[1]),
             'CHR2' = RSQLite::dbGetQuery(con, selectcomm.at2[2]),
             'CHR3' = RSQLite::dbGetQuery(con, selectcomm.at2[3]),
             'CHR4' = RSQLite::dbGetQuery(con, selectcomm.at2[4]),
             'CHR5' = RSQLite::dbGetQuery(con, selectcomm.at2[5])
             
      )
    }
    
    
  })
  #___________________________________###_____________________###________________
  
  
  
  
  output$text1 <- renderText({
    if(input$species == "All" || input$chromosome == "Select..."){
      "Summary of ZDRs in total speices:"
    }else {
      "Summary of selected species:"
    }
    
  })
  #
  
  output$summary <- renderTable({
    if(input$species == "All" || input$chromosome == "Select..."){
      data.frame(totalZDRs)
    }else {
      data.frame(summarydata1())
    }
    
  })
  #
  
  output$text2 <- renderText({
    if(input$species == "All" || input$chromosome == "Select..."){
      "Number of ZDRs in each speices:"
    }else {
      "Length of ZDRs in selected chromosome:"
    }
    
  })
  output$text3 <- renderText({
    if(input$species == "All" || input$chromosome == "Select..."){
      
    } else {
      "ZDRs Counts in Chromosomes:"}
  })
  #
  
  output$histplot1 <- renderPlot({
    if(input$chromosome == "Select..." || input$species == "All"){
      #barplot(totalZDRs$Total_, names.arg = totalZDRs$Species)
      ggplot(totalZDRs,aes(x= Species,y = TotalZDRs ,fill = TotalZDRs)) + geom_bar(stat="identity") + ggtitle("Total ZDRs in each Species") + ylab("Number of ZDRs") + scale_x_discrete(limits=totalZDRs$Species) + theme(axis.title.x= element_text(size=12),axis.title.y=element_text(size=12),plot.title = element_text(size=20),legend.position = "none")
    } else {
      ggplot(summarydata2(), aes(ZDRlength, fill = ..count..)) + geom_histogram(position = "identity",bins = 30) + theme(axis.title.x= element_text(size=12),axis.title.y=element_text(size=12),legend.position = "none") + ylab("Frequency")
    }
  })
  #
  
  output$histplot2 <- renderPlot({
    if(input$chromosome == "Select..." || input$species == "All") {
      return(NULL)
    } else {
      ggplot(summarydata1(),aes(x = Chromosome, y = ZDR_counts, fill = ZDR_counts)) + geom_bar(stat="identity")+ xlab(input$species) +
        # scale_x_discrete(limits= chromosomes[-1]) +
        theme(legend.position="none",axis.text.x=element_text(colour= "black",size = 10,angle=45,vjust=0.5), axis.title.x= element_text(size=12),axis.title.y=element_text(size=12))
    }
  })
  #
  
  output$downloadData1 <- downloadHandler(
    filename = function() { 
      paste(input$species, input$chromosome, '.csv', sep='_') 
    },
    content = function(file) {
      write.csv(data2(), file)
    }
  )
  ################Search page
  intr <- c("This program is used to predict the potential Z-DNA regions in given Fasta sequences.\n",
            "The program needs three arguments:\n",
            "1. Input file:  the file must be FASTA format of sequences without n or N base and gap.\n",
            "2. Sigma0 value:  it is negative supercoiling density always negative number.\n",
            "3. Output file name: you must assign a name for output file with ZDRs results. The default output filename is Sigma-0.07_out.txt\n",
            "Then get it!\n")
  
  #____________________________
  output$text <- renderPrint({
    #     
    inFile <- isolate(input$file1)
    if (is.null(inFile) | input$submit[1] == 0)
      return(cat(intr))
    ofn <- isolate(input$text)
    sigma <- isolate(input$slid1)
    if(sigma <= -0.07){
      a<-system2('perl', c("Z_catcher2.pl","-i",inFile$datapath, "-o",ofn,"-s",sigma),stdout = T)
      cat(a,sep = "\n")
      cat("Input file is:",inFile$name,"\n")
      cat("Output file is: ",ofn,sep = "","\n")
    } else{
      a<-system2('perl', c("Z_catcher2_lowenergy.pl","-i",inFile$datapath, "-o",ofn,"-s",sigma),stdout = T)
      cat(a,sep = "\n")
      cat("Input file is:",inFile$name,"\n")
      cat("Output file is: ",ofn,sep = "","\n")
    }
    #input$submit[1] =0
  })
  #   
  tmp <- reactive({
    input$submit
    isolate({read.delim(input$text, header = T,sep="\t")})
  })
  #
  output$contents <- DT::renderDataTable({
    #     
    if(input$submit == 0){
      return()
    }
    input$submit
    
    tmp()
  })
  #   
  

  #
  output$downloadData2 <- downloadHandler(
    
    filename = isolate(input$text),
    content = function(file) {
      write.csv(tmp(),file,quote=F,row.names = F)
    }
  )
  
  
  
  
})
  
  
  
  