keysFromSigList = function(sl) {
 snps = names(sl)
 lens = sapply(sl, nrow)
 snps = rep(snps,lens)
 pids = unlist(lapply(sl, function(x) x$"probeset_id"))
 syms = unlist(lapply(sl, function(x) x$"symbol"))
 paste(paste(snps, pids, sep=":"), syms, sep=":")
}


bceBrowse = function() {
#
# globals for quadpanel
#
 library(beeswarm)
 library(bceQTL)
 data(norsig)
 data(ernsig)
 data(erpsig)
 data(allpairedSig)

norkeys = paste("(NOR)", keysFromSigList( norsig ))
ernkeys = paste("(ER-)", keysFromSigList( ernsig ))
erpkeys = paste("(ER+)", keysFromSigList( erpsig ))
paikeys = paste("(PAIR)", keysFromSigList( allpairedSig ))
allkeys = c(norkeys, ernkeys, erpkeys, paikeys)
 
 if (!exists("nhscancpcns_ann")) data(nhscancpcns_ann)
 if (!exists("nhsnormpcns_ann")) data(nhsnormpcns_ann)
 if (!exists("tcann")) data(tcann)
 tumERP = nhscancpcns_ann[, which(nhscancpcns_ann$erpos=="1")]
 tumERN = nhscancpcns_ann[, which(nhscancpcns_ann$erpos!="1")]
 norERP = nhsnormpcns_ann[, which(nhsnormpcns_ann$erpos=="1")]
 norERN = nhsnormpcns_ann[, which(nhsnormpcns_ann$erpos!="1")]

# create combined syncd # tumor-normal

canc = nhscancpcns_ann
sampleNames(canc) = canc$nhsid
norm = nhsnormpcns_ann
sampleNames(norm) = norm$nhsid
common = intersect(sampleNames(canc), sampleNames(norm))
canc = canc[,common]
norm = norm[,common]
canc$isTumor = 1
norm$isTumor = 0
sampleNames(canc) = paste0(sampleNames(canc), "C")
sampleNames(norm) = paste0(sampleNames(norm), "N")
combTumNorm = combine(canc, norm)


layout = function(snp) {
  nt = table(pData(norERP)[[snp]])
  norpdf = data.frame(type = "norm ER+", numB = names(nt), freq=as.numeric(nt))
  nt = table(pData(norERN)[[snp]])
  norndf = data.frame(type = "norm ER-", numB = names(nt), freq=as.numeric(nt))
  nt = table(pData(tumERP)[[snp]])
  tumpdf = data.frame(type = "tum ER+", numB = names(nt), freq=as.numeric(nt))
  nt = table(pData(tumERP)[[snp]])
  tumndf = data.frame(type = "tum ER-", numB = names(nt), freq=as.numeric(nt))
  nt = table(pData(combTumNorm[,which(combTumNorm$isTumor==1)])[[snp]])
  combdf = data.frame(type = "paired", numB = names(nt), freq=as.numeric(nt))
  rbind(norpdf, norndf, tumpdf, tumndf, combdf)
}

layout_r = function(snp) {
  ran = function(x) round(as.numeric(x),0)
  nt = table(round(as.numeric(pData(norERP)[[snp]]),0))
  norpdf = data.frame(type = "norm ER+", numB = names(nt), freq=as.numeric(nt))
  nt = table(ran(pData(norERN)[[snp]]))
  norndf = data.frame(type = "norm ER-", numB = names(nt), freq=as.numeric(nt))
  nt = table(ran(pData(tumERP)[[snp]]))
  tumpdf = data.frame(type = "tum ER+", numB = names(nt), freq=as.numeric(nt))
  nt = table(ran(pData(tumERP)[[snp]]))
  tumndf = data.frame(type = "tum ER-", numB = names(nt), freq=as.numeric(nt))
  nt = table(ran(pData(combTumNorm[,which(combTumNorm$isTumor==1)])[[snp]]))
  combdf = data.frame(type = "paired", numB = names(nt), freq=as.numeric(nt))
  rbind(norpdf, norndf, tumpdf, tumndf, combdf)
}

ui = fluidPage(
 titlePanel("bceQTL browser"),
 sidebarPanel(
  helpText(paste("Select from among", length(allkeys), "SNP-gene pairs, stratified by normal, tumor ER+, tumor ER-, paired tumor:normal.  Available pairs have stratum/SNP-specific FDR < 0.10 in search over 26004 GlueArray transcript clusters using moderated t-statistics (or paired t-statistics for the tumor:normal pairs) in linear models with adjustment for age, plate, and a single expression PC.")),
  fluidRow(
   selectInput("selpair", "SNP-gene pair", choices=allkeys,
      selected=allkeys[1], selectize=TRUE)
   )
 ),
 mainPanel(
  tabsetPanel(
   tabPanel("swarm", plotOutput("swarmov")),
   tabPanel("quad", plotOutput("quadp")),
   tabPanel("stats", helpText("All probes with stratum/SNP-specific FDR < 0.10 are reported."), dataTableOutput("tab")),
   tabPanel("layout", dataTableOutput("lay")),
   tabPanel("layout_r", dataTableOutput("lay_r"))
   )
 )
)
server = function(input, output) {
  output$swarmov = renderPlot({
   strat = strsplit(input$selpair, " ")[[1]][1]
   intmp = strsplit(input$selpair, " ")[[1]][2]
   intmp = strsplit(intmp, ":")[[1]]
   if (strat == "(NOR)") beeswarm(split( exprs(nhsnormpcns_ann[intmp[2],]),
       round(pData(nhsnormpcns_ann)[[intmp[1]]],0) ), xlab=intmp[1], 
       ylab=intmp[2], main="Expression in normal adj.")
   else if (strat == "(ER+)") beeswarm(split( exprs(tumERP[intmp[2],]),
       round(pData(tumERP)[[intmp[1]]],0) ), xlab=intmp[1], 
       ylab=intmp[2], main="Expression in ER+ tumor")
   else if (strat == "(ER-)") beeswarm(split( exprs(tumERN[intmp[2],]),
       round(pData(tumERN)[[intmp[1]]],0) ), xlab=intmp[1], 
       ylab=intmp[2], main="Expression in ER- tumor")
   else if (strat == "(PAIR)") {
          pgt = round(pData(combTumNorm)[[intmp[1]]],0) 
          incol = combTumNorm$isTumor+1
          spl = split(incol, pgt)
          beeswarm(split( exprs(combTumNorm[intmp[2],]), pgt) , xlab=intmp[1], 
               ylab=intmp[2], main="Expression in paired tumor/normal",
               pwcol = split(incol, pgt), pch=19)
          } 
  })
  output$quadp = renderPlot({
   strat = strsplit(input$selpair, " ")[[1]][1]
   if (FALSE) plot(1,1)
   else {
         intmp = strsplit(input$selpair, " ")[[1]][2]
         intmp = strsplit(intmp, ":")[[1]]
         quadpanel(intmp[1], intmp[2], tumERN, tumERP, norERN, norERP)
        }
  })
  output$tab = renderDataTable({
   strat = strsplit(input$selpair, " ")[[1]][1]
   intmp = strsplit(input$selpair, " ")[[1]]
   type = intmp[1]
   intmp = intmp[2]
   intmp = strsplit(intmp, ":")[[1]]
   snp = intmp[1]
   if (type == "(NOR)") cbind(snp=snp, norsig[[ snp ]])
   else if (type == "(ER-)") cbind(snp=snp, ernsig[[ snp ]])
   else if (type == "(ER+)") cbind(snp=snp, erpsig[[ snp ]])
   else if (type == "(PAIR)") cbind(snp=snp, allpairedSig[[ snp ]])
  })
  output$lay = renderDataTable({
   strat = strsplit(input$selpair, " ")[[1]][1]
   intmp = strsplit(input$selpair, " ")[[1]]
   type = intmp[1]
   intmp = intmp[2]
   intmp = strsplit(intmp, ":")[[1]]
   snp = intmp[1]
   layout(snp)
  })
  output$lay_r = renderDataTable({
   strat = strsplit(input$selpair, " ")[[1]][1]
   intmp = strsplit(input$selpair, " ")[[1]]
   type = intmp[1]
   intmp = intmp[2]
   intmp = strsplit(intmp, ":")[[1]]
   snp = intmp[1]
   layout_r(snp)
  })
 }
shinyApp(ui, server)
}
